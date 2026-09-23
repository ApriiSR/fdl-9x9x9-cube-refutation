#!/usr/bin/env python3
"""catalogue.py -- assemble, audit and validate.

The catalogue is the concatenation of the per-shard payload files in increasing
shard index.  Since `enum` sorts the records inside a shard, and shard index is
the lexicographic rank of the shard's row 0, the catalogue is exactly the list
of all supports of [n]^3 in lexicographic order: a deterministic function of the
mathematics, with no dependence on the search order.  Its SHA-256 is therefore
a meaningful thing to publish and to compare across independent runs.

A matching SHA-256 is agreement with reference bytes, not a proof of
exhaustion, so the audit below is written to establish exhaustion separately:
given the shard universe, it requires the manifest to cover it exactly, every
payload to be present at its recorded length, every record to carry its own
shard's row 0, and the records to be in strict lexicographic order within each
shard and the shards in strict row order between them.

Subcommands (each has --help):

  pack N MANIFEST SHARDDIR OUT.bin --shards FILE
                                     assemble and checksum the catalogue
  audit N MANIFEST SHARDDIR --shards FILE
                                     full audit against the expected universe
  ledger IN.jsonl OUT.jsonl          drop timing fields, leaving a reproducible
                                     per-query result file, and summarize it
  validate N IN.jsonl --queries K    the strict result validator: what the
                                     ledger must SAY for the theorem to follow
  orbitstats IN.jsonl                orbit-size histogram and total
  sample N ORBITS K OUT.txt          K mapped shards, drawn uniformly, in the
                                     format `shards list` writes
  setcmpshards N ORBITS A B          compare two shard trees on the sampled
                                     shards, as sets and byte for byte
  reports KIND DIR ...               require one complete report per worker
  checksums WANT GOT --require ...   compare two checksum manifests strictly
"""
import argparse
import collections
import glob
import hashlib
import json
import math
import os
import re
import statistics
import sys


# ------------------------------------------------------------------ manifests
HEX16 = re.compile(r'\A[0-9a-f]{16}\Z')


class ManifestError(Exception):
    pass


def _require(cond, lineno, msg):
    if not cond:
        raise ManifestError(f'line {lineno}: {msg}')


def _check_record(r, lineno, n):
    """One completion record, against the single schema the writers produce."""
    w = n * n
    _require(isinstance(r, dict), lineno, 'not a JSON object')
    _require(isinstance(r.get('idx'), int) and r['idx'] >= 0, lineno, 'bad or missing idx')
    _require(isinstance(r.get('status'), str), lineno, 'bad or missing status')
    _require(isinstance(r.get('done'), bool), lineno, 'bad or missing done')
    if not r['done']:
        _require(r['status'] in ('BUDGET', 'BADROW0'), lineno,
                 f'unfinished shard with status {r["status"]!r}')
        return
    _require(r['status'] in ('EXHAUSTED', 'MAPPED'), lineno,
             f'finished shard with status {r["status"]!r}')
    row = r.get('row0')
    _require(isinstance(row, list) and len(row) == n
             and all(isinstance(v, int) and 0 <= v < n for v in row), lineno, 'bad row0')
    _require(sorted(row) == list(range(n)), lineno, 'row0 is not a permutation')
    _require(isinstance(r.get('count'), int) and r['count'] >= 0, lineno, 'bad count')
    _require(r.get('bytes') == r['count'] * w, lineno, 'bytes does not match count')
    _require(isinstance(r.get('digest'), str) and HEX16.match(r['digest']), lineno,
             'missing or malformed payload digest')
    _require(isinstance(r.get('nodes'), int) and r['nodes'] >= 0, lineno, 'bad nodes')
    _require(isinstance(r.get('wall'), (int, float)), lineno, 'bad wall')


def read_manifest(path, n):
    """{idx: record} for the shards marked done, plus the unfinished records.

    One schema, applied to every line.  A malformed line anywhere in the file is
    an error: silently dropping it is how a run comes to look complete when it
    is not.  The single exception is a final line with no newline, which is the
    signature of a kill during a write and costs at most the shard that was in
    flight; it is discarded with a warning.

    Two lines for the same shard are allowed only if they agree about what the
    shard IS -- its row, its record count, its length and its digest.  A worker
    legitimately redoes a shard whose payload went missing, and the second
    attempt will have taken a different length of time and a different number of
    search nodes; but two lines claiming different CONTENTS for the same shard
    are a real inconsistency and an error.
    """
    done, other = {}, []
    with open(path, 'rb') as f:
        raw = f.read()
    text = raw.decode('utf-8', errors='replace')
    torn = None
    if text and not text.endswith('\n'):
        cut = text.rfind('\n') + 1
        torn = text[cut:]
        text = text[:cut]
    for lineno, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except ValueError as e:
            raise ManifestError(f'line {lineno}: not JSON ({e})')
        _check_record(r, lineno, n)
        if not r['done']:
            other.append(r)
            continue
        prev = done.get(r['idx'])
        if prev is not None:
            key = ('row0', 'count', 'bytes', 'digest')
            if [prev[k] for k in key] != [r[k] for k in key]:
                raise ManifestError(f'line {lineno}: shard {r["idx"]} recorded twice, '
                                    f'with different contents')
        done[r['idx']] = r
    if torn is not None:
        print(f'{path}: discarding a torn final line of {len(torn)} bytes',
              file=sys.stderr)
    return done, other


def read_shardfile(path, n):
    """[(idx, row0)] from a `shards list` file, checked to be the shard universe.

    Every row must be an admissible permutation -- exactly one fixed and one
    reflected point -- and the indices must be exactly 0..K-1 in lexicographic
    row order, which is what makes a shard index a canonical name for a shard
    and the catalogue's order a property of the mathematics.
    """
    rows = {}
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != n + 1:
                raise ManifestError(f'{path} line {lineno}: expected {n + 1} fields')
            idx, row = int(parts[0]), [int(v) for v in parts[1:]]
            if sorted(row) != list(range(n)):
                raise ManifestError(f'{path} line {lineno}: row is not a permutation')
            f_pts = sum(1 for j, v in enumerate(row) if v == j)
            r_pts = sum(1 for j, v in enumerate(row) if v == n - 1 - j)
            if f_pts != 1 or r_pts != 1:
                raise ManifestError(f'{path} line {lineno}: shard {idx} row is not '
                                    f'admissible ({f_pts} fixed, {r_pts} reflected)')
            if idx in rows:
                raise ManifestError(f'{path} line {lineno}: shard {idx} repeats')
            rows[idx] = row
    ordered = sorted(rows)
    if ordered != list(range(len(rows))):
        raise ManifestError(f'{path}: shard indices are not 0..{len(rows) - 1}')
    for a, b in zip(ordered, ordered[1:]):
        if rows[a] >= rows[b]:
            raise ManifestError(f'{path}: shard {b} does not follow shard {a} in '
                                f'lexicographic row order')
    return rows


def shard_path(base, idx):
    return os.path.join(base, f'{idx // 1000:03d}', f's{idx:05d}.bin')


# ----------------------------------------------------------------------- pack
def cmd_pack(a):
    n = a.n
    w = n * n
    man, _ = read_manifest(a.manifest, n)
    if a.shards:
        want = read_shardfile(a.shards, n)
        missing = sorted(set(want) - set(man))
        extra = sorted(set(man) - set(want))
        if missing or extra:
            sys.exit(f'refusing to assemble: {len(missing)} shards of the universe are '
                     f'not in the manifest {missing[:10]}, and {len(extra)} manifest '
                     f'shards are not in the universe {extra[:10]}')
    idxs = sorted(man)
    h = hashlib.sha256()
    total = 0
    with open(a.out, 'wb') as out:
        for i in idxs:
            raw = open(shard_path(a.sharddir, i), 'rb').read()
            if len(raw) != man[i]['count'] * w:
                sys.exit(f'shard {i}: {len(raw)} bytes, manifest claims {man[i]["count"]} records')
            out.write(raw)
            h.update(raw)
            total += man[i]['count']
    meta = dict(n=n, shards=len(idxs), supports=total, bytes=total * w,
                sha256=h.hexdigest(), record_bytes=w,
                universe=a.shards or None,
                order='shard index, records sorted within a shard')
    with open(a.out + '.meta.json', 'w') as f:
        json.dump(meta, f, indent=1)
        f.write('\n')
    print(json.dumps(meta))
    return 0


# ---------------------------------------------------------------------- audit
def cmd_audit(a):
    """Manifest against disk, and both against the expected shard universe.

    What this establishes beyond "the bytes hash to what they hashed to before":

      * exhaustion -- the manifest covers the universe exactly, no shard missing
        and none invented;
      * that each shard is the shard it claims to be -- its manifest row and
        every record's own row 0 agree with the universe's row for that index;
      * canonical order -- records strictly increasing inside a shard (which
        also rules out duplicates), and shard index order agreeing with
        lexicographic row order, so the concatenation really is the supports in
        lexicographic order;
      * that every byte of every record is a coordinate in [n].
    """
    n = a.n
    w = n * n
    man, unfinished = read_manifest(a.manifest, n)
    problems = []
    universe = None
    if a.shards:
        universe = read_shardfile(a.shards, n)
        for i in sorted(set(universe) - set(man)):
            problems.append(f'shard {i}: in the universe, not finished in the manifest')
        for i in sorted(set(man) - set(universe)):
            problems.append(f'shard {i}: finished in the manifest, not in the universe')
    elif not a.partial:
        sys.exit('audit needs --shards FILE (the expected universe), or --partial '
                 'to say that an incomplete audit is intended')

    rows = set()
    dup_records = 0
    wrong_row0 = 0
    unsorted_shards = 0
    out_of_range = 0
    walls, nodes, counts = [], [], []
    for i, r in sorted(man.items()):
        p = shard_path(a.sharddir, i)
        if not os.path.exists(p):
            problems.append(f'shard {i}: payload missing')
            continue
        size = os.path.getsize(p)
        if size != r['count'] * w:
            problems.append(f'shard {i}: {size} bytes, manifest claims {r["count"]} records')
            continue
        raw = open(p, 'rb').read()
        recs = [raw[k * w:(k + 1) * w] for k in range(r['count'])]
        if any(b >= n for b in raw):
            out_of_range += 1
            problems.append(f'shard {i}: a record byte is not a coordinate in [{n}]')
        if any(recs[k] >= recs[k + 1] for k in range(len(recs) - 1)):
            unsorted_shards += 1
            problems.append(f'shard {i}: records are not in strict lexicographic order')
        if len(set(recs)) != len(recs):
            dup_records += 1
        row0 = bytes(r['row0'])
        if universe is not None and list(r['row0']) != universe[i]:
            problems.append(f'shard {i}: manifest row 0 is not the universe\'s row for it')
        wrong_row0 += sum(1 for rec in recs if rec[:n] != row0)
        if row0 in rows:
            problems.append(f'shard {i}: row 0 repeats an earlier shard')
        rows.add(row0)
        walls.append(r['wall'])
        nodes.append(r['nodes'])
        counts.append(r['count'])
    out = dict(shards=len(man), supports=sum(counts),
               universe=len(universe) if universe is not None else None,
               covers_the_universe=(universe is not None
                                    and set(man) == set(universe) and not problems),
               shards_with_duplicate_records=dup_records,
               shards_out_of_lexicographic_order=unsorted_shards,
               shards_with_out_of_range_bytes=out_of_range,
               records_with_the_wrong_row0=wrong_row0,
               distinct_shard_rows=len(rows),
               unfinished_records=len(unfinished),
               problems=problems[:20], problem_count=len(problems))
    if walls:
        out.update(core_seconds=round(sum(walls), 1),
                   core_hours=round(sum(walls) / 3600, 2),
                   nodes=sum(nodes),
                   wall_per_shard=dict(min=round(min(walls), 4),
                                       median=round(statistics.median(walls), 4),
                                       max=round(max(walls), 4)),
                   supports_per_shard=dict(min=min(counts),
                                           median=statistics.median(counts),
                                           mean=round(sum(counts) / len(counts), 1),
                                           max=max(counts),
                                           empty=sum(1 for c in counts if c == 0)))
    print(json.dumps(out, indent=1))
    bad = (problems or dup_records or wrong_row0 or unsorted_shards or out_of_range
           or (universe is not None and not out['covers_the_universe']))
    return 1 if bad else 0


# --------------------------------------------------------------------- sample
def cmd_sample(a):
    """A uniform sample of the shards the symmetric mode does NOT enumerate.

    Drawing from the mapped shards only is the point: a representative's payload
    is produced by the enumerator either way, so re-enumerating one would check
    nothing about the mapping."""
    import random
    rows = []
    with open(a.orbits) as f:
        for line in f:
            r = json.loads(line)
            if r['rep'] != r['idx']:
                rows.append((r['idx'], r['row0']))
    if a.k > len(rows):
        sys.exit(f'asked for {a.k} of {len(rows)} mapped shards')
    pick = sorted(random.Random(a.seed).sample(rows, a.k))
    with open(a.out, 'w') as f:
        for idx, row in pick:
            f.write(' '.join(str(v) for v in [idx] + row) + '\n')
    print(json.dumps(dict(mapped_shards=len(rows), sampled=len(pick), seed=a.seed,
                          first=pick[0][0], last=pick[-1][0])))
    return 0


def cmd_setcmpshards(a):
    """Compare the mapped payload and the directly enumerated payload of every
    shard in a sample, as sets (the substantive claim) and byte for byte (which
    also pins the ordering convention)."""
    w = a.n * a.n
    idxs = [int(line.split()[0]) for line in open(a.sample) if line.strip()]
    bad_set, bad_bytes, missing, records = [], [], [], 0
    for i in idxs:
        pa, pb = shard_path(a.a, i), shard_path(a.b, i)
        if not (os.path.exists(pa) and os.path.exists(pb)):
            missing.append(i)
            continue
        ra, rb = open(pa, 'rb').read(), open(pb, 'rb').read()
        sa = {ra[k:k + w] for k in range(0, len(ra), w)}
        sb = {rb[k:k + w] for k in range(0, len(rb), w)}
        records += len(sa)
        if sa != sb:
            bad_set.append(i)
        elif ra != rb:
            bad_bytes.append(i)
    ok = not bad_set and not bad_bytes and not missing
    print(json.dumps(dict(shards=len(idxs), records=records,
                          set_equal=len(idxs) - len(bad_set) - len(missing),
                          byte_identical=len(idxs) - len(bad_set) - len(bad_bytes) - len(missing),
                          differing_as_sets=bad_set[:20], differing_bytes_only=bad_bytes[:20],
                          missing=missing[:20], all_agree=ok)))
    return 0 if ok else 1


# --------------------------------------------------------------------- ledger
DROP = ('search_wall', 'clique_wall', 'wall')


def cmd_ledger(a):
    """Canonicalization only: drop the timings, sort, summarize.

    This deliberately makes no judgement about whether the results are the ones
    the theorem needs -- that is `validate`, kept separate so that neither can
    be mistaken for the other.
    """
    recs = [json.loads(l) for l in open(a.inp) if l.strip()]
    keep = [{k: v for k, v in r.items() if k not in DROP} for r in recs]
    with open(a.out, 'w') as f:
        for r in sorted(keep, key=lambda r: r['j']):
            f.write(json.dumps(r, sort_keys=True, separators=(',', ':')) + '\n')
    status = collections.Counter(r['status'] for r in keep)
    covers = collections.Counter(r['covers'] for r in keep)
    packing = collections.Counter(r['max_packing_with_query'] for r in keep)
    pools = [r['pool'] for r in keep]
    print(json.dumps(dict(
        queries=len(keep), status=dict(status), covers=dict(covers),
        max_packing_with_query=dict(sorted(packing.items())),
        triangle_free_pools=sum(1 for r in keep if r['triangles'] == 0),
        pool_members_meeting_their_query=sum(r['pool_meets_query'] for r in keep),
        max_search_depth=max(r['max_depth'] for r in keep),
        nodes=sum(r['nodes'] for r in keep),
        pool=dict(min=min(pools), median=statistics.median(pools),
                  mean=round(sum(pools) / len(pools), 1), max=max(pools))), indent=1))
    return 0


def cmd_validate(a):
    """The strict result validator: does the ledger say what the theorem needs?

    Every one of these is a way the run could finish and still not entitle
    anybody to the conclusion, so each is required explicitly rather than left
    to a reference digest to notice:

      * exactly the expected queries, each once -- a missing root orbit is the
        failure a matching count of *something* would hide;
      * every search completed; a BUDGET is a truncated search, not a result;
      * no covers found, and no branch surviving past the query plus one
        companion;
      * pools non-empty and disjoint from their query, so they really are
        companion pools;
      * the clique bound: omega(Gamma) <= n-2 < n-1, so no query lies in a
        partition;
      * and the elementary check behind it: an (n-1)-clique contains C(n-1,3)
        triangles, so a pool with fewer cannot have one.  At n = 9 that is 56
        triangles against the 0 or 8 recorded.
    """
    n = a.n
    recs = [json.loads(l) for l in open(a.inp) if l.strip()]
    problems = []
    ids = [r['j'] for r in recs]
    if sorted(ids) != list(range(a.queries)):
        dup = [j for j, c in collections.Counter(ids).items() if c > 1]
        missing = sorted(set(range(a.queries)) - set(ids))
        problems.append(f'query ids are not exactly 0..{a.queries - 1}: '
                        f'{len(recs)} records, {len(missing)} missing {missing[:10]}, '
                        f'{len(dup)} duplicated {dup[:10]}')
    need_triangles = math.comb(n - 1, 3)
    for r in recs:
        j = r['j']
        if r['status'] != 'EXHAUSTED':
            problems.append(f'query {j}: status {r["status"]}, not a completed search')
        if r['covers'] != 0:
            problems.append(f'query {j}: {r["covers"]} covers found')
        if r['pool'] <= 0:
            problems.append(f'query {j}: empty pool')
        if r['pool_meets_query'] != 0:
            problems.append(f'query {j}: {r["pool_meets_query"]} pool members meet the query')
        if r['max_companion_clique'] > n - 2:
            problems.append(f'query {j}: companion clique {r["max_companion_clique"]}, '
                            f'which does not rule out a partition')
        # The graph fields must describe a computed graph, not the placeholders
        # a run with the clique stage switched off leaves: a nonempty pool has a
        # clique of size at least 1, an edge makes it at least 2, a triangle 3.
        w = r['max_companion_clique']
        if r['pool'] > 0 and w < 1:
            problems.append(f'query {j}: clique number {w} for a pool of {r["pool"]}; '
                            f'the clique computation did not run')
        if (r['edges'] > 0) != (w >= 2) or (r['triangles'] > 0) != (w >= 3):
            problems.append(f'query {j}: clique number {w} is inconsistent with '
                            f'{r["edges"]} edges and {r["triangles"]} triangles')
        if r['max_packing_with_query'] != r['max_companion_clique'] + 1:
            problems.append(f'query {j}: packing and clique disagree')
        if r['max_depth'] > 2:
            problems.append(f'query {j}: a branch survived to depth {r["max_depth"]}')
        if r['triangles'] >= need_triangles:
            problems.append(f'query {j}: {r["triangles"]} triangles, enough for an '
                            f'{n - 1}-clique; the elementary bound does not apply')
        if math.comb(r['max_companion_clique'], 3) > r['triangles']:
            problems.append(f'query {j}: a {r["max_companion_clique"]}-clique needs '
                            f'{math.comb(r["max_companion_clique"], 3)} triangles, '
                            f'{r["triangles"]} recorded')
    cl = collections.Counter(r['max_companion_clique'] for r in recs)
    tr = collections.Counter(r['triangles'] for r in recs)
    print(json.dumps(dict(
        n=n, queries=len(recs), expected=a.queries,
        all_searches_completed=all(r['status'] == 'EXHAUSTED' for r in recs),
        total_covers=sum(r['covers'] for r in recs),
        max_search_depth=max(r['max_depth'] for r in recs),
        companion_clique_histogram=dict(sorted(cl.items())),
        triangle_histogram=dict(sorted(tr.items())),
        triangles_needed_for_a_partition=need_triangles,
        valid=not problems, problems=problems[:20], problem_count=len(problems)),
        indent=1))
    if problems:
        print(f'{len(problems)} problems; the result does NOT support the theorem',
              file=sys.stderr)
    return 1 if problems else 0


def cmd_orbitstats(a):
    recs = [json.loads(l) for l in open(a.inp) if l.strip()]
    hist = collections.Counter(r['orbit_size'] for r in recs)
    print(json.dumps(dict(orbits=len(recs), sum_orbit_sizes=sum(r['orbit_size'] for r in recs),
                          size_histogram=dict(sorted(hist.items()))), indent=1))
    return 0


# -------------------------------------------------------------------- reports
def cmd_reports(a):
    """Require one complete terminal report per worker, and full coverage.

    A parallel stage whose workers write their reports to files fails open by
    default: a worker that is killed, or runs out of memory, leaves a short file
    or no file, and a grep for the bad news then finds nothing and says so.  So
    the reports are counted, parsed, and required to cover the whole universe of
    work between them.
    """
    paths = sorted(glob.glob(os.path.join(a.dir, a.pattern)))
    problems = []
    if len(paths) != a.parts:
        problems.append(f'{len(paths)} report files, expected {a.parts}: '
                        f'{[os.path.basename(p) for p in paths]}')
    seen_parts, ids, records = set(), [], 0
    for p in paths:
        text = open(p).read().strip()
        if not text:
            problems.append(f'{os.path.basename(p)}: empty -- the worker left no report')
            continue
        try:
            r = json.loads(text.splitlines()[-1])
        except ValueError as e:
            problems.append(f'{os.path.basename(p)}: no terminal JSON report ({e})')
            continue
        if not isinstance(r, dict) or 'part' not in r or 'parts' not in r:
            problems.append(f'{os.path.basename(p)}: report does not say which slice it is')
            continue
        if r['parts'] != a.parts:
            problems.append(f'{os.path.basename(p)}: claims {r["parts"]} slices, driver ran {a.parts}')
        if r['part'] in seen_parts:
            problems.append(f'{os.path.basename(p)}: slice {r["part"]} reported twice')
        seen_parts.add(r['part'])
        if a.kind == 'pools':
            if not isinstance(r.get('disagreements'), int) or r['disagreements'] != 0:
                problems.append(f'{os.path.basename(p)}: {r["disagreements"]} pools '
                                f'disagree with the re-derivation {r.get("disagreeing")}')
            if not r.get('complete'):
                problems.append(f'{os.path.basename(p)}: report is not marked complete')
            qi = r.get('query_ids')
            if not isinstance(qi, list):
                problems.append(f'{os.path.basename(p)}: no query ids')
            else:
                if any(j % a.parts != r['part'] for j in qi):
                    problems.append(f'{os.path.basename(p)}: query ids are not its own slice')
                ids.extend(qi)
        else:
            bad = (r.get('records_without_n2_distinct_cells', 1)
                   + r.get('records_missing_a_line_exactly_once', 1))
            if bad:
                problems.append(f'{os.path.basename(p)}: {bad} records failed the definition')
            records += r.get('checked', 0)
    if seen_parts != set(range(a.parts)):
        problems.append(f'slices reported: {sorted(seen_parts)}, expected 0..{a.parts - 1}')
    if a.kind == 'pools':
        if sorted(ids) != list(range(a.expect)):
            missing = sorted(set(range(a.expect)) - set(ids))
            problems.append(f'the workers between them checked {len(ids)} queries, '
                            f'not the expected {a.expect}; missing {missing[:10]}')
    elif records != a.expect:
        problems.append(f'{records} records checked, expected {a.expect}')
    print(json.dumps(dict(kind=a.kind, reports=len(paths), workers=a.parts,
                          queries_checked=len(ids) if a.kind == 'pools' else None,
                          records_checked=records if a.kind != 'pools' else None,
                          expected=a.expect, ok=not problems,
                          problems=problems[:20]), indent=1))
    return 1 if problems else 0


# ------------------------------------------------------------------ checksums
SHA256 = re.compile(r'\A[0-9a-f]{64}\Z')


def read_checksums(path):
    d = {}
    for lineno, line in enumerate(open(path), 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) != 2:
            sys.exit(f'{path} line {lineno}: expected "<sha256>  <name>"')
        h, name = parts[0], os.path.basename(parts[1])
        if not SHA256.match(h):
            sys.exit(f'{path} line {lineno}: {h!r} is not a SHA-256 digest')
        if name in d:
            sys.exit(f'{path} line {lineno}: {name} appears twice')
        d[name] = h
    return d


def cmd_checksums(a):
    """Compare two checksum manifests, requiring exactly the expected names.

    Silently overwriting a duplicate name, or comparing whatever names happen to
    be in both files, would let an artifact go unchecked without anyone
    noticing; the set of names is therefore required, not inferred.
    """
    want, got = read_checksums(a.want), read_checksums(a.got)
    bad = 0
    required = set(a.require)
    for missing in sorted(required - set(want)):
        print(f'NOT LISTED  {missing}  (expected in {a.want})'); bad += 1
    for extra in sorted(set(want) - required):
        print(f'UNEXPECTED  {extra}  (in {a.want}, not in the required set)'); bad += 1
    for name in sorted(want):
        if name not in got:
            print(f'MISSING   {name}'); bad += 1
        elif want[name] != got[name]:
            print(f'MISMATCH  {name}\n  expected {want[name]}\n  actual   {got[name]}'); bad += 1
        else:
            print(f'OK        {name}  {got[name]}')
    for extra in sorted(set(got) - set(want)):
        print(f'UNLISTED  {extra}  was produced but is not in {a.want}'); bad += 1
    print(f'\n{len(want) - bad} of {len(want)} artifacts match')
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('pack', help='assemble and checksum the catalogue')
    p.add_argument('n', type=int); p.add_argument('manifest')
    p.add_argument('sharddir'); p.add_argument('out')
    p.add_argument('--shards', help='the expected shard universe (`shards list` output)')
    p.set_defaults(fn=cmd_pack)

    p = sub.add_parser('audit', help='full audit against the expected shard universe')
    p.add_argument('n', type=int); p.add_argument('manifest'); p.add_argument('sharddir')
    p.add_argument('--shards', help='the expected shard universe (`shards list` output)')
    p.add_argument('--partial', action='store_true',
                   help='audit only what the manifest holds; completeness is NOT checked')
    p.set_defaults(fn=cmd_audit)

    p = sub.add_parser('ledger', help='drop timings, sort, and summarize a result file')
    p.add_argument('inp'); p.add_argument('out')
    p.set_defaults(fn=cmd_ledger)

    p = sub.add_parser('validate', help='the strict result validator')
    p.add_argument('n', type=int); p.add_argument('inp')
    p.add_argument('--queries', type=int, required=True)
    p.set_defaults(fn=cmd_validate)

    p = sub.add_parser('orbitstats', help='orbit-size histogram')
    p.add_argument('inp'); p.set_defaults(fn=cmd_orbitstats)

    p = sub.add_parser('sample', help='a uniform sample of the mapped shards')
    p.add_argument('n', type=int); p.add_argument('orbits')
    p.add_argument('k', type=int); p.add_argument('out')
    p.add_argument('--seed', type=int, default=20260906)
    p.set_defaults(fn=cmd_sample)

    p = sub.add_parser('setcmpshards', help='compare two shard trees on a sample')
    p.add_argument('n', type=int); p.add_argument('sample')
    p.add_argument('a'); p.add_argument('b')
    p.set_defaults(fn=cmd_setcmpshards)

    p = sub.add_parser('reports', help='require one complete report per worker')
    p.add_argument('kind', choices=('pools', 'verify'))
    p.add_argument('dir'); p.add_argument('pattern')
    p.add_argument('--parts', type=int, required=True)
    p.add_argument('--expect', type=int, required=True,
                   help='queries (pools) or records (verify) the workers must cover')
    p.set_defaults(fn=cmd_reports)

    p = sub.add_parser('checksums', help='compare two checksum manifests strictly')
    p.add_argument('want'); p.add_argument('got')
    p.add_argument('--require', nargs='+', required=True,
                   help='the exact set of artifact names the manifest must list')
    p.set_defaults(fn=cmd_checksums)

    a = ap.parse_args()
    try:
        sys.exit(a.fn(a))
    except ManifestError as e:
        sys.exit(f'manifest: {e}')


if __name__ == '__main__':
    main()
