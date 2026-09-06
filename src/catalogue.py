#!/usr/bin/env python3
"""catalogue.py -- assemble and audit the support catalogue.

The catalogue is the concatenation of the per-shard payload files in increasing
shard index.  Since `enum` sorts the records inside a shard, and shard index is
the lexicographic rank of the shard's row 0, the catalogue is exactly the list
of all supports of [n]^3 in lexicographic order: a deterministic function of the
mathematics, with no dependence on the search order.  Its SHA-256 is therefore
a meaningful thing to publish and to compare across independent runs.

Subcommands:

  pack N MANIFEST SHARDDIR OUT.bin   assemble and checksum the catalogue
  audit N MANIFEST SHARDDIR          manifest-vs-disk audit and sweep statistics
  ledger IN.jsonl OUT.jsonl          drop timing fields, leaving a reproducible
                                     per-query result file, and summarise it
  orbitstats IN.jsonl                orbit-size histogram and total
  sample N ORBITS K OUT.txt          K mapped shards, drawn uniformly, in the
                                     format `shards list` writes
  setcmpshards N ORBITS A B          compare two shard trees on the sampled
                                     shards, as sets and byte for byte
"""
import argparse
import collections
import hashlib
import json
import os
import statistics
import sys


def read_manifest(path):
    """{idx: record} for shards marked done.  A torn last line is ignored, so a
    kill during a write costs at most the shard that was in flight."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get('done'):
                out[r['idx']] = r
    return out


def shard_path(base, idx):
    return os.path.join(base, f'{idx // 1000:03d}', f's{idx:05d}.bin')


def cmd_pack(a):
    n = a.n
    w = n * n
    man = read_manifest(a.manifest)
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
                order='shard index, records sorted within a shard')
    with open(a.out + '.meta.json', 'w') as f:
        json.dump(meta, f, indent=1)
        f.write('\n')
    print(json.dumps(meta))
    return 0


def cmd_audit(a):
    n = a.n
    w = n * n
    man = read_manifest(a.manifest)
    problems = []
    rows = set()
    dup_records = 0
    wrong_row0 = 0
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
        if len(set(recs)) != len(recs):
            dup_records += 1
        row0 = bytes(r['row0'])
        wrong_row0 += sum(1 for rec in recs if rec[:n] != row0)
        if row0 in rows:
            problems.append(f'shard {i}: row 0 repeats an earlier shard')
        rows.add(row0)
        walls.append(r['wall'])
        nodes.append(r['nodes'])
        counts.append(r['count'])
    out = dict(shards=len(man), supports=sum(counts),
               shards_with_duplicate_records=dup_records,
               records_with_the_wrong_row0=wrong_row0,
               distinct_shard_rows=len(rows),
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
    return 0 if not problems and not dup_records and not wrong_row0 else 1


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


DROP = ('search_wall', 'clique_wall', 'wall')


def cmd_ledger(a):
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


def cmd_orbitstats(a):
    recs = [json.loads(l) for l in open(a.inp) if l.strip()]
    hist = collections.Counter(r['orbit_size'] for r in recs)
    print(json.dumps(dict(orbits=len(recs), sum_orbit_sizes=sum(r['orbit_size'] for r in recs),
                          size_histogram=dict(sorted(hist.items()))), indent=1))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('pack', help='assemble and checksum the catalogue')
    p.add_argument('n', type=int); p.add_argument('manifest')
    p.add_argument('sharddir'); p.add_argument('out')
    p.set_defaults(fn=cmd_pack)

    p = sub.add_parser('audit', help='manifest-vs-disk audit and sweep statistics')
    p.add_argument('n', type=int); p.add_argument('manifest'); p.add_argument('sharddir')
    p.set_defaults(fn=cmd_audit)

    p = sub.add_parser('ledger', help='drop timings, sort, and summarise a result file')
    p.add_argument('inp'); p.add_argument('out')
    p.set_defaults(fn=cmd_ledger)

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

    a = ap.parse_args()
    sys.exit(a.fn(a))


if __name__ == '__main__':
    main()
