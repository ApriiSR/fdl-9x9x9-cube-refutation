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

    a = ap.parse_args()
    sys.exit(a.fn(a))


if __name__ == '__main__':
    main()
