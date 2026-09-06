#!/usr/bin/env python3
"""The order-8 controls.

Order 8 is where a fully diagonalised Latin cube is known to exist and its
census is known independently, so every step the order-9 argument runs through
can be checked against an answer that was not produced by this code:

    13 056 supports, 198 624 exact covers, 6 orbits under the cell group,
    and, for every one of the 13 056 supports, the number of covers containing
    it -- computed once by the census and once, support by support, by the
    same root search that is run at order 9.

It also runs the whole of `verify.sh symmetric` at order 8: 25 of the 5 568
shards enumerated, the rest produced by the plane-fixing subgroup of Lemma 6,
and the result required to be the census byte for byte.

The suite also checks that the time caps fire, and that the mapping's guard
fires, which is the control a negative result most needs: a cap that never
expires looks exactly like an exhaustion, and a check that never rejects looks
exactly like agreement.

Usage:  python3 tests/test_n8.py [--workers W] [--keep]
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, 'bin')
DATA = os.path.join(ROOT, 'data')
PY = sys.executable

PASS, FAIL = [], []


def check(name, ok, detail=''):
    (PASS if ok else FAIL).append(name)
    print(f'{"PASS" if ok else "FAIL"}  {name}{"  -- " + detail if detail else ""}',
          flush=True)


def run(args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 4)
    ap.add_argument('--keep', action='store_true', help='keep the scratch directory')
    a = ap.parse_args()
    W = max(1, a.workers)
    tmp = tempfile.mkdtemp(prefix='fdlh8-')
    t_all = time.time()
    print(f'scratch: {tmp}\nworkers: {W}\n')

    # ---- 1. the line structure ------------------------------------------
    for n, want in ((8, 244), (9, 301)):
        r = run([f'{BIN}/enum', str(n), 'lines'])
        got = int(r.stdout.split('lines=')[1].split()[0])
        p = run([PY, f'{ROOT}/src/check.py', 'lines', str(n)])
        got2 = int(p.stdout.split('lines=')[1].split()[0])
        check(f'main lines at n={n}: {want}, C and numpy agree',
              r.returncode == 0 and p.returncode == 0 and got == want == got2,
              f'{got} / {got2}')

    # ---- 2. the shard universe ------------------------------------------
    r = run([f'{BIN}/shards', '9', 'count'])
    counts = {int(l.split('\t')[0]): int(l.split('\t')[1])
              for l in r.stdout.strip().splitlines()[1:]}
    closed = json.loads(run([PY, f'{ROOT}/src/check.py', 'a007016', '9']).stdout)['sequence']
    oeis = [1, 0, 0, 8, 20, 96, 656, 5568, 48912]      # OEIS A007016, a(1..9)
    check('shard universe is A007016 (brute force = closed form = OEIS)',
          closed == oeis and all(counts[n] == oeis[n - 1] for n in range(2, 10)),
          f'n=8: {counts[8]}, n=9: {counts[9]}')

    shardfile = f'{tmp}/n8_shards.txt'
    run([f'{BIN}/shards', '8', 'list', shardfile])

    # ---- 3. enumerate T(8) ----------------------------------------------
    t0 = time.time()
    procs = [subprocess.Popen([f'{BIN}/enum', '8', 'shards', shardfile,
                               f'{tmp}/shards', f'{tmp}/man_{k}.jsonl',
                               '--slice', str(k), str(W)],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
             for k in range(W)]
    rcs = [p.wait() for p in procs]
    t_enum = time.time() - t0
    with open(f'{tmp}/manifest.jsonl', 'w') as out:
        for k in range(W):
            out.write(open(f'{tmp}/man_{k}.jsonl').read())
    aud = run([PY, f'{ROOT}/src/catalogue.py', 'audit', '8',
               f'{tmp}/manifest.jsonl', f'{tmp}/shards'])
    A = json.loads(aud.stdout)
    check('T(8): all 5 568 shards exhausted, 13 056 supports',
          all(r == 0 for r in rcs) and A['shards'] == 5568 and A['supports'] == 13056,
          f'{A["supports"]} supports in {t_enum:.1f} s wall, '
          f'{A.get("core_seconds")} core-s')
    check('no duplicate records, every record carries its own shard row 0',
          A['shards_with_duplicate_records'] == 0 and A['records_with_the_wrong_row0'] == 0
          and A['distinct_shard_rows'] == 5568)

    cat8 = f'{tmp}/n8_supports.bin'
    meta = json.loads(run([PY, f'{ROOT}/src/catalogue.py', 'pack', '8',
                           f'{tmp}/manifest.jsonl', f'{tmp}/shards', cat8]).stdout)
    ref = os.path.join(DATA, 'n8_supports.bin')
    if os.path.exists(ref):
        s = run([PY, f'{ROOT}/src/check.py', 'setcmp', '8', cat8, ref])
        eq = json.loads(s.stdout)
        check('the enumeration equals the shipped order-8 census, as a set and byte for byte',
              eq['set_equal'] and sha256(cat8) == sha256(ref), meta['sha256'])
    else:
        check('data/n8_supports.bin is present', False, 'missing')

    v = json.loads(run([PY, f'{ROOT}/src/check.py', 'verify', '8', cat8]).stdout)
    check('all 13 056 records are supports by the definition',
          v['records_without_n2_distinct_cells'] == 0
          and v['records_missing_a_line_exactly_once'] == 0,
          f'{v["checked"]} checked against {v["lines"]} main lines')

    # ---- 4. the cell group ----------------------------------------------
    g = run([f'{BIN}/orbits', '8', 'group'])
    check('the cell group has order 9 216', g.returncode == 0 and '|G|=9216' in g.stdout)
    c = run([f'{BIN}/orbits', '8', 'closure', cat8])
    C = json.loads(c.stdout)
    check('T(8) is closed under the cell group', C['missing'] == 0,
          f'{C["images_checked"]} images checked')
    o = run([f'{BIN}/orbits', '8', 'classify', cat8, f'{tmp}/orbits.jsonl', f'{tmp}/reps.bin'])
    O = json.loads(o.stdout)
    hist = json.loads(run([PY, f'{ROOT}/src/catalogue.py', 'orbitstats',
                           f'{tmp}/orbits.jsonl']).stdout)
    check('T(8) is 6 orbits: two of size 768, three of 2 304, one of 4 608',
          O['orbits'] == 6 and O['sum_orbit_sizes'] == 13056
          and hist['size_histogram'] == {'768': 2, '2304': 3, '4608': 1},
          json.dumps(hist['size_histogram']))

    # ---- 5. the census ---------------------------------------------------
    t0 = time.time()
    cen = run([f'{BIN}/pack', '8', 'census', cat8, f'{tmp}/n8_covers.json'])
    t_cen = time.time() - t0
    CEN = json.loads(cen.stdout)
    profile = [1, 1632, 26016, 60672, 69513, 147292, 154660, 198624, 198624]
    check('198 624 exact covers of [8]^3, with the published node profile',
          CEN['covers'] == 198624 and CEN['nodes_by_depth'] == profile,
          f'{t_cen:.1f} s')
    tallies = json.load(open(f'{tmp}/n8_covers.json'))['covers_containing']
    check('cover incidences add up: 198 624 x 8', sum(tallies) == 198624 * 8)

    # ---- 6. every root case against the census ---------------------------
    t0 = time.time()
    run([f'{BIN}/pools', '8', 'build', cat8, cat8, f'{tmp}/pools', f'{tmp}/poolsizes.jsonl'])
    procs = [subprocess.Popen([f'{BIN}/pack', '8', 'roots', cat8, f'{tmp}/pools',
                               f'{tmp}/roots_{k}.jsonl', '--slice', str(k), str(W)],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
             for k in range(W)]
    [p.wait() for p in procs]
    t_roots = time.time() - t0
    res = {}
    for k in range(W):
        for line in open(f'{tmp}/roots_{k}.jsonl'):
            r = json.loads(line)
            res[r['j']] = r
    agree = sum(1 for j, r in res.items() if r['covers'] == tallies[j])
    check('root search agrees with the census on every one of the 13 056 supports',
          len(res) == 13056 and agree == 13056,
          f'{agree}/{len(res)} in {t_roots:.1f} s')
    eight = [j * 1632 for j in range(8)]
    check('the eight sampled root cases, spelled out',
          all(res[j]['covers'] == tallies[j] for j in eight),
          ' '.join(f'{j}:{tallies[j]}' for j in eight))
    check('no pool member meets its own query',
          all(r['pool_meets_query'] == 0 for r in res.values()))
    check('the packing ceiling agrees with the search: a support lies in a cover '
          'iff its pool graph has a 7-clique',
          all((r['max_companion_clique'] == 7) == (r['covers'] > 0) for r in res.values()))

    # ---- 7. pools re-derived independently in numpy ----------------------
    import numpy as np
    allrec = np.fromfile(cat8, dtype=np.uint8).reshape(-1, 64)
    allrec[eight].tofile(f'{tmp}/q8.bin')
    run([f'{BIN}/pools', '8', 'build', cat8, f'{tmp}/q8.bin', f'{tmp}/pools8',
         f'{tmp}/poolsizes8.jsonl'])
    p = run([PY, f'{ROOT}/src/check.py', 'pools', '8', cat8, f'{tmp}/q8.bin', f'{tmp}/pools8'])
    P = json.loads(p.stdout.strip().splitlines()[-1])
    check('companion pools re-derived in numpy agree exactly', P['disagreements'] == 0,
          f'sizes {P["pool_min"]}..{P["pool_max"]}')

    # ---- 8. an explicit witness -----------------------------------------
    j0 = next(j for j in eight if tallies[j] > 0)
    allrec[[j0]].tofile(f'{tmp}/q1.bin')
    run([f'{BIN}/pools', '8', 'build', cat8, f'{tmp}/q1.bin', f'{tmp}/pools1', f'{tmp}/ps1.jsonl'])
    run([f'{BIN}/pack', '8', 'witness', f'{tmp}/q1.bin', f'{tmp}/pools1', '0', f'{tmp}/w.json'])
    wj = run([PY, f'{ROOT}/src/check.py', 'witness', '8', f'{tmp}/w.json'])
    WI = json.loads(wj.stdout)
    check('the witness is 8 pairwise disjoint supports -- an order-8 cube',
          WI['squares'] == 8 and WI['pairwise_disjoint'])

    # ---- 9. the stopping rules fire --------------------------------------
    head = f'{tmp}/head20.txt'
    with open(head, 'w') as f:
        f.writelines(open(shardfile).readlines()[:20])
    capdir, capman = f'{tmp}/capshards', f'{tmp}/capman.jsonl'
    run([f'{BIN}/enum', '8', 'shards', head, capdir, capman, '--cap', '1e-9'])
    lines = [json.loads(l) for l in open(capman)]
    left = [f for _, _, fs in os.walk(capdir) for f in fs]
    check('enum: the wall-clock cap fires, banks nothing, and leaves no payload',
          len(lines) == 20 and all(l['status'] == 'BUDGET' and not l['done'] for l in lines)
          and not left, f'{len(lines)} shards capped, {len(left)} files left behind')
    run([f'{BIN}/enum', '8', 'shards', head, capdir, capman, '--cap', '60'])
    lines = [json.loads(l) for l in open(capman)]
    check('enum: a capped shard is redone on the next pass',
          sum(1 for l in lines if l.get('done')) == 20)

    cap = run([f'{BIN}/pack', '8', 'census', cat8, f'{tmp}/cap.json', '--cap', '0.05'])
    CAP = json.loads(cap.stdout)
    check('pack: the cap fires and is reported as BUDGET, never as EXHAUSTED',
          CAP['status'] == 'BUDGET' and 0 < CAP['covers'] < 198624,
          f'{CAP["covers"]} of 198 624 covers found before the cap')

    rc = run([f'{BIN}/pack', '8', 'roots', f'{tmp}/q1.bin', f'{tmp}/pools1',
              f'{tmp}/capr.jsonl', '--cap', '1e-9'])
    R = json.loads(open(f'{tmp}/capr.jsonl').read().strip())
    check('pack roots: the per-query cap fires', R['status'] == 'BUDGET', str(rc.returncode))

    # ---- 10. the plane-fixing subgroup, and the symmetric mode ----------
    g = run([f'{BIN}/symmetry', '8', 'group'])
    check('the plane-fixing subgroup has order 384, index 24, and is closed',
          g.returncode == 0 and '|H|=384' in g.stdout and 'index=24' in g.stdout
          and 'leaving the plane: 0' in g.stdout,
          g.stdout.splitlines()[0] if g.stdout else '')

    symorb, symreps = f'{tmp}/n8_shard_orbits.jsonl', f'{tmp}/n8_reps.txt'
    o = run([f'{BIN}/symmetry', '8', 'orbits', shardfile, symorb, symreps])
    SO = json.loads(o.stdout)
    check('the 5 568 shards are 25 orbits, and every image is an admissible '
          'row that is in the shard list',
          o.returncode == 0 and SO['orbits'] == 25 and SO['sum_orbit_sizes'] == 5568
          and SO['representatives'] == 25,
          json.dumps(SO['orbit_size_histogram']))

    symdir = f'{tmp}/symshards'
    os.makedirs(symdir, exist_ok=True)
    run([f'{BIN}/enum', '8', 'shards', symreps, symdir, f'{tmp}/symman_reps.jsonl'])
    e = run([f'{BIN}/symmetry', '8', 'expand', symorb, symdir, f'{tmp}/symman_map.jsonl'])
    with open(f'{tmp}/symman.jsonl', 'w') as f:
        f.write(open(f'{tmp}/symman_reps.jsonl').read())
        f.write(open(f'{tmp}/symman_map.jsonl').read())
    symcat = f'{tmp}/n8_supports_sym.bin'
    run([PY, f'{ROOT}/src/catalogue.py', 'pack', '8', f'{tmp}/symman.jsonl', symdir, symcat])
    check('the symmetric mode reproduces the whole order-8 census from 25 of '
          'the 5 568 shards, byte for byte',
          e.returncode == 0 and os.path.exists(symcat)
          and sha256(symcat) == sha256(cat8) == sha256(ref),
          f'{json.loads(e.stdout)["records"]} records mapped' if e.stdout else '')

    # A mapping that is never rejected is not a check.  Point one shard at a
    # different element of the subgroup and the guard must fire: in an orbit of
    # full size the stabiliser is trivial, so any other element lands the
    # representative's records on a different row 0.
    lines = [json.loads(l) for l in open(symorb)]
    big = {r['orbit'] for r in lines
           if sum(1 for x in lines if x['orbit'] == r['orbit']) == 384}
    def repsize(i):
        q = os.path.join(symdir, f'{i // 1000:03d}', f's{i:05d}.bin')
        return os.path.getsize(q) if os.path.exists(q) else 0
    victim = next(r for r in lines
                  if r['orbit'] in big and r['rep'] != r['idx'] and repsize(r['rep']))
    victim['g'] = (victim['g'] + 1) % 384
    with open(f'{tmp}/bad_orbits.jsonl', 'w') as f:
        f.write(json.dumps(victim, separators=(',', ':')) + '\n')
    b = run([f'{BIN}/symmetry', '8', 'expand', f'{tmp}/bad_orbits.jsonl', symdir,
             f'{tmp}/badman.jsonl'])
    check('symmetry expand rejects a record image that does not carry the '
          "shard's own row 0",
          b.returncode != 0 and 'row 0' in b.stderr,
          b.stderr.strip().splitlines()[0] if b.stderr.strip() else 'no message')

    # ---- 11. determinism -------------------------------------------------
    run([f'{BIN}/enum', '8', 'one', f'{tmp}/d1.bin'] + open(shardfile).readline().split()[1:])
    run([f'{BIN}/enum', '8', 'one', f'{tmp}/d2.bin'] + open(shardfile).readline().split()[1:])
    check('a shard is byte-identical when re-enumerated',
          sha256(f'{tmp}/d1.bin') == sha256(f'{tmp}/d2.bin'))

    print(f'\n{len(PASS)} passed, {len(FAIL)} failed, {time.time() - t_all:.1f} s')
    if FAIL:
        print('failed: ' + ', '.join(FAIL))
    if not a.keep:
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print(f'kept {tmp}')
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
