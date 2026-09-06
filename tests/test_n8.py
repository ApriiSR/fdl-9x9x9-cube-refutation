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


def failing(args, **kw):
    """Run something that MUST refuse: nonzero exit, and a message saying why."""
    r = subprocess.run(args, capture_output=True, text=True, **kw)
    return r.returncode != 0, r


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
    ap.add_argument('--tmpdir', help='make the scratch directory here '
                                     '(default: the system temporary directory)')
    a = ap.parse_args()
    W = max(1, a.workers)
    if a.tmpdir:
        os.makedirs(a.tmpdir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='fdlh8-', dir=a.tmpdir)
    t_all = time.time()
    print(f'scratch: {tmp}\nworkers: {W}\n')

    # ---- 1. the line structure ------------------------------------------
    for n, want in ((8, 244), (9, 301)):
        cdump, pdump = f'{tmp}/lines_c_{n}.txt', f'{tmp}/lines_py_{n}.txt'
        r = run([f'{BIN}/enum', str(n), 'lines', '--dump', cdump])
        got = int(r.stdout.split('lines=')[1].split()[0])
        p = run([PY, f'{ROOT}/src/check.py', 'lines', str(n), '--dump', pdump])
        got2 = int(p.stdout.split('lines=')[1].split()[0])
        # The two constructions are compared as SETS of lines, in canonical form
        # -- each line's cells in increasing order, the lines sorted -- not by
        # their counts, which two different line sets can share.
        same = (os.path.exists(cdump) and os.path.exists(pdump)
                and open(cdump).read() == open(pdump).read())
        check(f'main lines at n={n}: {want}, and the C and numpy line SETS are identical',
              r.returncode == 0 and p.returncode == 0 and got == want == got2 and same,
              f'{got} / {got2}, {"same" if same else "DIFFERENT"} incidence structure')

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
               f'{tmp}/manifest.jsonl', f'{tmp}/shards', '--shards', shardfile])
    A = json.loads(aud.stdout)
    check('T(8): all 5 568 shards exhausted, 13 056 supports',
          all(r == 0 for r in rcs) and A['shards'] == 5568 and A['supports'] == 13056,
          f'{A["supports"]} supports in {t_enum:.1f} s wall, '
          f'{A.get("core_seconds")} core-s')
    check('the audit covers the shard universe exactly, with no duplicate record, '
          'every record carrying its own shard row 0, and every shard in '
          'lexicographic order',
          aud.returncode == 0 and A['covers_the_universe']
          and A['universe'] == 5568 and A['problem_count'] == 0
          and A['shards_with_duplicate_records'] == 0
          and A['shards_out_of_lexicographic_order'] == 0
          and A['shards_with_out_of_range_bytes'] == 0
          and A['records_with_the_wrong_row0'] == 0
          and A['distinct_shard_rows'] == 5568)

    cat8 = f'{tmp}/n8_supports.bin'
    meta = json.loads(run([PY, f'{ROOT}/src/catalogue.py', 'pack', '8',
                           f'{tmp}/manifest.jsonl', f'{tmp}/shards', cat8,
                           '--shards', shardfile]).stdout)
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
    # The caps are driven by an INJECTED clock (FDLH_CLOCK_STEP: every reading
    # advances a virtual clock by exactly that much), so a cap fires after an
    # exact number of readings rather than "eventually".  A cap test that
    # depends on wall time is a test that a loaded machine can turn green or red
    # for reasons that have nothing to do with the code.
    def clocked(step):
        e = dict(os.environ)
        e['FDLH_CLOCK_STEP'] = str(step)
        return e

    head = f'{tmp}/head20.txt'
    with open(head, 'w') as f:
        f.writelines(open(shardfile).readlines()[:20])
    capdir, capman = f'{tmp}/capshards', f'{tmp}/capman.jsonl'
    run([f'{BIN}/enum', '8', 'shards', head, capdir, capman, '--cap', '0.5'],
        env=clocked(1.0))
    lines = [json.loads(l) for l in open(capman)]
    left = [f for _, _, fs in os.walk(capdir) for f in fs]
    check('enum: the wall-clock cap fires, banks nothing, and leaves no payload',
          len(lines) == 20 and all(l['status'] == 'BUDGET' and not l['done'] for l in lines)
          and not left, f'{len(lines)} shards capped, {len(left)} files left behind')
    run([f'{BIN}/enum', '8', 'shards', head, capdir, capman, '--cap', '3600'])
    lines = [json.loads(l) for l in open(capman)]
    check('enum: a capped shard is redone on the next pass',
          sum(1 for l in lines if l.get('done')) == 20)

    # 5 000 readings into the census: enough that some covers have been found,
    # and the same number every time, because the clock is the injected one.
    cap = run([f'{BIN}/pack', '8', 'census', cat8, f'{tmp}/cap.json', '--cap', '0.5'],
              env=clocked(0.001))
    CAP = json.loads(cap.stdout)
    check('pack census: the cap fires at a deterministic point and is reported as '
          'BUDGET, never as EXHAUSTED',
          CAP['status'] == 'BUDGET' and 0 < CAP['covers'] < 198624,
          f'{CAP["covers"]} of 198 624 covers found before the cap')

    rc = run([f'{BIN}/pack', '8', 'roots', f'{tmp}/q1.bin', f'{tmp}/pools1',
              f'{tmp}/capr.jsonl', '--cap', '0.5'], env=clocked(1.0))
    R = json.loads(open(f'{tmp}/capr.jsonl').read().strip())
    check('pack roots: the per-query cap fires', R['status'] == 'BUDGET', str(rc.returncode))

    # The cap is a deadline on the WHOLE query, so it must also be able to fire
    # while the clique bound is being computed rather than the exact cover.
    rc = run([f'{BIN}/pack', '8', 'roots', f'{tmp}/q1.bin', f'{tmp}/pools1',
              f'{tmp}/capr2.jsonl', '--no-search', '--cap', '0.5'], env=clocked(1.0))
    R2 = json.loads(open(f'{tmp}/capr2.jsonl').read().strip())
    check('pack roots: the cap covers the clique stage too, not only the search',
          R2['status'] == 'BUDGET', f'status {R2["status"]} with the search skipped')

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
    run([PY, f'{ROOT}/src/catalogue.py', 'pack', '8', f'{tmp}/symman.jsonl', symdir,
         symcat, '--shards', shardfile])
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

    # ---- 11. the failure paths -------------------------------------------
    # Everything above checks that good input produces the right answer.  This
    # section checks that bad input produces a refusal -- which is the half a
    # negative result actually depends on, since a check that cannot fail says
    # nothing when it passes.

    # (a) a record byte that is not a coordinate in [n] would index outside the
    # cell map, the masks and the image buffers of every program that read it.
    bad = bytearray(open(cat8, 'rb').read(64 * 4))
    bad[70] = 9                                    # a symbol of 9 in [8]^3
    open(f'{tmp}/corrupt.bin', 'wb').write(bytes(bad))
    refusals = []
    for name, args in (
            ('pools build', [f'{BIN}/pools', '8', 'build', f'{tmp}/corrupt.bin',
                             f'{tmp}/q1.bin', f'{tmp}/badpools', f'{tmp}/badsizes.jsonl']),
            ('pack census', [f'{BIN}/pack', '8', 'census', f'{tmp}/corrupt.bin',
                             f'{tmp}/badcensus.json']),
            ('pack roots', [f'{BIN}/pack', '8', 'roots', f'{tmp}/corrupt.bin',
                            f'{tmp}/pools1', f'{tmp}/badroots.jsonl']),
            ('orbits closure', [f'{BIN}/orbits', '8', 'closure', f'{tmp}/corrupt.bin']),
            ('orbits classify', [f'{BIN}/orbits', '8', 'classify', f'{tmp}/corrupt.bin',
                                 f'{tmp}/bad.jsonl', f'{tmp}/bad.bin'])):
        ok, r = failing(args)
        refusals.append((name, ok, 'coordinate' in r.stderr))
    check('a corrupt record byte is refused by every program that loads records, '
          'before it is used as an index',
          all(ok and said for _, ok, said in refusals),
          ', '.join(n for n, ok, said in refusals if not (ok and said)) or 'all five refuse')

    # (b) an order beyond the compiled capacities must be refused up front, not
    # after an array has already been overrun.
    over = []
    for name, args in (('orbits', [f'{BIN}/orbits', '10', 'group']),
                       ('symmetry', [f'{BIN}/symmetry', '10', 'group']),
                       ('shards', [f'{BIN}/shards', '10', 'list', f'{tmp}/s10.txt'])):
        ok, r = failing(args)
        over.append((name, ok, 'supported' in r.stderr or 'holds' in r.stderr))
    check('an unsupported order is refused before construction, with the required '
          'and available capacities named',
          all(ok and said for _, ok, said in over),
          ', '.join(n for n, ok, said in over if not (ok and said)) or 'all three refuse')

    # (c) a pool whose members meet their own query is not a companion pool;
    # searching it would answer a different question.
    # A support always meets itself, so a "pool" containing its own query is the
    # simplest thing that is not a companion pool.
    os.makedirs(f'{tmp}/selfpool', exist_ok=True)
    shutil.copyfile(f'{tmp}/q1.bin', f'{tmp}/selfpool/pool_0.bin')
    ok, r = failing([f'{BIN}/pack', '8', 'roots', f'{tmp}/q1.bin', f'{tmp}/selfpool',
                     f'{tmp}/meets.jsonl'])
    check('pack roots refuses a pool whose members meet the query',
          ok and 'companion pool' in r.stderr,
          (r.stderr.strip().splitlines() or ['no message'])[-1])

    # (d) resume after a real kill, and after corruption.  A shard that was in
    # flight leaves a torn manifest line; a shard whose payload has been damaged
    # is recorded as done but no longer matches its digest.  Both must be redone
    # -- and nothing else must be.
    rdir, rman = f'{tmp}/resume/shards', f'{tmp}/resume/man.jsonl'
    os.makedirs(f'{tmp}/resume', exist_ok=True)
    # The first 600 shards: enough for a kill to land in the middle of the work,
    # and their indices are 0..599, so they are a shard universe in their own
    # right for the audit to be given.
    rshards = f'{tmp}/resume/shards.txt'
    with open(rshards, 'w') as f:
        f.writelines(open(shardfile).readlines()[:600])
    killed = subprocess.Popen([f'{BIN}/enum', '8', 'shards', rshards, rdir, rman],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    while not (os.path.exists(rman) and os.path.getsize(rman) > 8000):
        if killed.poll() is not None:
            break
        time.sleep(0.02)
    killed.kill(); killed.wait()
    part_files = [f for _, _, fs in os.walk(rdir) for f in fs if f.endswith('.part')]
    banked = len(open(rman).read().splitlines())
    # A torn final line: exactly what a kill during a write leaves behind.
    with open(rman, 'a') as f:
        f.write('{"idx":99999,"count":3,"nod')
    # And one completed payload damaged after the fact.
    done_before = [json.loads(l) for l in open(rman).read().splitlines()
                   if l.endswith('}')]
    nonempty = [r for r in done_before if r['count'] > 0]
    victim_idx = nonempty[len(nonempty) // 2]['idx']
    vpath = os.path.join(rdir, f'{victim_idx // 1000:03d}', f's{victim_idx:05d}.bin')
    vb = bytearray(open(vpath, 'rb').read())
    vb[0] ^= 1
    open(vpath, 'wb').write(bytes(vb))

    r2 = run([f'{BIN}/enum', '8', 'shards', rshards, rdir, rman])
    regen = [l for l in r2.stderr.splitlines() if 'regenerating it' in l]
    torn = [l for l in r2.stderr.splitlines() if 'torn final line' in l]
    check('a killed sweep resumes: the torn final manifest line is discarded, and '
          'exactly the shard whose payload was damaged is regenerated',
          r2.returncode == 0 and len(torn) == 1 and len(regen) == 1
          and f'shard {victim_idx}:' in regen[0]
          and 'digest' in regen[0],
          f'killed after {banked} shards, {len(part_files)} .part files left; '
          f'{(regen or ["nothing regenerated"])[0].strip()}')

    left_part = [f for _, _, fs in os.walk(rdir) for f in fs if f.endswith('.part')]
    rcat = f'{tmp}/n8_resumed.bin'
    ra = run([PY, f'{ROOT}/src/catalogue.py', 'audit', '8', rman, rdir,
              '--shards', rshards])
    run([PY, f'{ROOT}/src/catalogue.py', 'pack', '8', rman, rdir, rcat,
         '--shards', rshards])
    # The same 600 shards, enumerated once through with nothing interrupted.
    cdir, cman = f'{tmp}/resume/clean', f'{tmp}/resume/clean.jsonl'
    run([f'{BIN}/enum', '8', 'shards', rshards, cdir, cman])
    ccat = f'{tmp}/n8_clean.bin'
    run([PY, f'{ROOT}/src/catalogue.py', 'pack', '8', cman, cdir, ccat,
         '--shards', rshards])
    check('the interrupted-then-resumed sweep reproduces a single-pass sweep of '
          'the same shards byte for byte, and leaves no .part files behind',
          ra.returncode == 0 and not left_part and os.path.exists(rcat)
          and sha256(rcat) == sha256(ccat),
          f'{len(left_part)} .part files left, {os.path.getsize(rcat)} bytes')

    # (e) the audit must reject what it is there to catch.
    audits = []
    empty_man = f'{tmp}/empty_man.jsonl'
    open(empty_man, 'w').close()
    audits.append(('an empty manifest', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'audit', '8', empty_man, f'{tmp}/shards',
         '--shards', shardfile])[0]))
    short_shards = f'{tmp}/short_shards.txt'
    with open(short_shards, 'w') as f:
        f.writelines(open(shardfile).readlines()[:-1])
    audits.append(('a manifest with a shard outside the universe', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'audit', '8', f'{tmp}/manifest.jsonl',
         f'{tmp}/shards', '--shards', short_shards])[0]))
    bad_man = f'{tmp}/bad_man.jsonl'
    with open(bad_man, 'w') as f:
        lines = open(f'{tmp}/manifest.jsonl').read().splitlines()
        f.write('\n'.join(lines[:5] + ['{"idx":3,"bogus":1}'] + lines[5:]) + '\n')
    audits.append(('a malformed line in the middle of a manifest', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'audit', '8', bad_man, f'{tmp}/shards',
         '--shards', shardfile])[0]))
    nodig_man = f'{tmp}/nodigest_man.jsonl'
    with open(nodig_man, 'w') as f:
        for line in open(f'{tmp}/manifest.jsonl'):
            r = json.loads(line)
            r.pop('digest', None)
            f.write(json.dumps(r) + '\n')
    audits.append(('a completion record with no payload digest', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'audit', '8', nodig_man, f'{tmp}/shards',
         '--shards', shardfile])[0]))
    audits.append(('an audit with no universe to compare against', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'audit', '8', f'{tmp}/manifest.jsonl',
         f'{tmp}/shards'])[0]))
    unsorted_dir = f'{tmp}/unsorted'
    shutil.copytree(f'{tmp}/shards', unsorted_dir)
    first = json.loads(open(f'{tmp}/manifest.jsonl').readline())
    up = os.path.join(unsorted_dir, f'{first["idx"] // 1000:03d}', f's{first["idx"]:05d}.bin')
    raw = open(up, 'rb').read()
    open(up, 'wb').write(raw[-64:] + raw[:-64])          # rotate: same set, wrong order
    audits.append(('records out of lexicographic order inside a shard', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'audit', '8', f'{tmp}/manifest.jsonl',
         unsorted_dir, '--shards', shardfile])[0]))
    check('the audit rejects an empty manifest, a shard outside the universe, a '
          'malformed line, a missing digest, a missing universe, and records out '
          'of order',
          all(ok for _, ok in audits),
          ', '.join(n for n, ok in audits if not ok) or 'all six rejected')

    # (f) the per-worker report validator, which is what stops a killed worker
    # from being read as a silence.
    rep = f'{tmp}/reports'
    os.makedirs(rep, exist_ok=True)
    good = dict(part=0, parts=2, queries=2, query_ids=[0, 2], disagreements=0,
                complete=True)
    open(f'{rep}/poolchk_0.json', 'w').write(json.dumps(good))
    open(f'{rep}/poolchk_1.json', 'w').write(json.dumps(
        dict(good, part=1, query_ids=[1, 3])))
    full_ok, _ = failing([PY, f'{ROOT}/src/catalogue.py', 'reports', 'pools', rep,
                          'poolchk_*.json', '--parts', '2', '--expect', '4'])
    open(f'{rep}/poolchk_1.json', 'w').close()          # the killed worker
    empty_ok, er = failing([PY, f'{ROOT}/src/catalogue.py', 'reports', 'pools', rep,
                            'poolchk_*.json', '--parts', '2', '--expect', '4'])
    os.remove(f'{rep}/poolchk_1.json')                  # or no file at all
    gone_ok, _ = failing([PY, f'{ROOT}/src/catalogue.py', 'reports', 'pools', rep,
                          'poolchk_*.json', '--parts', '2', '--expect', '4'])
    open(f'{rep}/poolchk_1.json', 'w').write(json.dumps(
        dict(good, part=1, query_ids=[1])))             # a worker that stopped short
    short_ok, _ = failing([PY, f'{ROOT}/src/catalogue.py', 'reports', 'pools', rep,
                           'poolchk_*.json', '--parts', '2', '--expect', '4'])
    check('the report validator accepts two complete reports and rejects an empty '
          'one, a missing one, and one that covers fewer queries than it should',
          not full_ok and empty_ok and gone_ok and short_ok,
          'complete set accepted, all three failures caught')

    # (g) the result validator must reject results that do not support the
    # theorem.  Order 8 is the case where a cube EXISTS, so the same ledger the
    # order-8 controls are happy with must be refused as a proof of
    # nonexistence: 7-cliques and covers are exactly what it is looking for.
    rawres = f'{tmp}/n8_results_raw.jsonl'
    with open(rawres, 'w') as f:
        for k in range(W):
            f.write(open(f'{tmp}/roots_{k}.jsonl').read())
    run([PY, f'{ROOT}/src/catalogue.py', 'ledger', rawres, f'{tmp}/n8_ledger.jsonl'])
    v8_ok, v8 = failing([PY, f'{ROOT}/src/catalogue.py', 'validate', '8',
                         f'{tmp}/n8_ledger.jsonl', '--queries', '13056'])
    short_ledger = f'{tmp}/short_ledger.jsonl'
    with open(short_ledger, 'w') as f:
        f.writelines(open(f'{ROOT}/data/n9_root_results.jsonl').read().splitlines(True)[:-1])
    miss_ok, _ = failing([PY, f'{ROOT}/src/catalogue.py', 'validate', '9',
                          short_ledger, '--queries', '2049'])
    good_ok, gr = failing([PY, f'{ROOT}/src/catalogue.py', 'validate', '9',
                           f'{ROOT}/data/n9_root_results.jsonl', '--queries', '2049'])
    check('the result validator passes the shipped order-9 ledger, rejects it with '
          'one orbit missing, and rejects the order-8 ledger, where a cube exists',
          v8_ok and miss_ok and not good_ok,
          f'order 8: {json.loads(v8.stdout)["problem_count"]} problems')

    # (h) the checksum manifest is required by name, not by whatever happens to
    # be in both files.
    ck = []
    good_ck = f'{tmp}/ck_good.txt'
    open(good_ck, 'w').write(f'{"0" * 64}  a.bin\n{"1" * 64}  b.bin\n')
    ck.append(('the matching pair', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'checksums', good_ck, good_ck,
         '--require', 'a.bin', 'b.bin'])[0]))
    ck.append(('a required name that is not listed', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'checksums', good_ck, good_ck,
         '--require', 'a.bin', 'b.bin', 'c.bin'])[0]))
    dup_ck = f'{tmp}/ck_dup.txt'
    open(dup_ck, 'w').write(f'{"0" * 64}  a.bin\n{"2" * 64}  a.bin\n')
    ck.append(('a duplicated name', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'checksums', dup_ck, good_ck,
         '--require', 'a.bin'])[0]))
    junk_ck = f'{tmp}/ck_junk.txt'
    open(junk_ck, 'w').write('not-a-digest  a.bin\n')
    ck.append(('a malformed digest', failing(
        [PY, f'{ROOT}/src/catalogue.py', 'checksums', junk_ck, good_ck,
         '--require', 'a.bin'])[0]))
    check('the checksum comparison requires the exact manifest, and rejects a '
          'missing name, a duplicated name and a malformed digest',
          not ck[0][1] and all(ok for _, ok in ck[1:]),
          ', '.join(n for n, ok in ck[1:] if not ok) or 'all three rejected')

    # ---- 12. determinism -------------------------------------------------
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
