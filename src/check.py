#!/usr/bin/env python3
"""check.py -- independent checks on the data the C tools produce.

Nothing here shares code with the C programs.  In particular the main lines are
re-derived from the definition in numpy, so a record file can be checked
against what a support *is* rather than against what the enumerator thinks it
built.  numpy is the only dependency.

Subcommands (each has --help):

  lines N                     report the main-line count and cell degrees
  verify N FILE [--sample K]  check records are supports, by the definition
  a007016 N                   the admissible-row count in closed form
  setcmp N A.bin B.bin        compare two record files as sets
  pools N CAT QUERIES DIR     re-derive every companion pool and compare
  witness N FILE.json         check an explicit packing of pairwise disjoint
                              supports, from the definition
  sha256 FILE...              print `sha256  path` lines
"""
import argparse
import hashlib
import itertools
import json
import os
import sys

import numpy as np


# ---------------------------------------------------------------- main lines
def main_lines(n):
    """Every main line of [n]^3, once, as an (L, n) array of flat cell indices.

    A main line is indexed by t = 0..n-1 with each coordinate either a constant,
    or t, or n-1-t, and at least one coordinate non-constant.  Reversing t gives
    the same set, so each line is emitted once by requiring the first
    non-constant coordinate to be t.
    """
    out = []
    for kinds in itertools.product(('c', 't', 'r'), repeat=3):
        if all(k == 'c' for k in kinds):
            continue
        first = next(i for i, k in enumerate(kinds) if k != 'c')
        if kinds[first] != 't':
            continue
        consts = [i for i, k in enumerate(kinds) if k == 'c']
        for vals in itertools.product(range(n), repeat=len(consts)):
            cmap = dict(zip(consts, vals))
            line = []
            for t in range(n):
                x = [cmap[i] if kinds[i] == 'c' else (t if kinds[i] == 't' else n - 1 - t)
                     for i in range(3)]
                line.append((x[0] * n + x[1]) * n + x[2])
            out.append(line)
    return np.array(out, dtype=np.int64)


def load_records(path, n):
    w = n * n
    size = os.path.getsize(path)
    if size % w:
        sys.exit(f'{path}: size {size} is not a multiple of {w}')
    return np.fromfile(path, dtype=np.uint8).reshape(size // w, w)


def cells_of(recs, n):
    """(M, n*n) records -> (M, n*n) flat cell indices."""
    return np.arange(n * n, dtype=np.int64)[None, :] * n + recs.astype(np.int64)


# --------------------------------------------------------------- subcommands
def cmd_lines(a):
    L = main_lines(a.n)
    deg = np.bincount(L.ravel(), minlength=a.n ** 3)
    expect = 3 * a.n * a.n + 6 * a.n + 4
    print(f'n={a.n} lines={len(L)} expected={expect}')
    print(f'cells={a.n ** 3} incidences={L.size} degree min={deg.min()} max={deg.max()}')
    ok = len(L) == expect and all(len(set(row.tolist())) == a.n for row in L)
    print('every line has n distinct cells' if ok else 'MALFORMED LINES')
    return 0 if ok else 1


def cmd_verify(a):
    n = a.n
    recs = load_records(a.file, n)
    M = recs.shape[0]
    idx = np.arange(M)
    if a.sample and a.sample < M:
        idx = np.sort(np.random.default_rng(a.seed).choice(M, a.sample, replace=False))
    L = main_lines(n)
    bad_cells = bad_lines = 0
    for lo in range(0, len(idx), 4096):
        blk = cells_of(recs[idx[lo:lo + 4096]], n)
        # n*n distinct cells
        srt = np.sort(blk, axis=1)
        bad_cells += int((np.diff(srt, axis=1) == 0).any(axis=1).sum())
        # each main line met exactly once
        flat = np.zeros((blk.shape[0], n ** 3), dtype=np.int8)
        np.put_along_axis(flat, blk, 1, axis=1)
        hits = flat[:, L].sum(axis=2)          # (block, lines)
        bad_lines += int((hits != 1).any(axis=1).sum())
    print(json.dumps(dict(file=a.file, records=int(M), checked=int(len(idx)),
                          lines=int(len(L)),
                          records_without_n2_distinct_cells=int(bad_cells),
                          records_missing_a_line_exactly_once=int(bad_lines))))
    return 0 if (bad_cells == 0 and bad_lines == 0) else 1


def a007016(n):
    """Permutations of [n] with exactly one fixed and one reflected point.

    Inclusion-exclusion on the two boards A = {(j,j)} and B = {(j, n-1-j)}.
    Their union is floor(n/2) disjoint 4-cycles plus, for odd n, the shared
    centre cell; a selection may use that cell as an A-cell, a B-cell or both,
    so the completion factor must count DISTINCT cells -- the u coordinate.
    """
    from math import comb, factorial
    P = {(0, 0, 0): 1}

    def mul(p, q):
        r = {}
        for k, u in p.items():
            for l, v in q.items():
                m = (k[0] + l[0], k[1] + l[1], k[2] + l[2])
                r[m] = r.get(m, 0) + u * v
        return r

    cyc = {(0, 0, 0): 1, (1, 0, 1): 2, (0, 1, 1): 2, (2, 0, 2): 1, (0, 2, 2): 1}
    for _ in range(n // 2):
        P = mul(P, cyc)
    if n % 2:
        P = mul(P, {(0, 0, 0): 1, (1, 0, 1): 1, (0, 1, 1): 1, (1, 1, 1): 1})
    tot = 0
    for (x, y, u), m in P.items():
        if x >= 1 and y >= 1 and u <= n:
            tot += (-1) ** (x - 1 + y - 1) * comb(x, 1) * comb(y, 1) * m * factorial(n - u)
    return tot


def cmd_a007016(a):
    vals = [a007016(k) for k in range(1, a.n + 1)]
    print(json.dumps(dict(n=a.n, sequence=vals, value=vals[-1])))
    return 0


def cmd_setcmp(a):
    A = load_records(a.a, a.n)
    B = load_records(a.b, a.n)
    sa = {r.tobytes() for r in A}
    sb = {r.tobytes() for r in B}
    print(json.dumps(dict(a=a.a, b=a.b, a_records=int(A.shape[0]), b_records=int(B.shape[0]),
                          a_distinct=len(sa), b_distinct=len(sb),
                          only_in_a=len(sa - sb), only_in_b=len(sb - sa),
                          set_equal=sa == sb)))
    return 0 if sa == sb else 1


def masks_of(recs, n, chunk=1 << 18):
    """(M, n*n) records -> (M, W) uint64 cell bitmasks, chunked to bound memory."""
    W = (n ** 3 + 63) // 64
    m = np.zeros((recs.shape[0], W), dtype=np.uint64)
    for lo in range(0, recs.shape[0], chunk):
        c = cells_of(recs[lo:lo + chunk], n)
        word = (c // 64).astype(np.int64)
        bit = np.left_shift(np.uint64(1), (c % 64).astype(np.uint64))
        for w in range(W):
            m[lo:lo + chunk, w] = np.bitwise_or.reduce(
                np.where(word == w, bit, np.uint64(0)), axis=1)
    return m


def cmd_pools(a):
    n = a.n
    cat = load_records(a.catalogue, n)
    qry = load_records(a.queries, n)
    cm = masks_of(cat, n)
    qms = masks_of(qry, n)
    bad = 0
    sizes = []
    for j in range(qry.shape[0]):
        qm = qms[j]
        sel = np.all((cm & qm) == 0, axis=1)
        want = cat[sel]
        path = os.path.join(a.pooldir, f'pool_{j}.bin')
        have = load_records(path, n)
        sizes.append(int(want.shape[0]))
        if have.shape != want.shape or not np.array_equal(have, want):
            hs = {r.tobytes() for r in have}
            ws = {r.tobytes() for r in want}
            print(f'pool {j}: disk {have.shape[0]} rescan {want.shape[0]} '
                  f'only_on_disk {len(hs - ws)} only_in_rescan {len(ws - hs)}')
            bad += 1
        if (j + 1) % 250 == 0:
            print(f'  {j + 1}/{qry.shape[0]}', flush=True)
    print(json.dumps(dict(catalogue=int(cat.shape[0]), queries=int(qry.shape[0]),
                          disagreements=bad, pool_min=min(sizes), pool_max=max(sizes),
                          pool_mean=round(sum(sizes) / len(sizes), 1))))
    return 0 if bad == 0 else 1


def cmd_witness(a):
    n = a.n
    w = json.load(open(a.file))
    squares = w['squares']
    L = main_lines(n)
    recs = np.array([[int(v) for v in s] for s in squares], dtype=np.uint8)
    cells = cells_of(recs, n)
    problems = []
    for i, row in enumerate(cells):
        flat = np.zeros(n ** 3, dtype=np.int8)
        flat[row] = 1
        if flat.sum() != n * n:
            problems.append(f'square {i} has repeated cells')
        hits = flat[L].sum(axis=1)
        if not np.all(hits == 1):
            problems.append(f'square {i} misses {int((hits != 1).sum())} lines exactly once')
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            if set(cells[i].tolist()) & set(cells[j].tolist()):
                problems.append(f'squares {i} and {j} are not disjoint')
    print(json.dumps(dict(file=a.file, squares=len(squares), cells_each=n * n,
                          lines=int(len(L)), pairwise_disjoint=not problems,
                          problems=problems)))
    return 0 if not problems else 1


def cmd_sha256(a):
    for path in a.files:
        h = hashlib.sha256()
        with open(path, 'rb') as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                h.update(b)
        print(f'{h.hexdigest()}  {path}')
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('lines', help='main-line count and cell degrees')
    p.add_argument('n', type=int); p.set_defaults(fn=cmd_lines)

    p = sub.add_parser('verify', help='check records are supports, by the definition')
    p.add_argument('n', type=int); p.add_argument('file')
    p.add_argument('--sample', type=int, default=0, help='check K random records, not all')
    p.add_argument('--seed', type=int, default=20260906)
    p.set_defaults(fn=cmd_verify)

    p = sub.add_parser('a007016', help='admissible row-0 count in closed form')
    p.add_argument('n', type=int); p.set_defaults(fn=cmd_a007016)

    p = sub.add_parser('setcmp', help='compare two record files as sets')
    p.add_argument('n', type=int); p.add_argument('a'); p.add_argument('b')
    p.set_defaults(fn=cmd_setcmp)

    p = sub.add_parser('pools', help='re-derive every companion pool and compare')
    p.add_argument('n', type=int); p.add_argument('catalogue')
    p.add_argument('queries'); p.add_argument('pooldir')
    p.set_defaults(fn=cmd_pools)

    p = sub.add_parser('witness', help='check an explicit packing from the definition')
    p.add_argument('n', type=int); p.add_argument('file')
    p.set_defaults(fn=cmd_witness)

    p = sub.add_parser('sha256', help='print sha256 checksums')
    p.add_argument('files', nargs='+'); p.set_defaults(fn=cmd_sha256)

    a = ap.parse_args()
    sys.exit(a.fn(a))


if __name__ == '__main__':
    main()
