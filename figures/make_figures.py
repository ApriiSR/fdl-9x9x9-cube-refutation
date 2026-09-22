#!/usr/bin/env python3
"""Draw the README figures, each in a light and a dark variant.

    python3 figures/make_figures.py extract   # needs work/pools/ from a run
    python3 figures/make_figures.py draw      # the static SVGs
    python3 figures/make_figures.py gif       # the two rotating GIFs; needs
                                              # Chrome, ImageMagick and ffmpeg

`extract` pulls the few supports the figures show out of a finished run and
writes them to figures/figure_data.json, which is committed, so `draw` works
from a fresh checkout.  $CHROME overrides the path to Chrome.  A record is 81 bytes: byte 9*x0 + x1 is x2.
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, 'figure_data.json')
N = 9

THEMES = {
    'light': dict(line='#6e7781', faint='#d0d7de', text='#1f2328', bg='#ffffff',
                  plane='#9a6700'),
    'dark': dict(line='#8b949e', faint='#30363d', text='#e6edf3', bg='#0d1117',
                 plane='#d29922'),
}
RED = '#cf222e'
COMPANIONS = ['#0969da', '#1a7f37', '#bc4c00', '#8250df']   # never red

SERIF = "font-family=\"Georgia, 'Times New Roman', serif\""
SANS = 'font-family="-apple-system, \'Segoe UI\', Helvetica, Arial, sans-serif"'


# ---- data -----------------------------------------------------------------

def extract():
    import numpy as np
    pk = json.load(open(os.path.join(ROOT, 'data', 'exceptional_packing.json')))
    q = pk['query']
    pool = np.fromfile(os.path.join(ROOT, 'work', 'pools', f'pool_{q}.bin'),
                       np.uint8).reshape(-1, 81)
    cells = np.arange(81) * N
    mask = np.zeros((len(pool), N ** 3), np.int32)
    for k, rec in enumerate(pool):
        mask[k, cells + rec] = 1
    adj = (mask @ mask.T) == 0
    for k in range(len(pool)):
        adj[k, k] = False
    # the connected component of the pool graph containing the packing
    comp, stack = set(pk['members'][1:]), list(pk['members'][1:])
    while stack:
        v = stack.pop()
        for w in map(int, adj[v].nonzero()[0]):
            if w not in comp:
                comp.add(w)
                stack.append(w)
    comp = sorted(comp)
    deg = adj.sum(1)
    isolated = int((deg == 0).nonzero()[0][0])
    pair = next([v, int(adj[v].nonzero()[0][0])] for v in map(int, (deg == 1).nonzero()[0])
                if deg[adj[v].nonzero()[0][0]] == 1)
    out = dict(
        query=q,
        root=pk['squares'][0],
        pool_size=int(len(pool)),
        pool_edges=int(adj.sum() // 2),
        pool_isolated=int((adj.sum(1) == 0).sum()),
        packing=pk['members'][1:],
        component={str(v): [int(b) for b in pool[v]] for v in comp},
        edges=[[a, b] for a in comp for b in comp if a < b and adj[a, b]],
        # a component that is a single edge, and a vertex with no edges
        pair=pair,
        isolated=isolated,
        others={str(v): [int(b) for b in pool[v]] for v in pair + [isolated]},
        single_edges=int(sum(1 for v in range(len(pool)) if deg[v] == 1
                             and deg[adj[v].nonzero()[0][0]] == 1) // 2),
    )
    with open(DATA, 'w') as f:
        json.dump(out, f, indent=1)
    print(f'wrote {DATA}: component of {len(comp)} supports, '
          f'{len(out["edges"])} edges')


def cells_of(rec):
    return [(a, b, rec[N * a + b]) for a in range(N) for b in range(N)]


def rev(t):
    return N - 1 - t


def apply(rec, pi, sigma):
    """g(x)_i = sigma_i(x_{pi(i)}), as in Lemma 3."""
    R = [rev(t) for t in range(N)]
    for sg in sigma:
        assert all(sg[R[t]] == R[sg[t]] for t in range(N))          # in C(r)
        assert sg == sigma[0] or sg == [sigma[0][R[t]] for t in range(N)]
    out = [None] * 81
    for x in cells_of(rec):
        y = [sigma[i][x[pi[i]]] for i in range(3)]
        assert out[N * y[0] + y[1]] is None
        out[N * y[0] + y[1]] = y[2]
    return out


# ---- the 3D view ----------------------------------------------------------

def shade(color, f):
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    return '#%02x%02x%02x' % tuple(round(v * f) for v in (r, g, b))


def projector(u, az=0.52, el=0.42):
    """Orthographic camera above the cube: x1 east, x2 north, x0 up."""
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)
    h = (N - 1) / 2

    def P(x0, x1, x2):
        X, Y, Z = x1 - h, x2 - h, x0 - h
        d = X * sa + Y * ca                 # distance away from the viewer
        return ((X * ca - Y * sa) * u, -(Z * ce + d * se) * u, d * ce - Z * se)
    return P


def cube(layers, u, theme, axes=False, plane=None, s=0.62, az=0.52):
    """One cube.  layers: [(color, record, opacity)].  Returns (svg, box); the
    box does not depend on az, so frames of a rotation line up."""
    T = THEMES[theme]
    P = projector(u, az)
    lo, hi = -0.5, N - 0.5
    out = []
    if plane is not None:     # the plane x0 = plane, as a tinted slab
        pts = [P(plane, a, b)[:2] for a, b in ((lo, lo), (hi, lo), (hi, hi), (lo, hi))]
        out.append('<polygon points="%s" fill="%s" fill-opacity="0.16" stroke="%s" '
                   'stroke-width="1" stroke-dasharray="3 2"/>'
                   % (' '.join('%.1f,%.1f' % p for p in pts), T['plane'], T['plane']))
    for a in (lo, hi):
        for b in (lo, hi):
            for p, q in (((a, b, lo), (a, b, hi)), ((a, lo, b), (a, hi, b)),
                         ((lo, a, b), (hi, a, b))):
                A, B = P(*p), P(*q)
                out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                           'stroke-width="%.2f"/>' % (A[0], A[1], B[0], B[1], T['line'],
                                                      0.6 + u / 20))
    items = [(x, col, op) for col, rec, op in layers for x in cells_of(rec)]
    items.sort(key=lambda it: -P(*it[0])[2])
    h = s / 2
    for (x0, x1, x2), col, op in items:
        def V(d0, d1, d2):
            return P(x0 + d0, x1 + d1, x2 + d2)[:2]
        g = ''.join('<polygon points="%s" fill="%s"/>'
                    % (' '.join('%.1f,%.1f' % p for p in f), shade(col, k))
                    for f, k in faces(V, h, az))
        out.append(g if op == 1 else '<g opacity="%.2f">%s</g>' % (op, g))
    R = math.sqrt(2) * N / 2 * u
    top = (N / 2 * math.cos(0.42) + math.sqrt(2) * N / 2 * math.sin(0.42)) * u
    box = [-R - 3, -top - 3, R + 3, top + 3]
    if axes:
        # a small tripod below and left of the cube, rotating with it
        ox, oy, L = -R + 30, top + 26, 26
        O = P(0, 0, 0)
        for i, name in enumerate('012'):
            e = [0, 0, 0]
            e[i] = 1
            E = P(*e)
            dx, dy = E[0] - O[0], E[1] - O[1]
            m = math.hypot(dx, dy)
            if m < 1e-6:
                continue
            dx, dy = dx / m, dy / m
            k = min(1, m / u + 0.25)            # foreshortened axes look shorter
            tx, ty = ox + dx * L * k, oy + dy * L * k
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="1.3" marker-end="url(#arrow-%s)"/>'
                       % (ox, oy, tx, ty, T['text'], theme))
            out.append(var(ox + dx * (L * k + 11), oy + dy * (L * k + 11) + 5, name, theme))
        box[3] = oy + L + 24
        box[0] = min(box[0], ox - L - 24)
    return ''.join(out), tuple(box)


def faces(V, h, az):
    """The three faces of a small cube the camera sees, with their shading."""
    sx1 = -1 if math.sin(az) > 0 else 1        # which x1 face points at the camera
    sx2 = -1 if math.cos(az) > 0 else 1
    return [([V(h, -h, -h), V(h, h, -h), V(h, h, h), V(h, -h, h)], 1.0),
            ([V(-h, -h, sx2 * h), V(h, -h, sx2 * h), V(h, h, sx2 * h), V(-h, h, sx2 * h)], 0.8),
            ([V(-h, sx1 * h, -h), V(h, sx1 * h, -h), V(h, sx1 * h, h), V(-h, sx1 * h, h)], 0.62)]


def var(x, y, sub, theme, size=15):
    """An italic x with a numeric subscript, without relying on Unicode subscripts."""
    return ('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="%d" fill="%s" %s>'
            '<tspan font-style="italic">x</tspan><tspan dy="4" font-size="%d">%s</tspan></text>'
            % (x, y, size, THEMES[theme]['text'], SERIF, round(size * 0.7), sub))


def defs(theme):
    c = THEMES[theme]['text']
    return ('<defs><marker id="arrow-%s" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            '<path d="M0 0L10 5L0 10z" fill="%s"/></marker></defs>' % (theme, c))


def svg(body, box, theme, title):
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.0f %.0f %.0f %.0f" '
            'width="%.0f" height="%.0f" role="img"><title>%s</title>%s%s</svg>\n'
            % (x0, y0, w, h, w, h, title, defs(theme), body))


def place(frag, box, dx, dy):
    return '<g transform="translate(%.1f %.1f)">%s</g>' % (dx - box[0], dy - box[1], frag)


def label(x, y, text, theme, size=14, anchor='middle', italic=False):
    return ('<text x="%.1f" y="%.1f" text-anchor="%s" font-size="%d"%s fill="%s" %s>%s</text>'
            % (x, y, anchor, size, ' font-style="italic"' if italic else '',
               THEMES[theme]['text'], SERIF if italic else SANS, text))


# ---- the figures ----------------------------------------------------------

def fig_root(d, theme, az=0.52):
    frag, box = cube([(RED, d['root'], 1)], 17, theme, axes=True, az=az)
    return svg(frag, box, theme, 'A support of [9]^3 containing the center cell')


def fig_pair(d, theme, az=0.52):
    blue = d['component'][str(d['packing'][0])]
    frag, box = cube([(RED, d['root'], 1), (COMPANIONS[0], blue, 1)], 17, theme,
                     axes=True, az=az)
    return svg(frag, box, theme, 'Two disjoint supports in one cube')


def fig_orbit(d, theme):
    idt = list(range(N))
    R = [rev(t) for t in range(N)]
    tau = idt[:]
    tau[0], tau[1], tau[N - 2], tau[N - 1] = 1, 0, N - 1, N - 2
    T = d['root']
    panels = [
        (T, ['the root T']),
        (apply(T, (1, 2, 0), (idt, idt, idt)), ['cycle the axes', 'x₀ → x₂ → x₁ → x₀']),
        (apply(T, (0, 1, 2), (R, idt, idt)), ['reverse x₀', '(turn the cube over)']),
        (apply(T, (0, 1, 2), (tau, tau, tau)), ['relabel values', '0 ↔ 1 and 7 ↔ 8']),
    ]
    out, x = [], 0
    for rec, lines in panels:
        frag, box = cube([(RED, rec, 1)], 11, theme)
        w, h = box[2] - box[0], box[3] - box[1]
        out.append(place(frag, box, x, 0))
        for k, t in enumerate(lines):
            out.append(label(x + w / 2, h + 20 + 18 * k, t, theme, 14))
        x += w + 26
    return svg(''.join(out), (-4, -4, x - 22, h + 50), theme,
               'The root and three of its images under the symmetry group')


def greedy_colors(nodes, edges):
    nbr = {v: set() for v in nodes}
    for a, b in edges:
        nbr[a].add(b)
        nbr[b].add(a)
    col = {}
    for v in nodes:
        used = {col[w] for w in nbr[v] if w in col}
        col[v] = min(c for c in range(len(COMPANIONS)) if c not in used)
    return col


def fig_graph(d, theme):
    Th = THEMES[theme]
    A = d['packing']                                   # one 4-clique
    edges = [tuple(e) for e in d['edges']]
    adj = {}
    for a, b in edges:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    # the rest of the component: each clique member's one outside neighbor
    partner = {a: next(iter(adj[a] - set(A))) for a in A}
    # left square TL TR / BL BR, right square mirrored so two matching edges
    # are short and the other two arc over and under
    TL, TR, BL, BR = A
    pos, gap, side = {}, 190, 150
    for v, (i, j) in zip((TL, TR, BL, BR), ((0, 0), (1, 0), (0, 1), (1, 1))):
        pos[v] = (i * side, j * side)
    for v, (i, j) in zip((partner[TR], partner[TL], partner[BR], partner[BL]),
                         ((0, 0), (1, 0), (0, 1), (1, 1))):
        pos[v] = (side + gap + i * side, j * side)
    nodes = list(A) + [partner[a] for a in A]
    col = greedy_colors(nodes, edges)
    assert all(col[a] != col[b] for a, b in edges)
    assert col[A[0]] == 0                               # blue, as in the pair figure
    out = []
    for a, b in edges:
        (x0, y0), (x1, y1) = pos[a], pos[b]
        clique = a in A and b in A or a not in A and b not in A
        w = 2.2 if clique else 1.4
        if y0 == y1 and abs(x1 - x0) > side + gap:     # the long matching edges
            # control point far enough out that the curve clears the boxes
            # it passes over (a quadratic reaches half its control offset)
            bend = -170 if y0 == 0 else 170
            out.append('<path d="M%.1f %.1f Q%.1f %.1f %.1f %.1f" fill="none" '
                       'stroke="%s" stroke-width="%.1f"/>'
                       % (x0, y0, (x0 + x1) / 2, y0 + bend, x1, y1, Th['line'], w))
        else:
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="%.1f"/>' % (x0, y0, x1, y1, Th['line'], w))
    # below: a component that is one edge, and a vertex with no edges
    sup = dict(d['component'], **d['others'])
    a, b = d['pair']
    row = side + 190
    pos[a], pos[b] = (0, row), (side, row)             # under the left square
    pos[d['isolated']] = (side + gap + side / 2, row)  # under the right one
    col[a], col[b], col[d['isolated']] = 0, 1, 0
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
               'stroke-width="2.2"/>' % (pos[a][0], row, pos[b][0], row, Th['line']))
    nodes = nodes + [a, b, d['isolated']]
    r = 44
    for v in nodes:
        frag, box = cube([(COMPANIONS[col[v]], sup[str(v)], 1)], 4.2, theme)
        x, y = pos[v]
        out.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" rx="10" fill="%s" '
                   'stroke="%s" stroke-width="1.2"/>'
                   % (x - r, y - r, 2 * r, 2 * r, Th['bg'], Th['line']))
        w, h = box[2] - box[0], box[3] - box[1]
        out.append(place(frag, box, x - w / 2, y - h / 2))
    return svg(''.join(out), (-r - 6, -r - 48, 2 * side + gap + r + 6, row + r + 8),
               theme, 'The component of the companion graph containing its 4-cliques')


def matrix_projector(u, az=0.36, el=0.30):
    """Camera above and to the right, for the Latin-square figures: x0 runs top
    to bottom and x1 left to right, as in a matrix, and x2 runs away from the
    viewer, so the front face x2 = 0 is where the Latin square sits."""
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)
    h = (N - 1) / 2

    def P(x0, x1, x2):
        X, U, D = x1 - h, h - x0, x2 - h
        Xp, Dp = X * ca + D * sa, -X * sa + D * ca
        return Xp * u, -(U * ce + Dp * se) * u, Dp * ce - U * se
    return P


def matrix_cube(rec, u, theme, opacity=None, segments=(), plane=False, s=0.62,
                seg_width=1.6, numbers=False):
    """The support rec in the matrix view.  opacity(x0) fades cells by row;
    segments lists the rows whose cells get a line back to the front face."""
    Th = THEMES[theme]
    P = matrix_projector(u)
    lo, hi = -0.5, N - 0.5

    def poly(pts, attrs):
        return '<polygon points="%s" %s/>' % (' '.join('%.1f,%.1f' % q[:2] for q in pts), attrs)
    out = []
    if plane:
        out.append(poly([P(lo, a, b) for a, b in ((lo, lo), (hi, lo), (hi, hi), (lo, hi))],
                        'fill="%s" fill-opacity="0.16" stroke="%s" stroke-width="1" '
                        'stroke-dasharray="3 2"' % (Th['plane'], Th['plane'])))
    for a in (lo, hi):
        for b in (lo, hi):
            for p0, p1 in (((a, b, lo), (a, b, hi)), ((a, lo, b), (a, hi, b)),
                           ((lo, a, b), (hi, a, b))):
                A, B = P(*p0), P(*p1)
                out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                           'stroke-width="1.35"/>' % (A[0], A[1], B[0], B[1], Th['line']))
    for i in range(1, N):     # the front face as a faint 9 x 9 grid
        for p0, p1 in (((i - 0.5, lo, lo), (i - 0.5, hi, lo)), ((lo, i - 0.5, lo), (hi, i - 0.5, lo))):
            A, B = P(*p0), P(*p1)
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="0.7"/>' % (A[0], A[1], B[0], B[1], Th['faint']))
    h = s / 2
    items = []
    for x0, x1, x2 in cells_of(rec):
        op = opacity(x0) if opacity else 1

        def V(d0, d1, d2):
            return P(x0 + d0, x1 + d1, x2 + d2)
        g = ''.join(poly(f, 'fill="%s"' % shade(RED, k)) for f, k in (
            ([V(-h, -h, -h), V(-h, h, -h), V(-h, h, h), V(-h, -h, h)], 1.0),    # top
            ([V(-h, -h, -h), V(h, -h, -h), V(h, h, -h), V(-h, h, -h)], 0.8),    # front
            ([V(-h, h, -h), V(h, h, -h), V(h, h, h), V(-h, h, h)], 0.62)))      # right
        items.append((P(x0, x1, x2)[2], g if op == 1 else '<g opacity="%.2f">%s</g>' % (op, g)))
        if x0 in segments:
            # one piece per cell the segment crosses, so it sorts with the cubes
            ends = [lo] + [k + 0.5 for k in range(x2)]
            ends[-1:] = [ends[-1]] if x2 == 0 else ends[-1:]
            stops = [k + 0.5 for k in range(x2 - 1)] + [x2 - h]
            for z0, z1 in zip([lo] + stops[:-1], stops):
                A, B = P(x0, x1, z0), P(x0, x1, z1)
                items.append((P(x0, x1, (z0 + z1) / 2)[2],
                              '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                              'stroke-width="%.1f" stroke-linecap="round"/>'
                              % (A[0], A[1], B[0], B[1], Th['text'], seg_width)))
    items.sort(key=lambda t: -t[0])
    out += [g for _, g in items]
    for x0, x1, x2 in cells_of(rec):
        A = P(x0, x1, lo)
        if numbers:      # the Latin square written on the front face
            out.append('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="%.1f" '
                       'fill="%s" stroke="%s" stroke-width="2.5" paint-order="stroke" %s>%d</text>'
                       % (A[0], A[1] + u * 0.2, u * 0.55, Th['text'], Th['bg'], SANS, x2))
        elif x0 in segments:
            out.append('<circle cx="%.1f" cy="%.1f" r="2.2" fill="%s"/>'
                       % (A[0], A[1], Th['text']))
    # axes: x0 down the left of the front face, x1 along its bottom, x2 away
    # along the bottom right edge, each labelled beside its middle
    g = 1.0
    arrows = (((lo, lo - g, lo), (hi, lo - g, lo), '0', 0.5, (-14, 5)),
              ((hi + g, lo, lo), (hi + g, hi, lo), '1', 0.5, (0, 20)),
              ((hi + g, hi + g * 0.6, lo), (hi + g, hi + g * 0.6, hi), '2', 0.5, (13, 14)))
    for p0, p1, name, t, (dx, dy) in arrows:
        A, B = P(*p0), P(*p1)
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="1.2" marker-end="url(#arrow-%s)"/>'
                   % (A[0], A[1], B[0], B[1], Th['text'], theme))
        out.append(var(A[0] + t * (B[0] - A[0]) + dx, A[1] + t * (B[1] - A[1]) + dy,
                       name, theme))
    pts = [P(a, b, c) for a in (lo, hi + 2.2) for b in (lo - 1.8, hi + 1.6) for c in (lo, hi)]
    xs, ys = [q[0] for q in pts], [q[1] for q in pts]
    return ''.join(out), (min(xs) - 6, min(ys) - 6, max(xs) + 6, max(ys) + 6)


def latin_panel(T, X, Y, c, theme, rows=(), rings=()):
    """The Latin square of T as a grid of numbers, some rows shaded."""
    Th = THEMES[theme]
    out = []
    for a in rows:
        out.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" fill="%s" '
                   'fill-opacity="0.18"/>' % (X, Y + a * c, N * c, c, Th['plane']))
    for i in range(N + 1):
        out.append('<path d="M%.1f %.1fv%dM%.1f %.1fh%d" stroke="%s" stroke-width="%.1f"/>'
                   % (X + i * c, Y, N * c, X, Y + i * c, N * c,
                      Th['line'] if i in (0, N) else Th['faint'], 1.2 if i in (0, N) else 1))
    for a in range(N):
        for b in range(N):
            out.append('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="13" '
                       'fill="%s" %s>%d</text>' % (X + b * c + c / 2, Y + a * c + c / 2 + 4.5,
                                                   Th['text'], SANS, T[N * a + b]))
    for t, dash in rings:
        out.append('<circle cx="%.1f" cy="%.1f" r="10" fill="none" stroke="%s" '
                   'stroke-width="1.6"%s/>' % (X + t * c + c / 2, Y + c / 2, Th['text'], dash))
    return ''.join(out)


def fig_depth(d, theme):
    T = d['root']
    frag, box = matrix_cube(T, 21, theme, segments=range(N), seg_width=0.9, numbers=True)
    w0, ht = box[2] - box[0], box[3] - box[1]
    c = 24
    X, Y = w0 + 36, (ht - N * c) / 2
    out = [place(frag, box, 0, 0), latin_panel(T, X, Y, c, theme)]
    out.append(label(w0 / 2, ht + 18, 'the support, seen from the front', theme, 14))
    out.append(label(X + N * c / 2, ht + 18, 'its Latin square L(x₀, x₁) = x₂', theme, 14))
    return svg(''.join(out), (-4, -4, X + N * c + 6, ht + 28), theme,
               'A support as a Latin square read from the front face')


def fig_lemma2(d, theme):
    Th = THEMES[theme]
    T = d['root']
    p = T[:N]
    fixed = [t for t in range(N) if p[t] == t]
    refl = [t for t in range(N) if p[t] == rev(t)]
    assert len(fixed) == 1 and len(refl) == 1
    rings = ((fixed[0], ''), (refl[0], ' stroke-dasharray="3 2"'))
    out = []
    # (a) the cube, with the plane x0 = 0 at full strength and the rest faded
    frag, box = matrix_cube(T, 15, theme, opacity=lambda x0: 1 if x0 == 0 else 0.25,
                            segments=(0,), plane=True)
    w0, ht = box[2] - box[0], box[3] - box[1]
    out.append(place(frag, box, 0, 0))
    cap = max(ht, (ht - N * 24) / 2 + N * 24 + 40)
    out.append(label(w0 / 2, cap + 18, 'the support T, with the plane x₀ = 0', theme, 14))

    # (b) the plane x0 = 0 seen from above: x1 across, x2 away from the front face
    c = 24
    X = w0 + 44
    Y = (ht - N * c) / 2
    out.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" fill="%s" fill-opacity="0.16"/>'
               % (X, Y, N * c, N * c, Th['plane']))
    for i in range(N + 1):
        out.append('<path d="M%.1f %.1fv%dM%.1f %.1fh%d" stroke="%s" stroke-width="1"/>'
                   % (X + i * c, Y, N * c, X, Y + i * c, N * c,
                      Th['line'] if i in (0, N) else Th['faint']))

    def cellxy(x1, x2):
        return X + x1 * c + c / 2, Y + (N - 1 - x2) * c + c / 2
    for (a, b), (e, f) in (((0, 0), (N - 1, N - 1)), ((0, N - 1), (N - 1, 0))):
        (sx, sy), (ex, ey) = cellxy(a, b), cellxy(e, f)
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="1.2" opacity="0.45"/>' % (sx, sy, ex, ey, Th['text']))
    for x1 in range(N):
        cx, cy = cellxy(x1, p[x1])
        if p[x1] > 0:
            out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                       'stroke-width="1.6" stroke-linecap="round"/>'
                       % (cx, Y + N * c, cx, cy + c / 2 - 3, Th['text']))
        out.append('<circle cx="%.1f" cy="%.1f" r="2.2" fill="%s"/>' % (cx, Y + N * c, Th['text']))
        out.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" fill="%s"/>'
                   % (cx - c / 2 + 3, cy - c / 2 + 3, c - 6, c - 6, RED))
    for t, dash in rings:
        cx, cy = cellxy(t, p[t])
        out.append('<circle cx="%.1f" cy="%.1f" r="14" fill="none" stroke="%s" '
                   'stroke-width="1.6"%s/>' % (cx, cy, Th['text'], dash))
    for (x0_, y0_, x1_, y1_), name, (lx, ly) in (
            ((X, Y + N * c + 12, X + N * c, Y + N * c + 12), '1', (X + N * c / 2, Y + N * c + 32)),
            ((X - 12, Y + N * c, X - 12, Y), '2', (X - 26, Y + N * c / 2 + 5))):
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                   'stroke-width="1.2" marker-end="url(#arrow-%s)"/>'
                   % (x0_, y0_, x1_, y1_, Th['text'], theme))
        out.append(var(lx, ly, name, theme))
    out.append(label(X + N * c / 2, cap + 18, 'the plane x₀ = 0 from above', theme, 14))

    # (c) the Latin square, row 0 highlighted
    X2 = X + N * c + 44
    out.append(latin_panel(T, X2, Y, c, theme, rows=(0,), rings=rings))
    out.append(label(X2 + N * c + 8, Y + c / 2 + 5, 'row 0', theme, 12, anchor='start'))
    out.append(label(X2 + N * c / 2, cap + 18, 'its Latin square L(x₀, x₁) = x₂', theme, 14))
    W = X2 + N * c + 48
    t0, t1 = fixed[0], refl[0]
    out.append(label(W / 2, cap + 44,
                     'solid ring: the one fixed point, p(%d) = %d;  dashed ring: the one '
                     'reflected point, p(%d) = %d = 8 − %d' % (t0, p[t0], t1, p[t1], t1),
                     theme, 13))
    return svg(''.join(out), (-4, min(-4, Y - 26), W + 4, cap + 56), theme,
               'A support, its plane x0 = 0, and its Latin square')


ROTATING = {
    'support': fig_root,
    'two-supports': fig_pair,
}
FIGURES = {
    'orbit': fig_orbit,
    'companion-graph': fig_graph,
    'latin-depth': fig_depth,
    'latin-square': fig_lemma2,
}


def draw():
    d = json.load(open(DATA))
    for name, fn in FIGURES.items():
        for theme in THEMES:
            path = os.path.join(HERE, '%s-%s.svg' % (name, theme))
            with open(path, 'w') as f:
                f.write(fn(d, theme))
    print('wrote %d figures' % (2 * len(FIGURES)))


def gif(frames=90, fps=15, scale=2, cols=10):
    """One full turn about the x0 axis per figure and theme.  Frames are laid
    out on one page, rendered by headless Chrome in a single screenshot, cut
    apart with ImageMagick and assembled by ffmpeg with a fitted palette."""
    import base64
    import re
    import subprocess
    import tempfile
    chrome = os.environ.get(
        'CHROME', '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
    d = json.load(open(DATA))
    for name, fn in ROTATING.items():
        for theme in THEMES:
            svgs = [fn(d, theme, az=0.52 + 2 * math.pi * k / frames) for k in range(frames)]
            w, h = (int(v) for v in re.search(r'width="(\d+)" height="(\d+)"',
                                              svgs[0]).groups())
            rows = -(-frames // cols)
            with tempfile.TemporaryDirectory() as tmp:
                page = os.path.join(tmp, 'frames.html')
                with open(page, 'w') as f:
                    f.write('<html><body style="margin:0;background:%s;line-height:0">'
                            % THEMES[theme]['bg'])
                    for k, sv in enumerate(svgs):
                        f.write('<img style="display:inline-block" width="%d" height="%d" '
                                'src="data:image/svg+xml;base64,%s">'
                                % (w, h, base64.b64encode(sv.encode()).decode()))
                        if k % cols == cols - 1:
                            f.write('<br>')
                    f.write('</body></html>')
                shot = os.path.join(tmp, 'shot.png')
                subprocess.run([chrome, '--headless=new', '--disable-gpu',
                                '--hide-scrollbars', '--force-device-scale-factor=%d' % scale,
                                '--window-size=%d,%d' % (cols * w, rows * h),
                                '--screenshot=' + shot, 'file://' + page],
                               check=True, capture_output=True)
                subprocess.run(['magick', shot, '-crop', '%dx%d' % (w * scale, h * scale),
                                '+repage', os.path.join(tmp, 'tile_%03d.png')], check=True)
                # tiles run row by row, which is frame order; drop the padding tiles
                for k in range(frames, rows * cols):
                    t = os.path.join(tmp, 'tile_%03d.png' % k)
                    if os.path.exists(t):
                        os.remove(t)
                out = os.path.join(HERE, '%s-%s.gif' % (name, theme))
                subprocess.run(['ffmpeg', '-v', 'error', '-y', '-framerate', str(fps),
                                '-i', os.path.join(tmp, 'tile_%03d.png'), '-vf',
                                'split[a][b];[a]palettegen=stats_mode=full[p];'
                                '[b][p]paletteuse=dither=none', '-loop', '0', out],
                               check=True)
                print('wrote', out, os.path.getsize(out), 'bytes')


if __name__ == '__main__':
    {'extract': extract, 'draw': draw, 'gif': gif}[sys.argv[1] if len(sys.argv) > 1
                                                   else 'draw']()
