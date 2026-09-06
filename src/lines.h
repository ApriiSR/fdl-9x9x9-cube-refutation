/* lines.h -- the main lines of [n]^3, built from the definition.
 *
 * A *main line* of [n]^d is a set of n cells indexed by t = 0..n-1 in which
 * each coordinate is either a constant, or t, or n-1-t, with at least one
 * coordinate non-constant.  Reversing t maps a line to itself as a set, so we
 * enumerate each line once by requiring the first non-constant coordinate to
 * follow the pattern t.
 *
 * For d = 3 there are 3n^2 + 6n + 4 of them: 3n^2 axis lines, 6n face
 * diagonals, 4 space diagonals.  At n = 9 that is 301.
 *
 * Cells are indexed (i,j,k) -> (i*n + j)*n + k.
 */
#ifndef FDLH_LINES_H
#define FDLH_LINES_H

#define MAXN 12
#define MAXLINES (3 * MAXN * MAXN + 6 * MAXN + 4)
#define MAXCELLS (MAXN * MAXN * MAXN)
/* 3 axis + 6 face + 4 space lines can meet in one cell */
#define MAXDEG 16

/* Fills lines[][] with the main lines of [n]^3 and returns how many there are.
 * Deterministic: kinds are enumerated in the order c < t < r per coordinate,
 * constants in increasing order. */
static inline int fdlh_build_lines(int n, int lines[MAXLINES][MAXN])
{
    int nl = 0;
    int kind[3];
    /* 0 = constant, 1 = t, 2 = n-1-t */
    for (kind[0] = 0; kind[0] < 3; kind[0]++)
    for (kind[1] = 0; kind[1] < 3; kind[1]++)
    for (kind[2] = 0; kind[2] < 3; kind[2]++) {
        int nc = 0, first = -1;
        for (int a = 0; a < 3; a++)
            if (kind[a] != 0) { nc++; if (first < 0) first = a; }
        if (nc == 0) continue;               /* every coordinate constant */
        if (kind[first] != 1) continue;      /* canonical direction only */
        int nconst = 3 - nc;
        int total = 1;
        for (int a = 0; a < nconst; a++) total *= n;
        for (int code = 0; code < total; code++) {
            int c[3], rest = code;
            /* constants in increasing lexicographic order over the constant
             * positions, most significant first */
            for (int a = 2; a >= 0; a--)
                if (kind[a] == 0) { c[a] = rest % n; rest /= n; }
            for (int t = 0; t < n; t++) {
                int x[3];
                for (int a = 0; a < 3; a++)
                    x[a] = (kind[a] == 0) ? c[a] : (kind[a] == 1 ? t : n - 1 - t);
                lines[nl][t] = (x[0] * n + x[1]) * n + x[2];
            }
            nl++;
        }
    }
    return nl;
}

static inline int fdlh_expected_lines(int n) { return 3 * n * n + 6 * n + 4; }

#endif
