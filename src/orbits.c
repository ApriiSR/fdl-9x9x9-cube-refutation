/* orbits.c -- the cell group of [n]^3, and orbit classification of supports.
 *
 * The lines of [n]^3 are permuted among themselves by the maps
 *
 *     g(x)_i = sigma_i( x_{pi(i)} ),     pi in S_3,  sigma_0 in C_{S_n}(rev),
 *                                        sigma_1, sigma_2 in {sigma_0, sigma_0 . rev},
 *
 * where rev(t) = n-1-t and C(rev) is the set of permutations commuting with
 * it: those sending each pair {t, n-1-t} onto a pair.  Each output coordinate
 * is one input coordinate, relabeled by its own sigma_i; the three relabelings
 * must agree up to rev, which is what keeps the two varying coordinates of a
 * diagonal (always equal, or always partners) in step.  A constant coordinate
 * c goes to the constant sigma_i(c), and the patterns t and rev(t) go to
 * sigma_0(t) or rev(sigma_0(t)), so a line goes to a line.  (Reversing one axis
 * is sigma = (rev, id, id); relabeling every value by the same tau is
 * sigma_0 = sigma_1 = sigma_2 = tau.)  So the group maps supports to supports.
 *
 * pi and each sigma_i can be read back off the cell map, so every choice of
 * parameters gives a different map for n >= 2, and
 *     |G| = 6 * |C(rev)| * 2 * 2,
 * with |C(rev)| = 2^{floor(n/2)} * floor(n/2)!  -- 384 at n = 8 and n = 9 --
 * giving |G| = 9216 at both orders.  This program builds G, refuses a
 * duplicate map for n >= 2 rather than silently dropping it, and uses a fixed,
 * explicitly named generating set -- six elements chosen for a reason, not
 * searched for -- whose closure is computed and required to be the whole of G.
 *
 * symmetry.c constructs the same group again for its own use.  That is
 * duplicated code, deliberately left duplicated: the two programs are compared
 * by their outputs (|G| = 9216, and both report the same subgroup order), and
 * sharing the construction would remove that comparison.
 *
 * Every element of G fixes the center cell of an odd cube: rev fixes (n-1)/2
 * and so does every permutation in C(rev).  So G acts on the set of supports through
 * the center, which is what `classify` is used for.
 *
 * Records are n*n bytes: byte i*n+j is k for the cell (i,j,k).
 *
 * Build: cc -O2 -o orbits orbits.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "lines.h"
#include "util.h"

#define MAXG 20000
#define MAXTAU 512

static int N, NC, RB;                 /* side, n^3, record bytes = n^2 */
static int *G;                        /* NG group elements, NC ints each */
static int NG;
static int gen[8];                    /* indices into G */
static int NGEN;

#define now_s fdlh_now_s

/* ---- the centralizer of rev in S_n -------------------------------------- */
/* A permutation commutes with rev iff it permutes the pairs {t, n-1-t} as blocks,
 * possibly flipping each; the middle point of an odd n is fixed. */
static int taus[MAXTAU][MAXN];
static int ntau;

/* |C(rev)| = 2^{floor(n/2)} * floor(n/2)!, and |G| = 6 * |C(rev)| * 4.
 * Both arrays are fixed-size, so refuse an order that would overrun them
 * BEFORE anything is written, rather than after. */
static void check_capacity(int n)
{
    int h = n / 2;
    long long ctau = 1;
    for (int i = 2; i <= h; i++) ctau *= i;
    for (int i = 0; i < h; i++) ctau *= 2;
    long long cg = 6 * ctau * 4;
    if (ctau > MAXTAU || cg > MAXG) {
        fprintf(stderr,
                "order %d needs |C(rev)| = %lld and |G| = %lld; this build holds "
                "%d and %d.\nOrders up to 9 are supported; raise MAXTAU and MAXG "
                "in src/orbits.c to go further.\n", n, ctau, cg, MAXTAU, MAXG);
        exit(2);
    }
}

static void build_taus(void)
{
    int h = N / 2;
    int perm[MAXN];
    ntau = 0;

    /* all permutations of the h blocks, times all 2^h flips */
    long long nperm = 1;
    for (int i = 2; i <= h; i++) nperm *= i;
    for (long long pc = 0; pc < nperm; pc++) {
        /* Lehmer code -> permutation of the blocks */
        int avail[MAXN], na = h;
        for (int i = 0; i < h; i++) avail[i] = i;
        long long rest = pc;
        for (int i = 0; i < h; i++) {
            long long f = 1;
            for (int k = 2; k <= h - 1 - i; k++) f *= k;
            int q = (int)(rest / f); rest %= f;
            perm[i] = avail[q];
            for (int k = q; k < na - 1; k++) avail[k] = avail[k + 1];
            na--;
        }
        for (int flips = 0; flips < (1 << h); flips++) {
            if (ntau >= MAXTAU) { fprintf(stderr, "too many tau\n"); exit(1); }
            int *t = taus[ntau];
            for (int a = 0; a < h; a++) {
                int src0 = a, src1 = N - 1 - a;
                int b = perm[a];
                int d0 = b, d1 = N - 1 - b;
                if (flips & (1 << a)) { int s = d0; d0 = d1; d1 = s; }
                t[src0] = d0; t[src1] = d1;
            }
            if (N % 2) t[N / 2] = N / 2;
            ntau++;
        }
    }
}

/* ---- the cell group ------------------------------------------------------ */
static uint64_t fnv(const void *p, size_t nbytes)
{
    const unsigned char *b = p;
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < nbytes; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

/* g(x)_i = sigma_i(x_{pi(i)}), with sigma_0 = s0 and, for i = 1, 2, sigma_i =
 * s0 . rev when bit i-1 of b is set and s0 otherwise. */
static void build_element(const int *pi, const int *s0, int b, int *g)
{
    int sig[3][MAXN];
    for (int t = 0; t < N; t++) {
        sig[0][t] = s0[t];
        sig[1][t] = (b & 1) ? s0[N - 1 - t] : s0[t];
        sig[2][t] = (b & 2) ? s0[N - 1 - t] : s0[t];
    }
    for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
    for (int z = 0; z < N; z++) {
        int src[3] = {x, y, z}, img[3];
        for (int i = 0; i < 3; i++) img[i] = sig[i][src[pi[i]]];
        g[(x * N + y) * N + z] = (img[0] * N + img[1]) * N + img[2];
    }
}

static void build_group(void)
{
    build_taus();
    G = fdlh_alloc(sizeof(int) * (size_t)MAXG * NC, "the group");
    NG = 0;
    /* dedup by hash of the cell map */
    size_t hs = 1 << 16;
    int *tab = fdlh_alloc(sizeof(int) * hs, "group hash");
    for (size_t i = 0; i < hs; i++) tab[i] = -1;

    int pi[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    int *g = fdlh_alloc(sizeof(int) * (size_t)NC, "a group element");
    long long dups = 0;
    for (int a = 0; a < 6; a++)
    for (int ti = 0; ti < ntau; ti++)
    for (int b = 0; b < 4; b++) {
        build_element(pi[a], taus[ti], b, g);
        uint64_t h = fnv(g, sizeof(int) * (size_t)NC);
        size_t s = h & (hs - 1);
        int dup = 0;
        while (tab[s] >= 0) {
            if (!memcmp(G + (size_t)tab[s] * NC, g, sizeof(int) * (size_t)NC)) { dup = 1; break; }
            s = (s + 1) & (hs - 1);
        }
        if (dup) { dups++; continue; }
        if (NG >= MAXG) { fprintf(stderr, "group too large\n"); exit(1); }
        memcpy(G + (size_t)NG * NC, g, sizeof(int) * (size_t)NC);
        tab[s] = NG;
        NG++;
    }
    /* below n = 2 every map is the identity, so there duplicates are expected */
    if (N >= 2 && dups) {
        fprintf(stderr, "%lld parameter choices repeat a map already built\n", dups);
        exit(1);
    }
    free(g); free(tab);
}

/* ---- a hash index of the group elements, so composition can be looked up -- */
static int *gtab;
static size_t gmask;
static int g_identity;

static void index_group(void)
{
    size_t hs = 1;
    while (hs < (size_t)NG * 4) hs <<= 1;
    gmask = hs - 1;
    gtab = fdlh_alloc(sizeof(int) * hs, "group index");
    for (size_t i = 0; i < hs; i++) gtab[i] = -1;
    for (int i = 0; i < NG; i++) {
        size_t s = fnv(G + (size_t)i * NC, sizeof(int) * (size_t)NC) & gmask;
        while (gtab[s] >= 0) s = (s + 1) & gmask;
        gtab[s] = i;
    }
    g_identity = -1;
    for (int i = 0; i < NG && g_identity < 0; i++) {
        int ok = 1;
        for (int c = 0; c < NC; c++) if (G[(size_t)i * NC + c] != c) { ok = 0; break; }
        if (ok) g_identity = i;
    }
    if (g_identity < 0) { fprintf(stderr, "no identity in the group\n"); exit(1); }
}

static int find_group(const int *g)
{
    size_t s = fnv(g, sizeof(int) * (size_t)NC) & gmask;
    while (gtab[s] >= 0) {
        if (!memcmp(G + (size_t)gtab[s] * NC, g, sizeof(int) * (size_t)NC)) return gtab[s];
        s = (s + 1) & gmask;
    }
    return -1;
}

/* size of the subgroup generated by gens[0..ngens) */
static int closure_size(const int *gens, int ngens)
{
    static unsigned char *seen;
    static int *stack, *comp;
    if (!seen) {
        seen = fdlh_alloc((size_t)MAXG, "closure marks");
        stack = fdlh_alloc(sizeof(int) * (size_t)MAXG, "closure stack");
        comp = fdlh_alloc(sizeof(int) * (size_t)NC, "a composition");
    }
    memset(seen, 0, (size_t)NG);
    int sp = 0, cnt = 1;
    seen[g_identity] = 1; stack[sp++] = g_identity;
    while (sp) {
        int a = stack[--sp];
        for (int k = 0; k < ngens; k++) {
            const int *ga = G + (size_t)a * NC, *gb = G + (size_t)gens[k] * NC;
            for (int c = 0; c < NC; c++) comp[c] = gb[ga[c]];
            int j = find_group(comp);
            if (j < 0) { fprintf(stderr, "group is not closed under composition\n"); exit(1); }
            if (!seen[j]) { seen[j] = 1; stack[sp++] = j; cnt++; }
        }
    }
    return cnt;
}

/* An explicit generating set -- fixed, not searched for -- with a mathematical
 * reason for each element, and VERIFIED by closure rather than assumed:
 *
 *   A  sigma_i = tau for all i, tau the flip of the block {0, n-1}
 *   B  sigma_i = tau for all i, tau the transposition of the blocks {0,n-1}
 *      and {1,n-2}
 *   C  sigma_i = tau for all i, tau the cyclic shift of all floor(n/2) blocks
 *   D  pi = the transposition (x,y,z) -> (y,x,z)
 *   E  pi = the 3-cycle      (x,y,z) -> (y,z,x)
 *   F  sigma = (rev, id, id), reversal of the first coordinate
 *
 * A, B, C give every sigma_0 = sigma_1 = sigma_2 in C(rev) (a flip, plus the
 * symmetric group on the blocks); D, E generate S_3 on the coordinates; F,
 * conjugated by D and E, reverses any one axis, which with the equal-sigma
 * elements gives every sigma_i in {sigma_0, sigma_0 . rev}.  Together they must therefore be all of G, and the
 * closure computation checks that they are.
 */
static void add_generator(const int *pi, const int *s0, int b, int *g)
{
    build_element(pi, s0, b, g);
    int j = find_group(g);
    if (j < 0) { fprintf(stderr, "generator is not in the constructed group\n"); exit(1); }
    gen[NGEN++] = j;
}

static void pick_generators(void)
{
    int h = N / 2;
    int id3[3] = {0, 1, 2}, sw[3] = {1, 0, 2}, cy[3] = {1, 2, 0};
    int idt[MAXN], A[MAXN], B[MAXN], C[MAXN], R[MAXN];
    int *g = fdlh_alloc(sizeof(int) * (size_t)NC, "a generator");
    for (int t = 0; t < N; t++) { idt[t] = A[t] = B[t] = C[t] = t; R[t] = N - 1 - t; }
    A[0] = N - 1; A[N - 1] = 0;
    if (h >= 2) { B[0] = 1; B[1] = 0; B[N - 1] = N - 2; B[N - 2] = N - 1; }
    for (int a = 0; a < h; a++) { int b = (a + 1) % h; C[a] = b; C[N - 1 - a] = N - 1 - b; }

    NGEN = 0;
    add_generator(id3, A, 0, g);
    if (h >= 2) add_generator(id3, B, 0, g);
    if (h >= 3) add_generator(id3, C, 0, g);
    add_generator(sw, idt, 0, g);
    add_generator(cy, idt, 0, g);
    add_generator(id3, R, 3, g);           /* sigma = (rev, rev.rev, rev.rev) */
    free(g);

    int c = closure_size(gen, NGEN);
    if (c != NG) {
        fprintf(stderr, "the %d chosen generators generate only %d of %d elements\n",
                NGEN, c, NG);
        exit(1);
    }
}

/* ---- records and a hash index ------------------------------------------- */
static unsigned char *recs;
static long long NR;
static int *htab;
static size_t hmask;

static void index_records(void)
{
    size_t hs = 1;
    while (hs < (size_t)NR * 2) hs <<= 1;
    hmask = hs - 1;
    htab = fdlh_alloc(sizeof(int) * hs, "the record index");
    for (size_t i = 0; i < hs; i++) htab[i] = -1;
    for (long long i = 0; i < NR; i++) {
        const unsigned char *r = recs + i * RB;
        size_t s = fnv(r, (size_t)RB) & hmask;
        while (htab[s] >= 0) {
            if (!memcmp(recs + (size_t)htab[s] * RB, r, (size_t)RB)) {
                fprintf(stderr, "duplicate record at %lld\n", i); exit(1);
            }
            s = (s + 1) & hmask;
        }
        htab[s] = (int)i;
    }
}

static long long find_record(const unsigned char *r)
{
    size_t s = fnv(r, (size_t)RB) & hmask;
    while (htab[s] >= 0) {
        if (!memcmp(recs + (size_t)htab[s] * RB, r, (size_t)RB)) return htab[s];
        s = (s + 1) & hmask;
    }
    return -1;
}

static void image(const unsigned char *r, const int *g, unsigned char *out)
{
    memset(out, 0xff, (size_t)RB);
    for (int c = 0; c < RB; c++) {
        int cell = c * N + r[c];
        int img = g[cell];
        out[img / N] = (unsigned char)(img % N);
    }
}

/* Every byte is checked to be a coordinate in [n] on the way in: image() and
 * the mask builders use them as indices. */
static unsigned char *load(const char *path, long long *count)
{
    return fdlh_load_records(path, RB, N, count);
}

/* image() writes every one of the RB positions exactly once when g is a
 * bijection of the cells, but the destination is cleared first so that a
 * partially written buffer can never be mistaken for a record. */

/* ---- union-find ---------------------------------------------------------- */
static int *par;
static int find(int a) { while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; } return a; }
static void uni(int a, int b) { int ra = find(a), rb = find(b); if (ra != rb) par[ra] = rb; }

static void usage(int rc)
{
    fprintf(stderr,
"orbits -- the order-9216 cell group of [n]^3 and its action on supports\n"
"\n"
"usage:\n"
"  orbits <n> group\n"
"        Report |C(rev)|, |G|, and the deterministically chosen generators\n"
"        (checked by closure to generate the whole group).\n"
"  orbits <n> closure <recs.bin>\n"
"        Check that the given set of supports is closed under the generators:\n"
"        every image of every record must be in the set.  Exits nonzero if not.\n"
"  orbits <n> classify <recs.bin> <out.jsonl> <reps.bin>\n"
"        Decompose the records into G-orbits.  Writes one JSON line per orbit\n"
"        (rep index, orbit size) sorted by representative, and the\n"
"        representatives -- the lexicographically least record of each orbit,\n"
"        so the choice does not depend on the input order -- to reps.bin,\n"
"        sorted lexicographically.\n"
"  orbits <n> center <recs.bin> <out.bin>\n"
"        Select the records containing the center cell ((n-1)/2 three times).\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 3) usage(2);
    N = atoi(argv[1]);
    if (N < 0 || N > MAXN) { fprintf(stderr, "n out of range 0..%d\n", MAXN); return 2; }
    check_capacity(N);
    NC = N * N * N; RB = N * N;
    const char *mode = argv[2];

    if (!strcmp(mode, "center")) {
        if (argc < 5) usage(2);
        int m = (N - 1) / 2;
        if (N % 2 == 0) { fprintf(stderr, "no center cell for even n\n"); return 2; }
        recs = load(argv[3], &NR);
        FILE *o = fopen(argv[4], "wb");
        if (!o) { perror(argv[4]); return 1; }
        long long k = 0;
        for (long long i = 0; i < NR; i++)
            if (recs[i * RB + m * N + m] == m) { fwrite(recs + i * RB, 1, (size_t)RB, o); k++; }
        if (fclose(o)) { perror(argv[4]); return 1; }
        printf("{\"records\":%lld,\"center\":[%d,%d,%d],\"selected\":%lld}\n", NR, m, m, m, k);
        return 0;
    }

    double t0 = now_s();
    build_group();
    index_group();
    pick_generators();
    if (!strcmp(mode, "group")) {
        printf("n=%d |C(rev)|=%d |G|=%d generators=%d (%.1f s)\n",
               N, ntau, NG, NGEN, now_s() - t0);
        int expect = N >= 2 ? 6 * ntau * 4 : 1;   /* below n = 2 every map is the identity */
        printf("expected |G| = %s = %d -- %s\n", N >= 2 ? "6*|C(rev)|*4" : "1 (n < 2)", expect,
               expect == NG ? "ok" : "MISMATCH");
        if (N % 2) {
            int m = (N - 1) / 2, bad = 0;
            int c0 = (m * N + m) * N + m;
            for (int i = 0; i < NG; i++) if (G[(size_t)i * NC + c0] != c0) bad++;
            printf("elements not fixing the center cell: %d\n", bad);
            if (bad) return 1;
        }
        return expect == NG ? 0 : 1;
    }

    if (!strcmp(mode, "closure")) {
        if (argc < 4) usage(2);
        recs = load(argv[3], &NR);
        index_records();
        unsigned char *img = fdlh_alloc((size_t)RB, "an image");
        long long checked = 0, missing = 0;
        for (long long i = 0; i < NR; i++)
            for (int k = 0; k < NGEN; k++) {
                image(recs + i * RB, G + (size_t)gen[k] * NC, img);
                checked++;
                if (find_record(img) < 0) { missing++; if (missing < 5) fprintf(stderr, "missing image of record %lld under generator %d\n", i, k); }
            }
        printf("{\"records\":%lld,\"generators\":%d,\"images_checked\":%lld,"
               "\"missing\":%lld,\"wall\":%.1f}\n", NR, NGEN, checked, missing, now_s() - t0);
        return missing ? 1 : 0;
    }

    if (!strcmp(mode, "classify")) {
        if (argc < 6) usage(2);
        recs = load(argv[3], &NR);
        index_records();
        par = fdlh_alloc(sizeof(int) * (size_t)NR, "union-find");
        for (long long i = 0; i < NR; i++) par[i] = (int)i;
        unsigned char *img = fdlh_alloc((size_t)RB, "an image");
        for (long long i = 0; i < NR; i++)
            for (int k = 0; k < NGEN; k++) {
                image(recs + i * RB, G + (size_t)gen[k] * NC, img);
                long long j = find_record(img);
                if (j < 0) { fprintf(stderr, "the set is not closed under the group\n"); return 1; }
                uni((int)i, (int)j);
            }
        /* orbit sizes, and the lexicographically least member of each orbit */
        int *size = fdlh_calloc((size_t)NR, sizeof(int), "orbit sizes");
        int *least = fdlh_alloc(sizeof(int) * (size_t)NR, "orbit minima");
        for (long long i = 0; i < NR; i++) least[i] = -1;
        for (long long i = 0; i < NR; i++) {
            int r = find((int)i);
            size[r]++;
            if (least[r] < 0 || memcmp(recs + (size_t)i * RB, recs + (size_t)least[r] * RB, (size_t)RB) < 0)
                least[r] = (int)i;
        }
        long long norb = 0, tot = 0;
        int *reps = fdlh_alloc(sizeof(int) * (size_t)NR, "representatives");
        for (long long i = 0; i < NR; i++)
            if (find((int)i) == (int)i) { reps[norb++] = least[i]; tot += size[i]; }
        /* sort representatives lexicographically so the output is canonical */
        for (long long a = 1; a < norb; a++) {         /* insertion sort: few orbits */
            int v = reps[a]; long long b = a - 1;
            while (b >= 0 && memcmp(recs + (size_t)reps[b] * RB, recs + (size_t)v * RB, (size_t)RB) > 0) {
                reps[b + 1] = reps[b]; b--;
            }
            reps[b + 1] = v;
        }
        FILE *jo = fopen(argv[4], "w");
        FILE *bo = fopen(argv[5], "wb");
        if (!jo || !bo) { perror("output"); return 1; }
        for (long long a = 0; a < norb; a++) {
            int r = find(reps[a]);
            fprintf(jo, "{\"j\":%lld,\"orbit_size\":%d}\n", a, size[r]);
            fwrite(recs + (size_t)reps[a] * RB, 1, (size_t)RB, bo);
        }
        if (fclose(jo) || fclose(bo)) { perror("output"); return 1; }
        printf("{\"records\":%lld,\"group_order\":%d,\"generators\":%d,\"orbits\":%lld,"
               "\"sum_orbit_sizes\":%lld,\"wall\":%.1f}\n",
               NR, NG, NGEN, norb, tot, now_s() - t0);
        if (tot != NR) { fprintf(stderr, "orbit sizes do not sum to the record count\n"); return 1; }
        for (long long a = 0; a < norb; a++)
            if (NG % size[find(reps[a])]) { fprintf(stderr, "orbit size does not divide |G|\n"); return 1; }
        return 0;
    }

    usage(2);
    return 2;
}
