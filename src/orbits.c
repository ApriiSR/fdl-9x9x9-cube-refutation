/* orbits.c -- the cell group of [n]^3, and orbit classification of supports.
 *
 * The main lines of [n]^3 are permuted among themselves by
 *
 *     g(x)_i = tau( eps_i( x_{pi(i)} ) ),      pi in S_3,
 *                                              eps in {id, rev}^3,
 *                                              tau in C_{S_n}(rev),
 *
 * where rev(t) = n-1-t.  A coordinate permutation pi maps axis lines to axis
 * lines and diagonals to diagonals; a per-axis reversal eps_i maps each family
 * to itself; and a symbol permutation tau applied to ALL THREE coordinates at
 * once maps a line {a + t b} to {tau(a) + t b} only when tau commutes with rev,
 * which is exactly the condition defining C(rev).  So the group maps supports
 * to supports.
 *
 * The parametrisation is 2-to-1: eps = (rev,rev,rev) with tau = id is the same
 * cell map as eps = id with tau = rev.  Hence
 *     |G| = 6 * 8 * |C(rev)| / 2,
 * and |C(rev)| = 2^{floor(n/2)} * floor(n/2)!  -- 384 at n = 8 and n = 9 --
 * giving |G| = 9216 at both orders.  This program builds G, checks that count
 * by deduplication, and picks a generating set greedily and deterministically
 * (verified by closure, not assumed).
 *
 * Every element of G fixes the centre cell of an odd cube: rev fixes (n-1)/2
 * and so does every tau in C(rev).  So G acts on the set of supports through
 * the centre, which is what `classify` is used for.
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

#define MAXG 20000

static int N, NC, RB;                 /* side, n^3, record bytes = n^2 */
static int *G;                        /* NG group elements, NC ints each */
static int NG;
static int gen[8];                    /* indices into G */
static int NGEN;

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ---- the centraliser of rev in S_n -------------------------------------- */
/* tau commutes with rev iff it permutes the pairs {t, n-1-t} as blocks,
 * possibly flipping each; the middle point of an odd n is fixed. */
static int taus[512][MAXN];
static int ntau;

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
            if (ntau > 512) { fprintf(stderr, "too many tau\n"); exit(1); }
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

static void build_group(void)
{
    build_taus();
    G = malloc(sizeof(int) * (size_t)MAXG * NC);
    if (!G) { fprintf(stderr, "out of memory\n"); exit(1); }
    NG = 0;
    /* dedup by hash of the cell map */
    size_t hs = 1 << 16;
    int *tab = malloc(sizeof(int) * hs);
    for (size_t i = 0; i < hs; i++) tab[i] = -1;

    int pi[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    int *g = malloc(sizeof(int) * NC);
    for (int a = 0; a < 6; a++)
    for (int e = 0; e < 8; e++)
    for (int ti = 0; ti < ntau; ti++) {
        const int *t = taus[ti];
        for (int x = 0; x < N; x++)
        for (int y = 0; y < N; y++)
        for (int z = 0; z < N; z++) {
            int src[3] = {x, y, z}, img[3];
            for (int i = 0; i < 3; i++) {
                int v = src[pi[a][i]];
                if (e & (1 << i)) v = N - 1 - v;
                img[i] = t[v];
            }
            g[(x * N + y) * N + z] = (img[0] * N + img[1]) * N + img[2];
        }
        uint64_t h = fnv(g, sizeof(int) * (size_t)NC);
        size_t s = h & (hs - 1);
        int dup = 0;
        while (tab[s] >= 0) {
            if (!memcmp(G + (size_t)tab[s] * NC, g, sizeof(int) * (size_t)NC)) { dup = 1; break; }
            s = (s + 1) & (hs - 1);
        }
        if (dup) continue;
        if (NG >= MAXG) { fprintf(stderr, "group too large\n"); exit(1); }
        memcpy(G + (size_t)NG * NC, g, sizeof(int) * (size_t)NC);
        tab[s] = NG;
        NG++;
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
    gtab = malloc(sizeof(int) * hs);
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
        seen = malloc((size_t)MAXG);
        stack = malloc(sizeof(int) * (size_t)MAXG);
        comp = malloc(sizeof(int) * (size_t)NC);
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

/* An explicit generating set, with a mathematical reason for each element,
 * VERIFIED by closure rather than assumed:
 *
 *   A  tau = the flip of the block {0, n-1}          (one reversal of symbols)
 *   B  tau = the transposition of the blocks {0,n-1} and {1,n-2}
 *   C  tau = the cyclic shift of all floor(n/2) blocks
 *   D  pi = the transposition (x,y,z) -> (y,x,z)
 *   E  pi = the 3-cycle      (x,y,z) -> (y,z,x)
 *   F  eps = reversal of the first coordinate
 *
 * A, B, C generate C(rev) (a flip, plus the symmetric group on the blocks);
 * D, E generate S_3 on the coordinates; F, conjugated by D and E, gives the
 * three per-axis reversals.  Together they must therefore be all of G, and the
 * closure computation checks that they are.
 */
static void build_element(const int *pi, int eps, const int *tau, int *g)
{
    for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
    for (int z = 0; z < N; z++) {
        int src[3] = {x, y, z}, img[3];
        for (int i = 0; i < 3; i++) {
            int v = src[pi[i]];
            if (eps & (1 << i)) v = N - 1 - v;
            img[i] = tau[v];
        }
        g[(x * N + y) * N + z] = (img[0] * N + img[1]) * N + img[2];
    }
}

static void add_generator(const int *pi, int eps, const int *tau, int *g)
{
    build_element(pi, eps, tau, g);
    int j = find_group(g);
    if (j < 0) { fprintf(stderr, "generator is not in the constructed group\n"); exit(1); }
    gen[NGEN++] = j;
}

static void pick_generators(void)
{
    int h = N / 2;
    int id3[3] = {0, 1, 2}, sw[3] = {1, 0, 2}, cy[3] = {1, 2, 0};
    int idt[MAXN], A[MAXN], B[MAXN], C[MAXN];
    int *g = malloc(sizeof(int) * (size_t)NC);
    for (int t = 0; t < N; t++) idt[t] = A[t] = B[t] = C[t] = t;
    A[0] = N - 1; A[N - 1] = 0;
    if (h >= 2) { B[0] = 1; B[1] = 0; B[N - 1] = N - 2; B[N - 2] = N - 1; }
    for (int a = 0; a < h; a++) { int b = (a + 1) % h; C[a] = b; C[N - 1 - a] = N - 1 - b; }

    NGEN = 0;
    add_generator(id3, 0, A, g);
    if (h >= 2) add_generator(id3, 0, B, g);
    if (h >= 3) add_generator(id3, 0, C, g);
    add_generator(sw, 0, idt, g);
    add_generator(cy, 0, idt, g);
    add_generator(id3, 1, idt, g);
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
    htab = malloc(sizeof(int) * hs);
    if (!htab) { fprintf(stderr, "out of memory for the record index\n"); exit(1); }
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
    for (int c = 0; c < RB; c++) {
        int cell = c * N + r[c];
        int img = g[cell];
        out[img / N] = (unsigned char)(img % N);
    }
}

static unsigned char *load(const char *path, long long *count)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    fseek(f, 0, SEEK_END);
    long long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz % RB) { fprintf(stderr, "%s: size is not a multiple of %d\n", path, RB); exit(1); }
    unsigned char *b = malloc((size_t)sz);
    if (!b) { fprintf(stderr, "out of memory reading %s\n", path); exit(1); }
    if (fread(b, 1, (size_t)sz, f) != (size_t)sz) { perror(path); exit(1); }
    fclose(f);
    *count = sz / RB;
    return b;
}

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
"  orbits <n> centre <recs.bin> <out.bin>\n"
"        Select the records containing the centre cell ((n-1)/2 three times).\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 3) usage(2);
    N = atoi(argv[1]);
    if (N < 2 || N > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    NC = N * N * N; RB = N * N;
    const char *mode = argv[2];

    if (!strcmp(mode, "centre")) {
        if (argc < 5) usage(2);
        int m = (N - 1) / 2;
        if (N % 2 == 0) { fprintf(stderr, "no centre cell for even n\n"); return 2; }
        recs = load(argv[3], &NR);
        FILE *o = fopen(argv[4], "wb");
        if (!o) { perror(argv[4]); return 1; }
        long long k = 0;
        for (long long i = 0; i < NR; i++)
            if (recs[i * RB + m * N + m] == m) { fwrite(recs + i * RB, 1, (size_t)RB, o); k++; }
        if (fclose(o)) { perror(argv[4]); return 1; }
        printf("{\"records\":%lld,\"centre\":[%d,%d,%d],\"selected\":%lld}\n", NR, m, m, m, k);
        return 0;
    }

    double t0 = now_s();
    build_group();
    index_group();
    pick_generators();
    if (!strcmp(mode, "group")) {
        printf("n=%d |C(rev)|=%d |G|=%d generators=%d (%.1f s)\n",
               N, ntau, NG, NGEN, now_s() - t0);
        int expect = 6 * 8 * ntau / 2;
        printf("expected |G| = 6*8*|C(rev)|/2 = %d -- %s\n", expect,
               expect == NG ? "ok" : "MISMATCH");
        if (N % 2) {
            int m = (N - 1) / 2, bad = 0;
            int c0 = (m * N + m) * N + m;
            for (int i = 0; i < NG; i++) if (G[(size_t)i * NC + c0] != c0) bad++;
            printf("elements not fixing the centre cell: %d\n", bad);
            if (bad) return 1;
        }
        return expect == NG ? 0 : 1;
    }

    if (!strcmp(mode, "closure")) {
        if (argc < 4) usage(2);
        recs = load(argv[3], &NR);
        index_records();
        unsigned char *img = malloc((size_t)RB);
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
        par = malloc(sizeof(int) * (size_t)NR);
        for (long long i = 0; i < NR; i++) par[i] = (int)i;
        unsigned char *img = malloc((size_t)RB);
        for (long long i = 0; i < NR; i++)
            for (int k = 0; k < NGEN; k++) {
                image(recs + i * RB, G + (size_t)gen[k] * NC, img);
                long long j = find_record(img);
                if (j < 0) { fprintf(stderr, "the set is not closed under the group\n"); return 1; }
                uni((int)i, (int)j);
            }
        /* orbit sizes, and the lexicographically least member of each orbit */
        int *size = calloc((size_t)NR, sizeof(int));
        int *least = malloc(sizeof(int) * (size_t)NR);
        for (long long i = 0; i < NR; i++) least[i] = -1;
        for (long long i = 0; i < NR; i++) {
            int r = find((int)i);
            size[r]++;
            if (least[r] < 0 || memcmp(recs + (size_t)i * RB, recs + (size_t)least[r] * RB, (size_t)RB) < 0)
                least[r] = (int)i;
        }
        long long norb = 0, tot = 0;
        int *reps = malloc(sizeof(int) * (size_t)NR);
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
