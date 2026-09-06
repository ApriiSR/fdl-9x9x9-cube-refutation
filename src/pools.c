/* pools.c -- companion pools.
 *
 * The *companion pool* of a support T is the set of all supports disjoint from
 * T.  Any packing of pairwise disjoint supports that contains T uses only
 * members of T's pool, so a complete catalogue turns "are there n-1 further
 * supports, disjoint from T and from each other?" into a question about a
 * finite, explicitly listed set -- which is the whole reason for paying to
 * enumerate the catalogue once.
 *
 * Completeness of a pool is a property of the *filter*, not of a second search:
 * given a complete catalogue, scanning it and keeping the disjoint records
 * cannot miss anything.  (check.py re-derives every pool independently in
 * numpy, so the filter itself is checked rather than believed.)
 *
 * Records are n*n bytes: byte i*n+j is k for the cell (i,j,k).  Both the
 * catalogue and the query file are in that format.  Pool j is written to
 * <outdir>/pool_<j>.bin in catalogue order.
 *
 * Build: cc -O2 -o pools pools.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/stat.h>

#define PMAXN 12
#define MAXW (PMAXN * PMAXN * PMAXN / 64 + 1)

static int N, NC, RB, W;              /* side, n^3, n^2 bytes, mask words */

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static unsigned char *load(const char *path, long long *count, int rb)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    fseek(f, 0, SEEK_END);
    long long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz % rb) { fprintf(stderr, "%s: size is not a multiple of %d\n", path, rb); exit(1); }
    unsigned char *b = malloc((size_t)sz ? (size_t)sz : 1);
    if (!b) { fprintf(stderr, "out of memory reading %s\n", path); exit(1); }
    if (sz && fread(b, 1, (size_t)sz, f) != (size_t)sz) { perror(path); exit(1); }
    fclose(f);
    *count = sz / rb;
    return b;
}

static void mask_of(const unsigned char *r, uint64_t *m)
{
    memset(m, 0, sizeof(uint64_t) * (size_t)W);
    for (int c = 0; c < RB; c++) {
        int cell = c * N + r[c];
        m[cell >> 6] |= (uint64_t)1 << (cell & 63);
    }
}

static int popcount_mask(const uint64_t *m)
{
    int s = 0;
    for (int w = 0; w < W; w++) s += __builtin_popcountll(m[w]);
    return s;
}

static void usage(int rc)
{
    fprintf(stderr,
"pools -- build the companion pool of each query support\n"
"\n"
"usage:\n"
"  pools <n> build <catalogue.bin> <queries.bin> <outdir> <sizes.jsonl>\n"
"        For each record of queries.bin, write the catalogue records disjoint\n"
"        from it to <outdir>/pool_<j>.bin, in catalogue order, and append one\n"
"        JSON line giving the pool size.  One catalogue load serves every\n"
"        query.  Every catalogue record is checked to have n^2 distinct cells\n"
"        on the way in.\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 7 || strcmp(argv[2], "build")) usage(2);
    N = atoi(argv[1]);
    if (N < 2 || N > PMAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    NC = N * N * N; RB = N * N; W = (NC + 63) / 64;

    double t0 = now_s();
    long long NT, NQ;
    unsigned char *cat = load(argv[3], &NT, RB);
    unsigned char *qry = load(argv[4], &NQ, RB);
    const char *outdir = argv[5];
    mkdir(outdir, 0777);

    uint64_t *cm = malloc(sizeof(uint64_t) * (size_t)NT * W);
    if (!cm) { fprintf(stderr, "out of memory for %lld catalogue masks\n", NT); return 1; }
    for (long long i = 0; i < NT; i++) {
        mask_of(cat + i * RB, cm + i * W);
        if (popcount_mask(cm + i * W) != RB) {
            fprintf(stderr, "catalogue record %lld does not have %d distinct cells\n", i, RB);
            return 1;
        }
    }
    fprintf(stderr, "catalogue %lld records, queries %lld (loaded in %.1f s)\n",
            NT, NQ, now_s() - t0);

    FILE *sz = fopen(argv[6], "w");
    if (!sz) { perror(argv[6]); return 1; }

    uint64_t qm[MAXW];
    size_t cap = 1 << 14;
    unsigned char *buf = malloc(cap * (size_t)RB);
    long long minp = -1, maxp = -1, totp = 0;
    for (long long j = 0; j < NQ; j++) {
        mask_of(qry + j * RB, qm);
        long long k = 0;
        for (long long i = 0; i < NT; i++) {
            const uint64_t *m = cm + i * W;
            int ok = 1;
            for (int w = 0; w < W; w++) if (m[w] & qm[w]) { ok = 0; break; }
            if (!ok) continue;
            if ((size_t)k == cap) { cap *= 2; buf = realloc(buf, cap * (size_t)RB); }
            memcpy(buf + k * RB, cat + i * RB, (size_t)RB);
            k++;
        }
        char path[1024];
        snprintf(path, sizeof path, "%s/pool_%lld.bin", outdir, j);
        FILE *f = fopen(path, "wb");
        if (!f) { perror(path); return 1; }
        if (k && fwrite(buf, RB, (size_t)k, f) != (size_t)k) { perror(path); return 1; }
        if (fclose(f)) { perror(path); return 1; }
        fprintf(sz, "{\"j\":%lld,\"pool\":%lld}\n", j, k);
        totp += k;
        if (minp < 0 || k < minp) minp = k;
        if (k > maxp) maxp = k;
        if ((j + 1) % 250 == 0)
            fprintf(stderr, "[%.0f s] %lld/%lld pools\n", now_s() - t0, j + 1, NQ);
    }
    if (fclose(sz)) { perror(argv[6]); return 1; }
    printf("{\"catalogue\":%lld,\"queries\":%lld,\"pool_min\":%lld,\"pool_max\":%lld,"
           "\"pool_mean\":%.1f,\"wall\":%.1f}\n",
           NT, NQ, minp, maxp, NQ ? (double)totp / (double)NQ : 0.0, now_s() - t0);
    return 0;
}
