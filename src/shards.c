/* shards.c -- the shard universe: admissible row-0 permutations.
 *
 * Row 0 of a support of [n]^3 (its cells with i = 0) is forced to be
 *
 *     {(0, j, p(j)) : j in [n]}
 *
 * for a permutation p of [n]: the n lines (0,j,*) and the n lines (0,*,k) each
 * meet the support once, so the plane i = 0 holds exactly one cell per j and
 * exactly one per k.  The two face diagonals lying in that plane, {(0,t,t)} and
 * {(0,t,n-1-t)}, are each met exactly once too, so p has exactly one fixed
 * point (p(j) = j) and exactly one reflected point (p(j) = n-1-j).
 *
 * This program enumerates every permutation of [n] in lexicographic order and
 * keeps the admissible ones.  It is a brute force -- n! is 362 880 at n = 9 --
 * so the shard universe is complete by exhaustion rather than by argument, and
 * a shard *index* is simply the lexicographic rank among admissible rows.
 *
 * The count is OEIS A007016 (Simpson 1995): 0, 0, 8, 20, 96, 656, 5568, 48912
 * for n = 2..9.  `count` prints the sequence so it can be checked against OEIS.
 *
 * Build: cc -O2 -o shards shards.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SMAXN 12
#define MAXROWS 65536

static int n_side;
static int rows[MAXROWS][SMAXN];
static long long kept;
static int cur[SMAXN], used[SMAXN];
static int store_rows;

static int admissible(const int *p, int n)
{
    int f = 0, a = 0;
    for (int j = 0; j < n; j++) {
        if (p[j] == j) f++;
        if (p[j] == n - 1 - j) a++;
    }
    return f == 1 && a == 1;
}

/* every permutation of [n] in lexicographic order */
static void dfs(int j, int n)
{
    if (j == n) {
        if (!admissible(cur, n)) return;
        if (store_rows) {
            if (kept >= MAXROWS) { fprintf(stderr, "more than %d shards\n", MAXROWS); exit(1); }
            memcpy(rows[kept], cur, sizeof(int) * (size_t)n);
        }
        kept++;
        return;
    }
    for (int v = 0; v < n; v++) {
        if (used[v]) continue;
        used[v] = 1; cur[j] = v;
        dfs(j + 1, n);
        used[v] = 0;
    }
}

/* xorshift64*, so a shuffled shard order is reproducible from the seed alone */
static unsigned long long rng_state;
static unsigned long long rng_next(void)
{
    unsigned long long x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return x * 2685821657736338717ULL;
}

static void usage(int rc)
{
    fprintf(stderr,
"shards -- the admissible row-0 permutations of a support of [n]^3\n"
"\n"
"usage:\n"
"  shards <n> list [out.txt|-] [--shuffle SEED]\n"
"        Write one line per shard: `idx p0 p1 ... p<n-1>`, idx being the\n"
"        lexicographic rank.  --shuffle permutes the ORDER of the lines only\n"
"        (indices and rows are unchanged), which keeps a partial sweep an\n"
"        unbiased sample of the shard universe.\n"
"  shards <n> count\n"
"        Print the admissible count for each of 2..n, i.e. OEIS A007016.\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 3) usage(2);
    n_side = atoi(argv[1]);
    if (n_side < 2 || n_side > 10) { fprintf(stderr, "n out of range 2..10\n"); return 2; }

    if (!strcmp(argv[2], "count")) {
        printf("n\tadmissible_rows\n");
        for (int n = 2; n <= n_side; n++) {
            kept = 0; store_rows = 0;
            memset(used, 0, sizeof used);
            dfs(0, n);
            printf("%d\t%lld\n", n, kept);
        }
        return 0;
    }
    if (strcmp(argv[2], "list")) usage(2);

    const char *path = (argc > 3 && strcmp(argv[3], "-")) ? argv[3] : NULL;
    unsigned long long seed = 0;
    int shuffled = 0;
    for (int a = 3; a < argc; a++)
        if (!strcmp(argv[a], "--shuffle") && a + 1 < argc) {
            seed = strtoull(argv[a + 1], NULL, 10); shuffled = 1; a++;
        }

    /* Count first, and refuse an order that would not fit, BEFORE anything is
     * written into the fixed-size table. */
    kept = 0; store_rows = 0;
    memset(used, 0, sizeof used);
    dfs(0, n_side);
    if (kept > MAXROWS) {
        fprintf(stderr, "n=%d has %lld admissible rows; this build holds %d.\n"
                        "Orders up to 9 are supported; raise MAXROWS in "
                        "src/shards.c to go further.\n", n_side, kept, MAXROWS);
        return 2;
    }

    kept = 0; store_rows = 1;
    memset(used, 0, sizeof used);
    dfs(0, n_side);

    long long *order = malloc(sizeof(long long) * (size_t)kept);
    if (!order) { fprintf(stderr, "out of memory\n"); return 1; }
    for (long long i = 0; i < kept; i++) order[i] = i;
    if (shuffled) {
        rng_state = seed ? seed : 0x9E3779B97F4A7C15ULL;
        for (long long i = kept - 1; i > 0; i--) {
            long long j = (long long)(rng_next() % (unsigned long long)(i + 1));
            long long t = order[i]; order[i] = order[j]; order[j] = t;
        }
    }

    FILE *f = path ? fopen(path, "w") : stdout;
    if (!f) { perror(path); return 1; }
    for (long long i = 0; i < kept; i++) {
        long long idx = order[i];
        fprintf(f, "%lld", idx);
        for (int k = 0; k < n_side; k++) fprintf(f, " %d", rows[idx][k]);
        fputc('\n', f);
    }
    if (path && fclose(f)) { perror(path); return 1; }
    fprintf(stderr, "n=%d admissible rows: %lld%s\n", n_side, kept,
            shuffled ? " (order shuffled)" : "");
    free(order);
    return 0;
}
