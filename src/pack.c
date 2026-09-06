/* pack.c -- packing supports: exhaustive exact cover, and the packing ceiling.
 *
 * A fully diagonalised Latin cube of order n is a partition of [n]^3 into n
 * pairwise disjoint supports, so it is an exact cover of the n^3 cells by
 * elements of the support catalogue.  Two questions are asked here.
 *
 * EXHAUSTION.  Depth-first exact cover over a pool of supports.  At each node
 * the *cell* branched on is the uncovered cell lying in the fewest surviving
 * supports (lowest cell index winning ties); since that choice is a function of
 * the set already chosen and not of the order it was chosen in, each unordered
 * partition is reached along exactly one path and is counted once.  A branch
 * dies as soon as some uncovered cell has no surviving candidate.
 *
 *   census   all covers of [n]^3 from a complete catalogue (198 624 at n = 8),
 *            with the node count at each depth and, for every support, the
 *            number of covers containing it;
 *   roots    for each query support, the covers containing it, searched over
 *            that support's companion pool only.  A wall-clock cap is
 *            available and is recorded in the output when it fires, so a
 *            truncated search can never be read as an exhaustion.
 *
 * THE PACKING CEILING.  A second instrument on the same data: build the graph
 * on a pool with an edge between disjoint members.  A clique of size k there is
 * a packing of k+1 pairwise disjoint supports containing the query, and a cube
 * of order n needs n.  Edges, triangles and the maximum clique are all
 * reported.
 *
 * Two of those three are counts rather than searches, and that is the point.  A
 * k-clique contains C(k,3) triangles, so the triangle count alone bounds the
 * clique number: at n = 9 a witness would need an 8-clique and therefore 56
 * triangles, and a pool graph with fewer cannot contain one.  The maximum
 * clique is computed as well, by the recursive branch and bound bk() below --
 * a different algorithm from the exact cover, though the two share this file's
 * record loader, cell masks and disjoint().
 *
 * Records are n*n bytes: byte i*n+j is k for the cell (i,j,k).
 *
 * Build: cc -O2 -o pack pack.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "util.h"

#define KMAXN 12
#define MAXW (KMAXN * KMAXN * KMAXN / 64 + 1)

static int N, NC, RB, W;

#define now_s fdlh_now_s

static void mask_of(const unsigned char *r, uint64_t *m)
{
    memset(m, 0, sizeof(uint64_t) * (size_t)W);
    for (int c = 0; c < RB; c++) {
        int cell = c * N + r[c];
        m[cell >> 6] |= (uint64_t)1 << (cell & 63);
    }
}

static int disjoint(const uint64_t *a, const uint64_t *b)
{
    for (int w = 0; w < W; w++) if (a[w] & b[w]) return 0;
    return 1;
}

/* ---- exact cover --------------------------------------------------------- */
static uint64_t *M;              /* pool masks, W words each */
static long long *cover_count;   /* per-pool-member cover tally (census mode) */
static long long nodes_by_depth[KMAXN + 2];
static long long nsol;
static int max_depth_seen;
static double cap_seconds, t_start;
static int timed_out;
static long long node_check;
static int chosen[KMAXN + 1];

/* The cap is a deadline on the whole query, not on the search alone: the graph
 * construction, the triangle count and the clique search all consult it, so a
 * query that hits it is reported as BUDGET whichever stage was running. */
static int over_budget(long long every)
{
    if (timed_out) return 1;
    if (cap_seconds <= 0.0) return 0;
    if (++node_check < every) return 0;
    node_check = 0;
    if (now_s() - t_start > cap_seconds) timed_out = 1;
    return timed_out;
}

/* candidate lists per depth, so recursion needs no allocation */
static int *cand_store;
static int cand_stride;

static void search(uint64_t *cov, int depth, const int *surv, int nsurv)
{
    if (timed_out) return;
    nodes_by_depth[depth]++;
    if (depth > max_depth_seen) max_depth_seen = depth;
    if (depth == N) { nsol++; if (cover_count) for (int a = 0; a < N; a++) cover_count[chosen[a]]++; return; }
    if (over_budget(64)) return;
    /* the uncovered cell lying in the fewest survivors */
    int bestcell = -1, bestn = 1 << 30;
    for (int c = 0; c < NC; c++) {
        if (cov[c >> 6] & ((uint64_t)1 << (c & 63))) continue;
        int cnt = 0;
        for (int a = 0; a < nsurv; a++)
            if (M[(size_t)surv[a] * W + (c >> 6)] & ((uint64_t)1 << (c & 63))) cnt++;
        if (cnt < bestn) { bestn = cnt; bestcell = c; if (!cnt) break; }
    }
    if (bestcell < 0 || bestn == 0) return;

    int *cand = cand_store + (size_t)depth * cand_stride;
    int nc = 0;
    for (int a = 0; a < nsurv; a++)
        if (M[(size_t)surv[a] * W + (bestcell >> 6)] & ((uint64_t)1 << (bestcell & 63)))
            cand[nc++] = surv[a];

    int *next = cand_store + (size_t)(depth + N + 1) * cand_stride;
    for (int a = 0; a < nc; a++) {
        int s = cand[a];
        const uint64_t *ms = M + (size_t)s * W;
        uint64_t ncov[MAXW];
        for (int w = 0; w < W; w++) ncov[w] = cov[w] | ms[w];
        int nn = 0;
        for (int b = 0; b < nsurv; b++) {
            int u = surv[b];
            if (u == s) continue;
            if (disjoint(M + (size_t)u * W, ms)) next[nn++] = u;
        }
        chosen[depth] = s;
        /* deeper levels use slots depth+1 and depth+N+2, so `cand` (slot depth)
         * and `next` (slot depth+N+1) both survive the recursive call */
        search(ncov, depth + 1, next, nn);
        if (timed_out) break;
    }
}

/* ---- maximum clique in the disjointness graph ---------------------------- */
static uint64_t *adj;            /* AW words per vertex */
static int AW;

static int popcnt(const uint64_t *a)
{
    int s = 0;
    for (int w = 0; w < AW; w++) s += __builtin_popcountll(a[w]);
    return s;
}

static int best_clique;
static int cur_set[64], best_set[64];
static void bk(uint64_t *cand, int size)
{
    if (over_budget(64)) return;
    int nc = popcnt(cand);
    if (size + nc <= best_clique) return;
    if (!nc) {
        if (size > best_clique) { best_clique = size; memcpy(best_set, cur_set, sizeof(int) * (size_t)size); }
        return;
    }
    for (int w = 0; w < AW; w++) {
        uint64_t bits = cand[w];
        while (bits) {
            int b = __builtin_ctzll(bits);
            bits &= bits - 1;
            int v = w * 64 + b;
            if (size + popcnt(cand) <= best_clique) return;
            uint64_t *nx = fdlh_alloc(sizeof(uint64_t) * (size_t)AW, "clique candidates");
            for (int u = 0; u < AW; u++) nx[u] = cand[u] & adj[(size_t)v * AW + u];
            if (size < 64) cur_set[size] = v;
            bk(nx, size + 1);
            free(nx);
            cand[w] &= ~((uint64_t)1 << b);      /* v exhausted; drop it */
            if (timed_out) return;
        }
    }
    if (size > best_clique) { best_clique = size; memcpy(best_set, cur_set, sizeof(int) * (size_t)size); }
}

static void usage(int rc)
{
    fprintf(stderr,
"pack -- exhaustive exact cover over a support pool, and the packing ceiling\n"
"\n"
"usage:\n"
"  pack <n> census <catalogue.bin> <out.json> [--cap SECONDS]\n"
"        Enumerate every partition of [n]^3 into n supports drawn from the\n"
"        catalogue.  Prints the total, the node count at each depth, and\n"
"        writes per-support cover tallies to out.json.\n"
"  pack <n> roots <queries.bin> <pooldir> <out.jsonl> [options]\n"
"        For each query j, search its pool <pooldir>/pool_<j>.bin for a cover\n"
"        containing it, and independently compute its pool graph's edges,\n"
"        triangles and maximum clique.  One JSON line per query.\n"
"        options:\n"
"          --slice K W    take only the queries with j = K (mod W)\n"
"          --cap SECONDS  wall-clock deadline for the WHOLE query -- the graph,\n"
"                         the triangle count and the clique search as well as\n"
"                         the exact cover.  A query that hits it is reported as\n"
"                         BUDGET, never as EXHAUSTED.  Set FDLH_CLOCK_STEP to a\n"
"                         positive number of seconds to drive the deadline from\n"
"                         a deterministic virtual clock instead.\n"
"          --no-clique    skip the clique computation (exhaustion only)\n"
"          --no-search    skip the exact cover (clique only)\n"
"  pack <n> witness <queries.bin> <pooldir> <j> <out.json>\n"
"        Write query j together with a maximum clique of its pool graph, as an\n"
"        explicit packing of pairwise disjoint supports (check.py witness\n"
"        re-checks it from the definition).\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 5) usage(2);
    N = atoi(argv[1]);
    if (N < 2 || N > KMAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    NC = N * N * N; RB = N * N; W = (NC + 63) / 64;

    if (!strcmp(argv[2], "census")) {
        for (int a = 5; a < argc; a++)
            if (!strcmp(argv[a], "--cap") && a + 1 < argc) { cap_seconds = atof(argv[a + 1]); a++; }
            else usage(2);
        long long NT;
        unsigned char *cat = fdlh_load_records(argv[3], RB, N, &NT);
        M = fdlh_alloc(sizeof(uint64_t) * (size_t)NT * W, "catalogue masks");
        for (long long i = 0; i < NT; i++) mask_of(cat + i * RB, M + i * W);
        cover_count = fdlh_calloc((size_t)NT, sizeof(long long), "cover tallies");
        int *surv = fdlh_alloc(sizeof(int) * (size_t)NT, "survivors");
        for (long long i = 0; i < NT; i++) surv[i] = (int)i;
        cand_stride = (int)NT;
        cand_store = fdlh_alloc(sizeof(int) * (size_t)cand_stride * (2 * N + 2), "candidate lists");
        uint64_t cov[MAXW];
        memset(cov, 0, sizeof cov);
        node_check = 63;
        t_start = now_s();
        search(cov, 0, surv, (int)NT);
        double wall = now_s() - t_start;
        printf("{\"n\":%d,\"catalogue\":%lld,\"covers\":%lld,\"status\":\"%s\","
               "\"wall\":%.3f,\"nodes_by_depth\":[", N, NT, nsol,
               timed_out ? "BUDGET" : "EXHAUSTED", wall);
        for (int d = 0; d <= N; d++) printf("%lld%s", nodes_by_depth[d], d < N ? "," : "");
        printf("]}\n");
        FILE *f = fopen(argv[4], "w");
        if (!f) { perror(argv[4]); return 1; }
        fprintf(f, "{\"n\":%d,\"catalogue\":%lld,\"covers\":%lld,\"status\":\"%s\",\n",
                N, NT, nsol, timed_out ? "BUDGET" : "EXHAUSTED");
        fprintf(f, " \"nodes_by_depth\":[");
        for (int d = 0; d <= N; d++) fprintf(f, "%lld%s", nodes_by_depth[d], d < N ? "," : "");
        fprintf(f, "],\n \"covers_containing\":[");
        for (long long i = 0; i < NT; i++) fprintf(f, "%lld%s", cover_count[i], i + 1 < NT ? "," : "");
        fprintf(f, "]}\n");
        if (fclose(f)) { perror(argv[4]); return 1; }
        return timed_out ? 1 : 0;
    }

    if (!strcmp(argv[2], "witness")) {
        if (argc < 7) usage(2);
        long long NQ2;
        unsigned char *qry2 = fdlh_load_records(argv[3], RB, N, &NQ2);
        long long j = atoll(argv[5]);
        if (j < 0 || j >= NQ2) { fprintf(stderr, "query index out of range\n"); return 2; }
        char path[1024];
        snprintf(path, sizeof path, "%s/pool_%lld.bin", argv[4], j);
        long long P;
        unsigned char *pool = fdlh_load_records(path, RB, N, &P);
        M = fdlh_alloc(sizeof(uint64_t) * (size_t)(P ? P : 1) * W, "pool masks");
        for (long long i = 0; i < P; i++) mask_of(pool + i * RB, M + i * W);
        AW = (int)((P + 63) / 64); if (!AW) AW = 1;
        adj = fdlh_calloc((size_t)(P ? P : 1) * AW, sizeof(uint64_t), "adjacency");
        for (long long x = 0; x < P; x++)
            for (long long y = x + 1; y < P; y++)
                if (disjoint(M + (size_t)x * W, M + (size_t)y * W)) {
                    adj[(size_t)x * AW + (y >> 6)] |= (uint64_t)1 << (y & 63);
                    adj[(size_t)y * AW + (x >> 6)] |= (uint64_t)1 << (x & 63);
                }
        best_clique = 0;
        uint64_t *all = fdlh_calloc((size_t)AW, sizeof(uint64_t), "clique candidates");
        for (long long x = 0; x < P; x++) all[x >> 6] |= (uint64_t)1 << (x & 63);
        bk(all, 0);
        FILE *f = fopen(argv[6], "w");
        if (!f) { perror(argv[6]); return 1; }
        fprintf(f, "{\"n\":%d,\"query\":%lld,\"pool\":%lld,\"clique\":%d,"
                   "\"packing\":%d,\"members\":[%lld", N, j, P, best_clique, best_clique + 1, j);
        for (int i = 0; i < best_clique; i++) fprintf(f, ",%d", best_set[i]);
        fprintf(f, "],\n \"squares\":[\n");
        for (int i = 0; i <= best_clique; i++) {
            const unsigned char *r = (i == 0) ? qry2 + j * RB : pool + (size_t)best_set[i - 1] * RB;
            fprintf(f, "  [");
            for (int c = 0; c < RB; c++) fprintf(f, "%d%s", r[c], c + 1 < RB ? "," : "");
            fprintf(f, "]%s\n", i < best_clique ? "," : "");
        }
        fprintf(f, " ]}\n");
        if (fclose(f)) { perror(argv[6]); return 1; }
        printf("{\"query\":%lld,\"pool\":%lld,\"max_companion_clique\":%d,\"packing\":%d}\n",
               j, P, best_clique, best_clique + 1);
        free(all);
        return 0;
    }

    if (strcmp(argv[2], "roots") || argc < 6) usage(2);
    const char *qpath = argv[3], *pooldir = argv[4], *outpath = argv[5];
    int slice_k = 0, slice_w = 1, do_clique = 1, do_search = 1;
    for (int a = 6; a < argc; a++) {
        if (!strcmp(argv[a], "--slice") && a + 2 < argc) { slice_k = atoi(argv[a+1]); slice_w = atoi(argv[a+2]); a += 2; }
        else if (!strcmp(argv[a], "--cap") && a + 1 < argc) { cap_seconds = atof(argv[a+1]); a++; }
        else if (!strcmp(argv[a], "--no-clique")) do_clique = 0;
        else if (!strcmp(argv[a], "--no-search")) do_search = 0;
        else usage(2);
    }
    if (slice_w < 1 || slice_k < 0 || slice_k >= slice_w) { fprintf(stderr, "bad --slice\n"); return 2; }

    long long NQ;
    unsigned char *qry = fdlh_load_records(qpath, RB, N, &NQ);
    FILE *out = fopen(outpath, "w");
    if (!out) { perror(outpath); return 1; }

    for (long long j = 0; j < NQ; j++) {
        if (j % slice_w != slice_k) continue;
        char path[1024];
        snprintf(path, sizeof path, "%s/pool_%lld.bin", pooldir, j);
        long long P;
        unsigned char *pool = fdlh_load_records(path, RB, N, &P);
        uint64_t qm[MAXW];
        mask_of(qry + j * RB, qm);

        /* The deadline covers the whole query: loading, the graph, the search
         * and the clique bound alike. */
        timed_out = 0; node_check = 0; t_start = now_s();

        M = fdlh_alloc(sizeof(uint64_t) * (size_t)(P ? P : 1) * W, "pool masks");
        long long meets = 0;
        for (long long i = 0; i < P; i++) {
            mask_of(pool + i * RB, M + i * W);
            if (!disjoint(M + i * W, qm)) meets++;
        }
        /* A pool member meeting its own query is not a companion.  Searching
         * such a pool would answer a different question, so refuse it. */
        if (meets) {
            fprintf(stderr, "query %lld: %lld of %lld pool members meet the query -- "
                            "this is not a companion pool\n", j, meets, P);
            return 1;
        }

        double w0 = now_s();
        long long covers = 0, nodes = 0;
        int maxd = 0, budget = 0;
        if (do_search) {
            memset(nodes_by_depth, 0, sizeof nodes_by_depth);
            nsol = 0; max_depth_seen = 0;
            cover_count = NULL;
            cand_stride = (int)(P ? P : 1);
            cand_store = fdlh_alloc(sizeof(int) * (size_t)cand_stride * (2 * N + 2), "candidate lists");
            int *surv = fdlh_alloc(sizeof(int) * (size_t)(P ? P : 1), "survivors");
            for (long long i = 0; i < P; i++) surv[i] = (int)i;
            uint64_t cov[MAXW];
            memcpy(cov, qm, sizeof(uint64_t) * (size_t)W);
            chosen[0] = -1;
            search(cov, 1, surv, (int)P);           /* depth 1: the query itself */
            covers = nsol; maxd = max_depth_seen;
            for (int d = 0; d <= N; d++) nodes += nodes_by_depth[d];
            free(surv); free(cand_store); cand_store = NULL;
        }
        double w1 = now_s();

        long long edges = 0, triangles = 0;
        int clique = 0;
        if (do_clique) {
            AW = (int)((P + 63) / 64); if (!AW) AW = 1;
            adj = fdlh_calloc((size_t)(P ? P : 1) * AW, sizeof(uint64_t), "adjacency");
            for (long long a = 0; a < P && !over_budget(1024); a++)
                for (long long b = a + 1; b < P; b++)
                    if (disjoint(M + (size_t)a * W, M + (size_t)b * W)) {
                        adj[(size_t)a * AW + (b >> 6)] |= (uint64_t)1 << (b & 63);
                        adj[(size_t)b * AW + (a >> 6)] |= (uint64_t)1 << (a & 63);
                        edges++;
                    }
            uint64_t *tmp = fdlh_alloc(sizeof(uint64_t) * (size_t)AW, "triangle scratch");
            for (long long a = 0; a < P && !over_budget(1024); a++)
                for (long long b = a + 1; b < P; b++)
                    if (adj[(size_t)a * AW + (b >> 6)] & ((uint64_t)1 << (b & 63))) {
                        int s = 0;
                        for (int u = 0; u < AW; u++) {
                            tmp[u] = adj[(size_t)a * AW + u] & adj[(size_t)b * AW + u];
                            s += __builtin_popcountll(tmp[u]);
                        }
                        triangles += s;
                    }
            triangles /= 3;
            best_clique = 0;
            if (P) {
                uint64_t *all = fdlh_calloc((size_t)AW, sizeof(uint64_t), "clique candidates");
                for (long long a = 0; a < P; a++) all[a >> 6] |= (uint64_t)1 << (a & 63);
                bk(all, 0);
                free(all);
            }
            clique = best_clique;
            free(tmp); free(adj); adj = NULL;
        }
        budget = timed_out;
        double w2 = now_s();

        fprintf(out, "{\"j\":%lld,\"pool\":%lld,\"pool_meets_query\":%lld,"
                     "\"covers\":%lld,\"status\":\"%s\",\"nodes\":%lld,\"max_depth\":%d,"
                     "\"edges\":%lld,\"triangles\":%lld,\"max_companion_clique\":%d,"
                     "\"max_packing_with_query\":%d,\"search_wall\":%.4f,\"clique_wall\":%.4f}\n",
                j, P, meets, covers, budget ? "BUDGET" : (do_search ? "EXHAUSTED" : "SKIPPED"),
                nodes, maxd, edges, triangles, clique, clique + 1, w1 - w0, w2 - w1);
        fflush(out);
        free(pool); free(M); M = NULL;
    }
    if (fclose(out)) { perror(outpath); return 1; }
    return 0;
}
