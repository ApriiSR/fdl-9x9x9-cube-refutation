/* enum.c -- enumerate the supports of [n]^3 as an exact cover.
 *
 * A *support* is a set of cells of [n]^3 meeting every main line exactly once.
 * (Equivalently, the set of cells carrying one symbol in a fully diagonalised
 * Latin cube.)  That is an exact-cover problem verbatim:
 *
 *     items   = the 3n^2 + 6n + 4 main lines   (301 at n = 9)
 *     options = the n^3 cells                  (729 at n = 9)
 *     a cell covers exactly the main lines through it -- the 3 axis lines
 *     always, then 0 to 6 of the face diagonals and 0 to 4 of the space
 *     diagonals, so degrees run 3..13 at n = 9 (the centre (4,4,4) is the
 *     unique degree-13 cell, on all six face diagonals and all four space
 *     diagonals at once).
 *
 * The search is Knuth's Algorithm X with dancing links and minimum-remaining-
 * values item selection.  Nothing here knows that a support is the graph of a
 * Latin square: the n^2 cells and the Latin structure come out of the cover,
 * they are not put in.
 *
 * SHARDING.  The n cells with i = 0 are forced to be one per (0,j) z-line and
 * one per symbol, and to meet the two face diagonals of the plane i = 0 once
 * each; so row 0 of a support is a permutation p of [n] with exactly one fixed
 * point and exactly one reflected point (p(j) = n-1-j).  Fixing that
 * permutation splits the problem into one independent shard per admissible p.
 * The shard file produced by `shards` lists them.
 *
 * OUTPUT.  n*n bytes per support: byte i*n+j is k for the support's cell
 * (i, j, k).  Records within a shard are sorted lexicographically, and shards
 * are meant to be concatenated in increasing shard index, so the catalogue is
 * a deterministic function of the mathematics and not of the search order.
 *
 * Build: cc -O2 -o enum enum.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

#include "lines.h"
#include "util.h"

#define MAXNODES (MAXLINES + 1 + MAXDEG * MAXCELLS)

static int N, NI, NC;

/* dancing links */
static int Lk[MAXNODES], Rk[MAXNODES], Uk[MAXNODES], Dk[MAXNODES];
static int Col[MAXNODES], Row[MAXNODES];
static int Sz[MAXLINES + 1];
static int nnode, root;
static int cellhead[MAXCELLS];
static int lines[MAXLINES][MAXN];
static int cellitems[MAXCELLS][MAXDEG];
static int cellndeg[MAXCELLS];

static void build_matrix(void)
{
    NC = N * N * N;
    root = NI;
    for (int c = 0; c <= NI; c++) {
        Lk[c] = (c == 0) ? root : c - 1;
        Rk[c] = (c == root) ? 0 : c + 1;
        Uk[c] = Dk[c] = c;
        Col[c] = c; Row[c] = -1; Sz[c] = 0;
    }
    Lk[0] = root; Rk[root] = 0;
    Lk[root] = NI - 1; Rk[NI - 1] = root;
    nnode = NI + 1;

    for (int c = 0; c < NC; c++) cellndeg[c] = 0;
    for (int l = 0; l < NI; l++)
        for (int t = 0; t < N; t++) {
            int c = lines[l][t];
            if (cellndeg[c] >= MAXDEG) {
                fprintf(stderr, "cell degree exceeds MAXDEG (%d)\n", MAXDEG);
                exit(1);
            }
            cellitems[c][cellndeg[c]++] = l;
        }
    for (int c = 0; c < NC; c++) {
        int first = -1, prev = -1;
        for (int d = 0; d < cellndeg[c]; d++) {
            int l = cellitems[c][d], x = nnode++;
            Col[x] = l; Row[x] = c;
            Uk[x] = Uk[l]; Dk[x] = l; Dk[Uk[l]] = x; Uk[l] = x;
            Sz[l]++;
            if (first < 0) { first = x; Lk[x] = Rk[x] = x; }
            else { Lk[x] = prev; Rk[x] = first; Rk[prev] = x; Lk[first] = x; }
            prev = x;
        }
        cellhead[c] = first;
    }
}

static void cover(int c)
{
    Rk[Lk[c]] = Rk[c]; Lk[Rk[c]] = Lk[c];
    for (int i = Dk[c]; i != c; i = Dk[i])
        for (int j = Rk[i]; j != i; j = Rk[j]) {
            Uk[Dk[j]] = Uk[j]; Dk[Uk[j]] = Dk[j];
            Sz[Col[j]]--;
        }
}
static void uncover(int c)
{
    for (int i = Uk[c]; i != c; i = Uk[i])
        for (int j = Lk[i]; j != i; j = Lk[j]) {
            Sz[Col[j]]++;
            Uk[Dk[j]] = j; Dk[Uk[j]] = j;
        }
    Rk[Lk[c]] = c; Lk[Rk[c]] = c;
}
static void select_cell(int c)
{
    int r = cellhead[c];
    cover(Col[r]);
    for (int j = Rk[r]; j != r; j = Rk[j]) cover(Col[j]);
}
static void deselect_cell(int c)
{
    int r = cellhead[c];
    for (int j = Lk[r]; j != r; j = Lk[j]) uncover(Col[j]);
    uncover(Col[r]);
}

/* ---- search ------------------------------------------------------------- */
static int sol[MAXCELLS], fixed_cells[MAXN], nfixed;
static long long nsol, nnodes;
static int rec_bytes;
static int do_collect;
static unsigned char *recs;         /* collected records, sorted before write */
static size_t nrecs, caprecs;

static double cap_seconds;          /* 0 = no cap */
static double t_start;
static int timed_out;
static long long node_check;

#define now_s fdlh_now_s

static void emit(int depth)
{
    nsol++;
    if (!do_collect) return;
    if (nrecs == caprecs) {
        caprecs = caprecs ? caprecs * 2 : 4096;
        unsigned char *grown = realloc(recs, caprecs * (size_t)rec_bytes);
        if (!grown) { fprintf(stderr, "out of memory collecting records\n"); exit(1); }
        recs = grown;
    }
    unsigned char *p = recs + nrecs * (size_t)rec_bytes;
    memset(p, 0xff, rec_bytes);
    for (int a = 0; a < nfixed; a++) { int c = fixed_cells[a]; p[c / N] = (unsigned char)(c % N); }
    for (int a = 0; a < depth; a++)  { int c = sol[a];         p[c / N] = (unsigned char)(c % N); }
    nrecs++;
}

static void search(int depth)
{
    if (timed_out) return;
    if (Rk[root] == root) { emit(depth); return; }
    if (cap_seconds > 0.0 && ++node_check >= 256) {
        node_check = 0;
        if (now_s() - t_start > cap_seconds) { timed_out = 1; return; }
    }
    nnodes++;
    int best = -1, bs = 1 << 30;
    for (int c = Rk[root]; c != root; c = Rk[c])
        if (Sz[c] < bs) { bs = Sz[c]; best = c; if (bs <= 1) break; }
    if (bs == 0) return;
    cover(best);
    for (int r = Dk[best]; r != best; r = Dk[r]) {
        sol[depth] = Row[r];
        for (int j = Rk[r]; j != r; j = Rk[j]) cover(Col[j]);
        search(depth + 1);
        for (int j = Lk[r]; j != r; j = Lk[j]) uncover(Col[j]);
        if (timed_out) break;
    }
    uncover(best);
}

static int reccmp(const void *a, const void *b)
{
    return memcmp(a, b, (size_t)rec_bytes);
}

/* lexicographic on the cell indices as numbers, not as bytes */
static int linecmp(const int *a, const int *b)
{
    for (int t = 0; t < N; t++) if (a[t] != b[t]) return a[t] < b[t] ? -1 : 1;
    return 0;
}

/* row 0 must be a permutation with exactly one fixed and one reflected point */
static int admissible_row0(const int *p)
{
    int f = 0, a = 0, seen[MAXN];
    memset(seen, 0, sizeof seen);
    for (int j = 0; j < N; j++) {
        if (p[j] < 0 || p[j] >= N || seen[p[j]]) return 0;
        seen[p[j]] = 1;
        if (p[j] == j) f++;
        if (p[j] == N - 1 - j) a++;
    }
    return f == 1 && a == 1;
}

/* Returns 0 exhausted, 1 hit the time cap, 2 inadmissible row 0. */
static int run_one(const int *p, int collect, long long *count,
                   long long *nodes, double *wall)
{
    build_matrix();
    nfixed = 0; nsol = 0; nnodes = 0; timed_out = 0; node_check = 255;
    nrecs = 0; do_collect = collect;
    if (p) {
        if (!admissible_row0(p)) return 2;
        for (int j = 0; j < N; j++) fixed_cells[nfixed++] = (0 * N + j) * N + p[j];
    }
    t_start = now_s();
    for (int a = 0; a < nfixed; a++) select_cell(fixed_cells[a]);
    search(0);
    for (int a = nfixed - 1; a >= 0; a--) deselect_cell(fixed_cells[a]);
    *wall = now_s() - t_start;
    *count = nsol; *nodes = nnodes;
    if (collect && !timed_out) qsort(recs, nrecs, rec_bytes, reccmp);
    return timed_out ? 1 : 0;
}

static int write_records(const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); return 1; }
    if (nrecs && fwrite(recs, rec_bytes, nrecs, f) != nrecs) { perror(path); fclose(f); return 1; }
    if (fflush(f) || fclose(f)) { perror(path); return 1; }
    return 0;
}

static void usage(int rc)
{
    fprintf(stderr,
"enum -- enumerate the supports of [n]^3 (exact cover over the main lines)\n"
"\n"
"usage:\n"
"  enum <n> lines [--dump FILE]\n"
"        Report the main-line count and the cell-degree range, and check that\n"
"        every line has n distinct in-range cells.  Exits nonzero on failure.\n"
"        --dump writes the incidence structure in canonical form -- one line\n"
"        per main line, its cells in increasing order, the lines sorted -- so\n"
"        that an independent construction can be compared as a SET of lines\n"
"        rather than as a count.\n"
"  enum <n> all <out.bin|->\n"
"        Enumerate every support of [n]^3 (no sharding).  Records are sorted.\n"
"  enum <n> one <out.bin|-> p0 p1 ... p<n-1>\n"
"        Enumerate one shard, given row 0 explicitly.\n"
"  enum <n> shards <shardfile> <outdir> <manifest.jsonl> [options]\n"
"        Enumerate a list of shards, one payload file per shard, appending one\n"
"        JSON line per finished shard to the manifest.  Resumable: a shard is\n"
"        skipped only if the manifest records it done AND its payload is still\n"
"        on disk at the recorded length and digest; otherwise it is redone and\n"
"        the reason is printed.  A torn final manifest line -- a kill during a\n"
"        write -- is discarded before appending, and a stale .part file is\n"
"        removed before the shard is redone.\n"
"        shardfile lines are: idx p0 p1 ... p<n-1>   (see the `shards` tool)\n"
"        options:\n"
"          --slice K W    take only shards whose position in the file is K mod W\n"
"          --done FILE    additionally read already-finished shards from FILE\n"
"                         (so a sweep can resume with a different worker count:\n"
"                          merge the workers' manifests and pass the merge)\n"
"          --cap SECONDS  give up on a shard after SECONDS (records BUDGET,\n"
"                         leaves no payload, and does not mark it done)\n"
"\n"
"A payload file holds n*n bytes per support: byte i*n+j is k for the cell\n"
"(i,j,k).  Concatenating payloads in increasing shard index gives the catalogue.\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 3) usage(2);
    N = atoi(argv[1]);
    if (N < 0 || N > MAXN) { fprintf(stderr, "n out of range 0..%d\n", MAXN); return 2; }
    NI = fdlh_build_lines(N, lines);
    rec_bytes = N * N;
    const char *mode = argv[2];

    if (!strcmp(mode, "lines")) {
        int expect = fdlh_expected_lines(N);
        build_matrix();
        long long inc = 0; int mind = N ? 1 << 30 : 0, maxd = 0;
        for (int c = 0; c < N * N * N; c++) {
            inc += cellndeg[c];
            if (cellndeg[c] < mind) mind = cellndeg[c];
            if (cellndeg[c] > maxd) maxd = cellndeg[c];
        }
        printf("n=%d lines=%d expected=%d\n", N, NI, expect);
        printf("cells=%d incidences=%lld (expect %d) degree min=%d max=%d\n",
               N * N * N, inc, NI * N, mind, maxd);
        static char seen[MAXCELLS];
        for (int l = 0; l < NI; l++) {
            memset(seen, 0, (size_t)N * N * N);
            for (int t = 0; t < N; t++) {
                int c = lines[l][t];
                if (c < 0 || c >= N * N * N || seen[c]) { printf("BAD line %d\n", l); return 1; }
                seen[c] = 1;
            }
        }
        printf("every line has %d distinct in-range cells\n", N);
        /* The canonical form of the incidence structure: each line as its cells
         * in increasing order, the lines themselves in increasing order.  Two
         * implementations that agree here agree as line SETS, which equal counts
         * would not establish. */
        if (argc >= 5 && !strcmp(argv[3], "--dump")) {
            int *ord = fdlh_alloc(sizeof(int) * (size_t)NI, "line order");
            for (int l = 0; l < NI; l++) ord[l] = l;
            for (int l = 0; l < NI; l++) {          /* cells of each line, sorted */
                for (int a = 1; a < N; a++) {
                    int v = lines[l][a], b = a - 1;
                    while (b >= 0 && lines[l][b] > v) { lines[l][b + 1] = lines[l][b]; b--; }
                    lines[l][b + 1] = v;
                }
            }
            for (int a = 1; a < NI; a++) {          /* lines, lexicographically */
                int v = ord[a], b = a - 1;
                while (b >= 0 && linecmp(lines[ord[b]], lines[v]) > 0) {
                    ord[b + 1] = ord[b]; b--;
                }
                ord[b + 1] = v;
            }
            FILE *d = fopen(argv[4], "w");
            if (!d) { perror(argv[4]); return 1; }
            for (int a = 0; a < NI; a++) {
                for (int t = 0; t < N; t++) fprintf(d, "%d%s", lines[ord[a]][t], t + 1 < N ? " " : "\n");
            }
            if (fclose(d)) { perror(argv[4]); return 1; }
            free(ord);
        }
        return (NI == expect && inc == (long long)NI * N) ? 0 : 1;
    }

    if (!strcmp(mode, "all")) {
        if (argc < 4) usage(2);
        int collect = strcmp(argv[3], "-") != 0;
        long long cnt, nds; double w;
        run_one(NULL, collect, &cnt, &nds, &w);
        if (collect && write_records(argv[3])) return 1;
        printf("n=%d unsharded: %lld supports, %lld nodes, %.3f s\n", N, cnt, nds, w);
        return 0;
    }

    if (!strcmp(mode, "one")) {
        if (argc < 4 + N) usage(2);
        int collect = strcmp(argv[3], "-") != 0;
        int p[MAXN];
        for (int j = 0; j < N; j++) p[j] = atoi(argv[4 + j]);
        long long cnt, nds; double w;
        int rc = run_one(p, collect, &cnt, &nds, &w);
        if (rc == 2) { fprintf(stderr, "row 0 is not admissible\n"); return 3; }
        if (rc == 0 && collect && write_records(argv[3])) return 1;
        printf("{\"count\":%lld,\"nodes\":%lld,\"wall\":%.4f,\"status\":\"%s\",\"done\":%s}\n",
               cnt, nds, w, rc == 1 ? "BUDGET" : "EXHAUSTED", rc == 1 ? "false" : "true");
        return 0;
    }

    if (!strcmp(mode, "shards")) {
        if (argc < 6) usage(2);
        const char *shardfile = argv[3], *outdir = argv[4], *manifest = argv[5];
        int slice_k = 0, slice_w = 1;
        const char *donefile = NULL;
        for (int a = 6; a < argc; a++) {
            if (!strcmp(argv[a], "--slice") && a + 2 < argc) { slice_k = atoi(argv[a+1]); slice_w = atoi(argv[a+2]); a += 2; }
            else if (!strcmp(argv[a], "--cap") && a + 1 < argc) { cap_seconds = atof(argv[a+1]); a += 1; }
            else if (!strcmp(argv[a], "--done") && a + 1 < argc) { donefile = argv[a+1]; a += 1; }
            else usage(2);
        }
        if (slice_w < 1 || slice_k < 0 || slice_k >= slice_w) { fprintf(stderr, "bad --slice\n"); return 2; }

        /* Resume: shards already marked done, in this worker's own manifest and
         * in the optional merged one.  A completion record is trusted only if
         * it has the exact shape a writer produces AND its payload is still on
         * disk at the recorded length and digest; anything else is redone. */
        struct fdlh_done done;
        fdlh_done_init(&done, 1 << 17);
        fdlh_done_read(&done, manifest);
        if (donefile) fdlh_done_read(&done, donefile);
        mkdir(outdir, 0777);
        FILE *sf = fopen(shardfile, "r");
        if (!sf) { perror(shardfile); return 1; }
        FILE *out_mf = fdlh_manifest_append(manifest);
        if (!out_mf) { perror(manifest); return 1; }

        char line[4096];
        long long seen = 0;
        while (fgets(line, sizeof line, sf)) {
            if (line[0] == '#' || line[0] == '\n') continue;
            int idx, p[MAXN], ok = 1;
            char *tok = strtok(line, " \t\n");
            if (!tok) continue;
            idx = atoi(tok);
            for (int j = 0; j < N; j++) {
                tok = strtok(NULL, " \t\n");
                if (!tok) { ok = 0; break; }
                p[j] = atoi(tok);
            }
            if (!ok) continue;
            if ((seen++ % slice_w) != slice_k) continue;

            char sub[1024], path[1200], tmp[1300];
            snprintf(sub, sizeof sub, "%s/%03d", outdir, idx / 1000);
            snprintf(path, sizeof path, "%s/s%05d.bin", sub, idx);
            snprintf(tmp, sizeof tmp, "%s.part", path);

            long long dbytes; uint64_t ddig; int dhas;
            if (fdlh_done_get(&done, idx, &dbytes, &ddig, &dhas)) {
                const char *why = "the completion record has no length or digest";
                if (dhas && fdlh_payload_ok(path, dbytes, ddig, &why)) continue;
                fprintf(stderr, "shard %d: %s -- regenerating it\n", idx, why);
            }
            mkdir(sub, 0777);
            unlink(tmp);            /* a .part left by an interrupted run */

            long long cnt, nds; double w;
            int rc = run_one(p, 1, &cnt, &nds, &w);
            if (rc == 2) {
                fprintf(out_mf, "{\"idx\":%d,\"status\":\"BADROW0\",\"done\":false}\n", idx);
                fflush(out_mf); continue;
            }
            if (rc == 1) {   /* time cap: leave nothing behind, do not mark done */
                unlink(tmp);
                fprintf(out_mf, "{\"idx\":%d,\"status\":\"BUDGET\",\"done\":false,"
                                "\"count\":%lld,\"wall\":%.4f}\n", idx, cnt, w);
                fflush(out_mf); continue;
            }
            if (write_records(tmp)) return 1;
            if (rename(tmp, path)) { perror("rename"); return 1; }
            uint64_t dig = fdlh_fnv64(recs, nrecs * (size_t)rec_bytes);
            fprintf(out_mf, "{\"idx\":%d,\"row0\":[", idx);
            for (int j = 0; j < N; j++) fprintf(out_mf, "%d%s", p[j], j + 1 < N ? "," : "");
            fprintf(out_mf, "],\"count\":%lld,\"nodes\":%lld,\"wall\":%.4f,\"bytes\":%lld,"
                            "\"digest\":\"%016llx\",\"status\":\"EXHAUSTED\",\"done\":true}\n",
                    cnt, nds, w, cnt * (long long)rec_bytes, (unsigned long long)dig);
            fflush(out_mf);
        }
        if (fclose(out_mf)) { perror(manifest); return 1; }
        fclose(sf);
        return 0;
    }

    usage(2);
    return 2;
}
