/* symmetry.c -- the subgroup of the cell group that fixes the plane x = 0,
 * and the orbits it induces on the shard universe.
 *
 * The enumeration is sharded by row 0, that is by the support's intersection
 * with the plane x = 0.  An element of the cell group
 *
 *     g(x)_i = tau( eps_i( x_{pi(i)} ) ),      pi in S_3,
 *                                              eps in {id, rev}^3,
 *                                              tau in C_{S_n}(rev),
 *
 * therefore respects the sharding exactly when it maps that plane to itself.
 * Which elements do?  The image's first coordinate is tau(eps_0(x_{pi(0)})).
 * If pi(0) != 0 it varies with a coordinate that is free on the plane, so it
 * cannot be constant; hence pi(0) = 0, and then the condition is
 * tau(eps_0(0)) = 0.  Write H for the set of such elements.  It is the setwise
 * stabiliser of a plane, hence a subgroup, and
 *
 *     |H| = |G| / (3 * 2*floor(n/2)) = 9216 / 24 = 384   at n = 8 and n = 9,
 *
 * the index being the number of planes {x_a = c} in the G-orbit of {x_0 = 0}:
 * three choices of axis, and c = tau(0) or tau(n-1) ranges over everything
 * except the middle symbol of an odd n.
 *
 * On the plane, such an h acts as (0,j,k) -> (0, u(j), v(k)) when pi fixes the
 * last two coordinates and as (0,j,k) -> (0, u(k), v(j)) when pi swaps them,
 * with u = tau . eps_1 and v = tau . eps_2.  So a shard's row 0, the
 * permutation p, is carried to v.p.u^{-1} or to v.p^{-1}.u^{-1}: the image
 * depends on p alone, and H permutes the admissible permutations among
 * themselves.  That is Lemma 6 of the README, and this program is its
 * computational half: `orbits` checks the closure explicitly and `expand`
 * produces each shard's contents as the image of its orbit representative's.
 *
 * Records are n*n bytes: byte i*n+j is k for the cell (i,j,k), the same format
 * `enum` writes, and images are sorted lexicographically so a mapped shard is
 * byte-identical to a directly enumerated one.
 *
 * Build: cc -O2 -o symmetry symmetry.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

#define YMAXN 12
#define MAXG 20000
#define MAXSHARDS 65536

static int N, NC, RB;                 /* side, n^3, record bytes = n^2 */
static int *G; static int NG;         /* the whole cell group */
static int *H; static int NH;         /* the plane-fixing subgroup */

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static uint64_t fnv(const void *p, size_t nbytes)
{
    const unsigned char *b = p;
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < nbytes; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

/* ---- C(rev), then the cell group ---------------------------------------- */
static int taus[512][YMAXN];
static int ntau;

static void build_taus(void)
{
    int h = N / 2, perm[YMAXN];
    ntau = 0;
    long long nperm = 1;
    for (int i = 2; i <= h; i++) nperm *= i;
    for (long long pc = 0; pc < nperm; pc++) {
        int avail[YMAXN], na = h;
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
                int b = perm[a], d0 = b, d1 = N - 1 - b;
                if (flips & (1 << a)) { int s = d0; d0 = d1; d1 = s; }
                t[a] = d0; t[N - 1 - a] = d1;
            }
            if (N % 2) t[N / 2] = N / 2;
            ntau++;
            if (ntau > 512) { fprintf(stderr, "too many tau\n"); exit(1); }
        }
    }
}

static void build_group(void)
{
    build_taus();
    G = malloc(sizeof(int) * (size_t)MAXG * NC);
    if (!G) { fprintf(stderr, "out of memory\n"); exit(1); }
    NG = 0;
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
        size_t s = fnv(g, sizeof(int) * (size_t)NC) & (hs - 1);
        int dup = 0;
        while (tab[s] >= 0) {
            if (!memcmp(G + (size_t)tab[s] * NC, g, sizeof(int) * (size_t)NC)) { dup = 1; break; }
            s = (s + 1) & (hs - 1);
        }
        if (dup) continue;
        if (NG >= MAXG) { fprintf(stderr, "group too large\n"); exit(1); }
        memcpy(G + (size_t)NG * NC, g, sizeof(int) * (size_t)NC);
        tab[s] = NG++;
    }
    free(g); free(tab);
}

/* h fixes the plane x = 0 setwise iff every cell of that plane has an image in
 * it; the plane is exactly the first n^2 cell indices. */
static int fixes_plane(const int *g)
{
    for (int c = 0; c < RB; c++) if (g[c] >= RB) return 0;
    return 1;
}

static int h_identity;

static void build_subgroup(void)
{
    H = malloc(sizeof(int) * (size_t)NG * NC);
    if (!H) { fprintf(stderr, "out of memory\n"); exit(1); }
    NH = 0;
    for (int i = 0; i < NG; i++)
        if (fixes_plane(G + (size_t)i * NC))
            memcpy(H + (size_t)NH++ * NC, G + (size_t)i * NC, sizeof(int) * (size_t)NC);
    h_identity = -1;
    for (int i = 0; i < NH && h_identity < 0; i++) {
        int ok = 1;
        for (int c = 0; c < NC; c++) if (H[(size_t)i * NC + c] != c) { ok = 0; break; }
        if (ok) h_identity = i;
    }
    if (h_identity < 0) { fprintf(stderr, "the subgroup does not contain the identity\n"); exit(1); }
}

/* ---- the shard universe -------------------------------------------------- */
static int nshard;
static unsigned char rowbytes[MAXSHARDS][YMAXN];
static int shard_idx[MAXSHARDS];       /* file order -> shard index */
static int *rowtab; static size_t rowmask;   /* row bytes -> position */

static int admissible(const unsigned char *p)
{
    int f = 0, a = 0, seen[YMAXN];
    memset(seen, 0, sizeof seen);
    for (int j = 0; j < N; j++) {
        if (p[j] >= N || seen[p[j]]) return 0;
        seen[p[j]] = 1;
        if (p[j] == j) f++;
        if (p[j] == N - 1 - j) a++;
    }
    return f == 1 && a == 1;
}

static void index_rows(void)
{
    size_t hs = 1;
    while (hs < (size_t)nshard * 4) hs <<= 1;
    rowmask = hs - 1;
    rowtab = malloc(sizeof(int) * hs);
    for (size_t i = 0; i < hs; i++) rowtab[i] = -1;
    for (int i = 0; i < nshard; i++) {
        size_t s = fnv(rowbytes[i], (size_t)N) & rowmask;
        while (rowtab[s] >= 0) {
            if (!memcmp(rowbytes[rowtab[s]], rowbytes[i], (size_t)N)) {
                fprintf(stderr, "the shard file repeats a row\n"); exit(1);
            }
            s = (s + 1) & rowmask;
        }
        rowtab[s] = i;
    }
}

static int find_row(const unsigned char *p)
{
    size_t s = fnv(p, (size_t)N) & rowmask;
    while (rowtab[s] >= 0) {
        if (!memcmp(rowbytes[rowtab[s]], p, (size_t)N)) return rowtab[s];
        s = (s + 1) & rowmask;
    }
    return -1;
}

static void read_shardfile(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    char line[4096];
    nshard = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *tok = strtok(line, " \t\n");
        if (!tok) continue;
        if (nshard >= MAXSHARDS) { fprintf(stderr, "too many shards\n"); exit(1); }
        shard_idx[nshard] = atoi(tok);
        for (int j = 0; j < N; j++) {
            tok = strtok(NULL, " \t\n");
            if (!tok) { fprintf(stderr, "%s: short line\n", path); exit(1); }
            rowbytes[nshard][j] = (unsigned char)atoi(tok);
        }
        if (!admissible(rowbytes[nshard])) {
            fprintf(stderr, "%s: row %d is not admissible\n", path, shard_idx[nshard]);
            exit(1);
        }
        nshard++;
    }
    fclose(f);
    index_rows();
}

/* the image of a row under h, read off the plane cell by cell */
static void row_image(const int *g, const unsigned char *p, unsigned char *out)
{
    int seen[YMAXN];
    memset(seen, 0, sizeof seen);
    memset(out, 0xff, (size_t)N);
    for (int j = 0; j < N; j++) {
        int cell = (0 * N + j) * N + p[j];
        int img = g[cell];
        if (img >= RB) { fprintf(stderr, "a subgroup element left the plane\n"); exit(1); }
        int jj = img / N, kk = img % N;
        if (seen[jj]) { fprintf(stderr, "a row image is not a function\n"); exit(1); }
        seen[jj] = 1;
        out[jj] = (unsigned char)kk;
    }
}

/* ---- records ------------------------------------------------------------- */
static void record_image(const unsigned char *r, const int *g, unsigned char *out)
{
    for (int c = 0; c < RB; c++) {
        int cell = c * N + r[c];
        int img = g[cell];
        out[img / N] = (unsigned char)(img % N);
    }
}

static int reccmp(const void *a, const void *b)
{
    return memcmp(a, b, (size_t)RB);
}

static void shard_path(char *buf, size_t nbuf, const char *dir, int idx)
{
    snprintf(buf, nbuf, "%s/%03d/s%05d.bin", dir, idx / 1000, idx);
}

/* ---- modes --------------------------------------------------------------- */
static int cmd_group(void)
{
    double t0 = now_s();
    build_group();
    build_subgroup();
    /* the subgroup really is one: closed under composition, and closed under
     * inverses because a finite closed non-empty set is */
    int *comp = malloc(sizeof(int) * (size_t)NC);
    long long outside = 0;
    for (int a = 0; a < NH; a++)
    for (int b = 0; b < NH; b++) {
        const int *ga = H + (size_t)a * NC, *gb = H + (size_t)b * NC;
        for (int c = 0; c < NC; c++) comp[c] = gb[ga[c]];
        if (!fixes_plane(comp)) outside++;
    }
    free(comp);
    /* the planes {x_a = c} the whole group can send {x_0 = 0} to */
    int planes = 0;
    for (int a = 0; a < 3; a++)
    for (int c = 0; c < N; c++) {
        for (int i = 0; i < NG; i++) {
            const int *g = G + (size_t)i * NC;
            int ok = 1;
            for (int j = 0; j < RB && ok; j++) {
                int img = g[j], x[3];
                x[2] = img % N; x[1] = (img / N) % N; x[0] = img / (N * N);
                if (x[a] != c) ok = 0;
            }
            if (ok) { planes++; break; }
        }
    }
    int expect = 3 * 2 * (N / 2);
    printf("n=%d |G|=%d |H|=%d index=%d (%.1f s)\n", N, NG, NH, NG / NH, now_s() - t0);
    printf("planes in the orbit of {x=0}: %d -- expected 3*2*floor(n/2) = %d -- %s\n",
           planes, expect, planes == expect ? "ok" : "MISMATCH");
    printf("expected |H| = |G|/%d = %d -- %s\n", expect, NG / expect,
           NH == NG / expect ? "ok" : "MISMATCH");
    printf("products of subgroup elements leaving the plane: %lld\n", outside);
    return (planes == expect && NH == NG / expect && outside == 0) ? 0 : 1;
}

static int cmd_orbits(const char *shardfile, const char *outjson, const char *repsfile)
{
    double t0 = now_s();
    build_group();
    build_subgroup();
    read_shardfile(shardfile);

    int *rep = malloc(sizeof(int) * (size_t)nshard);
    int *gel = malloc(sizeof(int) * (size_t)nshard);
    int *orb = malloc(sizeof(int) * (size_t)nshard);
    for (int i = 0; i < nshard; i++) { rep[i] = -1; gel[i] = -1; orb[i] = -1; }

    /* Shards in increasing index; the representative is the least-indexed
     * member of its orbit, so the choice depends on nothing but the group. */
    int *byidx = malloc(sizeof(int) * (size_t)nshard);
    for (int i = 0; i < nshard; i++) byidx[i] = i;
    for (int a = 1; a < nshard; a++) {          /* insertion sort on shard index */
        int v = byidx[a], b = a - 1;
        while (b >= 0 && shard_idx[byidx[b]] > shard_idx[v]) { byidx[b + 1] = byidx[b]; b--; }
        byidx[b + 1] = v;
    }

    unsigned char img[YMAXN];
    int norb = 0;
    long long *osize = malloc(sizeof(long long) * (size_t)nshard);
    for (int a = 0; a < nshard; a++) {
        int i = byidx[a];
        if (rep[i] >= 0) continue;
        int o = norb++;
        osize[o] = 0;
        for (int k = 0; k < NH; k++) {
            row_image(H + (size_t)k * NC, rowbytes[i], img);
            if (!admissible(img)) {
                fprintf(stderr, "a subgroup element maps an admissible row to an "
                                "inadmissible one -- shard %d, element %d\n", shard_idx[i], k);
                return 1;
            }
            int j = find_row(img);
            if (j < 0) {
                fprintf(stderr, "the shard universe is not closed under the subgroup "
                                "-- shard %d, element %d\n", shard_idx[i], k);
                return 1;
            }
            if (rep[j] < 0) { rep[j] = i; gel[j] = k; orb[j] = o; osize[o]++; }
        }
        if (rep[i] != i) {
            fprintf(stderr, "shard %d is not its own representative\n", shard_idx[i]);
            return 1;
        }
    }
    for (int i = 0; i < nshard; i++)
        if (rep[i] < 0) { fprintf(stderr, "shard %d was reached by nothing\n", shard_idx[i]); return 1; }

    FILE *jo = fopen(outjson, "w");
    if (!jo) { perror(outjson); return 1; }
    FILE *rf = fopen(repsfile, "w");
    if (!rf) { perror(repsfile); return 1; }
    long long nreps = 0;
    for (int a = 0; a < nshard; a++) {
        int i = byidx[a];
        fprintf(jo, "{\"idx\":%d,\"orbit\":%d,\"rep\":%d,\"g\":%d,\"row0\":[",
                shard_idx[i], orb[i], shard_idx[rep[i]], gel[i]);
        for (int j = 0; j < N; j++) fprintf(jo, "%d%s", rowbytes[i][j], j + 1 < N ? "," : "");
        fprintf(jo, "]}\n");
        if (rep[i] == i) {
            fprintf(rf, "%d", shard_idx[i]);
            for (int j = 0; j < N; j++) fprintf(rf, " %d", rowbytes[i][j]);
            fputc('\n', rf);
            nreps++;
        }
    }
    if (fclose(jo) || fclose(rf)) { perror("output"); return 1; }

    /* orbit-size histogram; every size must divide |H| */
    long long sizes[64], counts[64]; int nsz = 0, bad = 0;
    long long tot = 0;
    for (int o = 0; o < norb; o++) {
        tot += osize[o];
        if (NH % osize[o]) bad++;
        int f = -1;
        for (int s = 0; s < nsz; s++) if (sizes[s] == osize[o]) { f = s; break; }
        if (f < 0) { if (nsz >= 64) { fprintf(stderr, "too many distinct orbit sizes\n"); return 1; }
                     sizes[nsz] = osize[o]; counts[nsz] = 0; f = nsz++; }
        counts[f]++;
    }
    for (int a = 1; a < nsz; a++) {
        long long sv = sizes[a], cv = counts[a]; int b = a - 1;
        while (b >= 0 && sizes[b] > sv) { sizes[b+1] = sizes[b]; counts[b+1] = counts[b]; b--; }
        sizes[b+1] = sv; counts[b+1] = cv;
    }
    printf("{\"shards\":%d,\"subgroup_order\":%d,\"orbits\":%d,\"representatives\":%lld,"
           "\"sum_orbit_sizes\":%lld,\"reduction\":%.2f,\"orbit_size_histogram\":{",
           nshard, NH, norb, nreps, tot, nshard / (double)norb);
    for (int s = 0; s < nsz; s++)
        printf("%s\"%lld\":%lld", s ? "," : "", sizes[s], counts[s]);
    printf("},\"wall\":%.1f}\n", now_s() - t0);
    if (tot != nshard || nreps != norb || bad) {
        fprintf(stderr, "orbit bookkeeping is inconsistent\n");
        return 1;
    }
    return 0;
}

struct job { int idx, rep, g; unsigned char row[YMAXN]; };

static int cmd_expand(const char *orbfile, const char *sharddir, const char *manifest,
                      int slice_k, int slice_w, const char *donefile)
{
    double t0 = now_s();
    build_group();
    build_subgroup();

    unsigned char *done = calloc(1 << 20, 1);
    for (int pass = 0; pass < 2; pass++) {
        const char *src = pass ? donefile : manifest;
        if (!src) continue;
        FILE *mf = fopen(src, "r");
        if (!mf) continue;
        char line[4096];
        while (fgets(line, sizeof line, mf)) {
            char *q = strstr(line, "\"idx\":");
            if (q && strstr(line, "\"done\":true")) {
                int v = atoi(q + 6);
                if (v >= 0 && v < (1 << 20)) done[v] = 1;
            }
        }
        fclose(mf);
    }

    FILE *of = fopen(orbfile, "r");
    if (!of) { perror(orbfile); return 1; }
    FILE *out_mf = fopen(manifest, "a");
    if (!out_mf) { perror(manifest); return 1; }

    char line[8192];
    long long seen = 0, made = 0, records = 0;
    int cached_rep = -1;
    unsigned char *src = NULL, *dst = NULL;
    size_t nsrc = 0, capdst = 0;

    while (fgets(line, sizeof line, of)) {
        struct job J;
        char *q;
        if (!(q = strstr(line, "\"idx\":"))) continue;
        J.idx = atoi(q + 6);
        if (!(q = strstr(line, "\"rep\":"))) continue;
        J.rep = atoi(q + 6);
        if (!(q = strstr(line, "\"g\":"))) continue;
        J.g = atoi(q + 4);
        if (!(q = strstr(line, "\"row0\""))) continue;
        if (!(q = strchr(q, '['))) continue;
        q++;
        for (int j = 0; j < N; j++) { J.row[j] = (unsigned char)atoi(q); q = strchr(q, ','); if (q) q++; else if (j + 1 < N) return 1; }
        if (J.rep == J.idx) continue;                 /* enumerated directly */
        if ((seen++ % slice_w) != slice_k) continue;
        if (J.idx >= 0 && J.idx < (1 << 20) && done[J.idx]) continue;
        if (J.g < 0 || J.g >= NH) { fprintf(stderr, "shard %d: bad element index\n", J.idx); return 1; }

        if (J.rep != cached_rep) {
            char p[1200];
            shard_path(p, sizeof p, sharddir, J.rep);
            FILE *f = fopen(p, "rb");
            if (!f) { perror(p); return 1; }
            fseek(f, 0, SEEK_END); long long sz = ftell(f); fseek(f, 0, SEEK_SET);
            if (sz % RB) { fprintf(stderr, "%s: not a multiple of %d\n", p, RB); return 1; }
            free(src);
            src = malloc((size_t)sz ? (size_t)sz : 1);
            if (!src) { fprintf(stderr, "out of memory\n"); return 1; }
            if (sz && fread(src, 1, (size_t)sz, f) != (size_t)sz) { perror(p); return 1; }
            fclose(f);
            nsrc = (size_t)(sz / RB);
            cached_rep = J.rep;
        }
        if (nsrc > capdst) {
            free(dst);
            capdst = nsrc;
            dst = malloc(capdst ? capdst * (size_t)RB : 1);
            if (!dst) { fprintf(stderr, "out of memory\n"); return 1; }
        }

        const int *g = H + (size_t)J.g * NC;
        double w0 = now_s();
        for (size_t r = 0; r < nsrc; r++) {
            record_image(src + r * (size_t)RB, g, dst + r * (size_t)RB);
            if (memcmp(dst + r * (size_t)RB, J.row, (size_t)N)) {
                fprintf(stderr, "shard %d: an image does not carry the shard's row 0\n", J.idx);
                return 1;
            }
        }
        qsort(dst, nsrc, (size_t)RB, reccmp);
        for (size_t r = 1; r < nsrc; r++)
            if (!memcmp(dst + (r - 1) * (size_t)RB, dst + r * (size_t)RB, (size_t)RB)) {
                fprintf(stderr, "shard %d: the image has a repeated record\n", J.idx);
                return 1;
            }
        double w = now_s() - w0;

        char sub[1024], path[1200], tmp[1300];
        snprintf(sub, sizeof sub, "%s/%03d", sharddir, J.idx / 1000);
        mkdir(sub, 0777);
        shard_path(path, sizeof path, sharddir, J.idx);
        snprintf(tmp, sizeof tmp, "%s.part", path);
        FILE *f = fopen(tmp, "wb");
        if (!f) { perror(tmp); return 1; }
        if (nsrc && fwrite(dst, (size_t)RB, nsrc, f) != nsrc) { perror(tmp); return 1; }
        if (fflush(f) || fsync(fileno(f)) || fclose(f)) { perror(tmp); return 1; }
        if (rename(tmp, path)) { perror("rename"); return 1; }

        fprintf(out_mf, "{\"idx\":%d,\"row0\":[", J.idx);
        for (int j = 0; j < N; j++) fprintf(out_mf, "%d%s", J.row[j], j + 1 < N ? "," : "");
        fprintf(out_mf, "],\"count\":%zu,\"nodes\":0,\"wall\":%.4f,\"bytes\":%zu,"
                        "\"status\":\"MAPPED\",\"done\":true,\"from\":%d,\"g\":%d}\n",
                nsrc, w, nsrc * (size_t)RB, J.rep, J.g);
        fflush(out_mf);
        made++; records += (long long)nsrc;
    }
    fclose(of);
    if (fclose(out_mf)) { perror(manifest); return 1; }
    printf("{\"mapped_shards\":%lld,\"records\":%lld,\"wall\":%.1f}\n", made, records, now_s() - t0);
    free(done); free(src); free(dst);
    return 0;
}

static void usage(int rc)
{
    fprintf(stderr,
"symmetry -- the plane-fixing subgroup of the cell group, and shard orbits\n"
"\n"
"usage:\n"
"  symmetry <n> group\n"
"        Build the cell group, select the elements fixing the plane x = 0,\n"
"        and check the order against |G| / (3 * 2*floor(n/2)), the index being\n"
"        the number of planes the whole group can send {x = 0} to.\n"
"  symmetry <n> orbits <shardfile> <out.jsonl> <reps.txt>\n"
"        Decompose the shard universe into orbits of that subgroup.  Writes\n"
"        one JSON line per shard -- its orbit, its representative (the least\n"
"        shard index in the orbit) and the index of a subgroup element taking\n"
"        the representative's row 0 to its own -- and a reps.txt in the format\n"
"        `shards list` uses, so it can be fed straight to `enum`.\n"
"        Fails if any image is inadmissible or missing from the shard file.\n"
"  symmetry <n> expand <orbits.jsonl> <sharddir> <manifest.jsonl> [options]\n"
"        For every non-representative shard, read its representative's payload,\n"
"        apply the recorded subgroup element to each record, sort, and write\n"
"        the shard -- byte-identical to a direct enumeration.  Each image is\n"
"        checked to carry the shard's own row 0, and the images are checked to\n"
"        be distinct.  Appends one manifest line per shard, in the same format\n"
"        `enum` uses, so the two manifests concatenate.\n"
"        options:\n"
"          --slice K W    take only shards whose position is K mod W\n"
"          --done FILE    additionally read already-finished shards from FILE\n"
"\n"
"Records are n*n bytes: byte i*n+j is k for the cell (i,j,k).\n");
    exit(rc);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(2);
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) usage(0);
    if (argc < 3) usage(2);
    N = atoi(argv[1]);
    if (N < 2 || N > YMAXN) { fprintf(stderr, "n out of range 2..%d\n", YMAXN); return 2; }
    NC = N * N * N; RB = N * N;
    const char *mode = argv[2];

    if (!strcmp(mode, "group")) return cmd_group();

    if (!strcmp(mode, "orbits")) {
        if (argc < 6) usage(2);
        return cmd_orbits(argv[3], argv[4], argv[5]);
    }

    if (!strcmp(mode, "expand")) {
        if (argc < 6) usage(2);
        int slice_k = 0, slice_w = 1;
        const char *donefile = NULL;
        for (int a = 6; a < argc; a++) {
            if (!strcmp(argv[a], "--slice") && a + 2 < argc) { slice_k = atoi(argv[a+1]); slice_w = atoi(argv[a+2]); a += 2; }
            else if (!strcmp(argv[a], "--done") && a + 1 < argc) { donefile = argv[a+1]; a += 1; }
            else usage(2);
        }
        if (slice_w < 1 || slice_k < 0 || slice_k >= slice_w) { fprintf(stderr, "bad --slice\n"); return 2; }
        return cmd_expand(argv[3], argv[4], argv[5], slice_k, slice_w, donefile);
    }

    usage(2);
    return 2;
}
