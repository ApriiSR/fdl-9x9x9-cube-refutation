/* util.h -- loading, validating and timing, shared by the C programs.
 *
 * Three things live here because all six programs need them and getting any of
 * them subtly wrong would be a silent failure rather than a loud one:
 *
 *   LOADING.  A record file is n*n bytes per record, byte i*n+j being k for the
 *   cell (i,j,k).  Every byte must therefore be below n.  A byte that is not is
 *   an index out of the cell map, the mask allocation and the image buffer of
 *   every program that reads it, so the range is checked once, on the way in,
 *   before anything indexes with it.
 *
 *   TIMING.  A monotonic clock, with a deterministic injected substitute: if
 *   FDLH_CLOCK_STEP is set to a positive number of seconds, every reading
 *   advances a virtual clock by exactly that much.  A wall-clock cap can then be
 *   tested to fire at an exact point rather than "eventually", which is what a
 *   stopping rule needs if the test is to mean anything on a loaded machine.
 *
 *   MANIFESTS.  A manifest is append-only JSON lines.  A kill during a write
 *   can leave a final line without its newline; appending after it would join
 *   two records into one.  fdlh_manifest_append() truncates a torn final line
 *   before returning the handle, so the only record ever lost is the one that
 *   was in flight.  fdlh_fnv64() is the payload digest recorded alongside, so a
 *   resumed run can tell a completed shard from a corrupted one.  It is a
 *   64-bit non-cryptographic checksum: it catches corruption, it is not a
 *   defence against a forged payload.  The catalogue's SHA-256 is that.
 */
#ifndef FDLH_UTIL_H
#define FDLH_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

/* ---- timing -------------------------------------------------------------- */
static inline double fdlh_now_s(void)
{
    static int injected = -1;
    static double step, virt;
    if (injected < 0) {
        const char *e = getenv("FDLH_CLOCK_STEP");
        step = (e && *e) ? atof(e) : 0.0;
        injected = step > 0.0;
        virt = 0.0;
    }
    if (injected) { virt += step; return virt; }
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ---- allocation ---------------------------------------------------------- */
static inline void *fdlh_alloc(size_t nbytes, const char *what)
{
    void *p = malloc(nbytes ? nbytes : 1);
    if (!p) { fprintf(stderr, "out of memory allocating %s (%zu bytes)\n", what, nbytes); exit(1); }
    return p;
}

static inline void *fdlh_calloc(size_t count, size_t size, const char *what)
{
    void *p = calloc(count ? count : 1, size);
    if (!p) { fprintf(stderr, "out of memory allocating %s (%zu x %zu)\n", what, count, size); exit(1); }
    return p;
}

/* ---- digest -------------------------------------------------------------- */
static inline uint64_t fdlh_fnv64(const void *p, size_t nbytes)
{
    const unsigned char *b = (const unsigned char *)p;
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < nbytes; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

/* ---- loading ------------------------------------------------------------- */
/* Every byte of a record file is a coordinate in [n].  Reject anything else
 * before it is used as an index. */
static inline int fdlh_bytes_in_range(const unsigned char *b, long long nbytes, int n,
                               const char *path, int rb)
{
    for (long long i = 0; i < nbytes; i++)
        if (b[i] >= (unsigned char)n) {
            fprintf(stderr, "%s: record %lld byte %lld is %d, not a coordinate in [%d]\n",
                    path, i / rb, i % rb, b[i], n);
            return 0;
        }
    return 1;
}

/* Read a whole record file, checking the size is a multiple of rb and that
 * every byte is a coordinate in [n].  Exits on any failure. */
static inline unsigned char *fdlh_load_records(const char *path, int rb, int n, long long *count)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    if (fseek(f, 0, SEEK_END)) { perror(path); exit(1); }
    long long sz = ftell(f);
    if (sz < 0) { perror(path); exit(1); }
    if (fseek(f, 0, SEEK_SET)) { perror(path); exit(1); }
    if (sz % rb) { fprintf(stderr, "%s: size %lld is not a multiple of %d\n", path, sz, rb); exit(1); }
    unsigned char *b = (unsigned char *)fdlh_alloc((size_t)sz, path);
    if (sz && fread(b, 1, (size_t)sz, f) != (size_t)sz) { perror(path); exit(1); }
    fclose(f);
    if (!fdlh_bytes_in_range(b, sz, n, path, rb)) exit(1);
    *count = sz / rb;
    return b;
}

/* ---- manifests ----------------------------------------------------------- */
/* Open a JSON-lines manifest for appending, first discarding a final line that
 * has no newline -- the signature of a kill during a write.  Interior lines are
 * never touched: a corrupt line that is NOT the last one is a real
 * inconsistency and the readers reject it rather than repairing it. */
static inline FILE *fdlh_manifest_append(const char *path)
{
    FILE *f = fopen(path, "r+");
    if (f) {
        if (!fseek(f, 0, SEEK_END)) {
            long sz = ftell(f), keep = 0, pos = sz;
            while (pos > 0) {
                if (fseek(f, --pos, SEEK_SET)) { keep = sz; break; }
                if (fgetc(f) == '\n') { keep = pos + 1; break; }
            }
            if (sz > 0 && keep != sz) {
                fprintf(stderr, "%s: discarding %ld bytes of a torn final line\n", path, sz - keep);
                if (ftruncate(fileno(f), (off_t)keep)) { perror(path); fclose(f); return NULL; }
            }
        }
        fclose(f);
    }
    return fopen(path, "a");
}

/* A manifest line is accepted as a completion record only if it has the exact
 * shape the writers produce: it must begin with {"idx":<int>, end with }, and
 * carry "done":true.  A torn or truncated line fails the closing test, and a
 * line from some other producer fails the opening one. */
static inline int fdlh_manifest_done_line(const char *line, int *idx)
{
    const char *p = line;
    size_t len = strlen(p);
    while (len && (p[len - 1] == '\n' || p[len - 1] == '\r' || p[len - 1] == ' ')) len--;
    if (len < 10 || p[len - 1] != '}') return 0;
    if (strncmp(p, "{\"idx\":", 7)) return 0;
    const char *q = p + 7;
    if (*q < '0' || *q > '9') return 0;
    long v = strtol(q, (char **)&q, 10);
    if (*q != ',') return 0;
    if (!strstr(p, "\"done\":true")) return 0;
    if (v < 0 || v > 0x7fffffffL) return 0;
    *idx = (int)v;
    return 1;
}

/* Pull a "key":<integer> out of a manifest line.  Returns 0 if absent. */
static inline int fdlh_manifest_field(const char *line, const char *key, long long *out)
{
    char pat[64];
    snprintf(pat, sizeof pat, "\"%s\":", key);
    const char *q = strstr(line, pat);
    if (!q) return 0;
    q += strlen(pat);
    char *end;
    long long v = strtoll(q, &end, 10);
    if (end == q) return 0;
    *out = v;
    return 1;
}

/* Pull a "key":"<hex>" out of a manifest line. */
static inline int fdlh_manifest_hex(const char *line, const char *key, uint64_t *out)
{
    char pat[64];
    snprintf(pat, sizeof pat, "\"%s\":\"", key);
    const char *q = strstr(line, pat);
    if (!q) return 0;
    q += strlen(pat);
    char *end;
    uint64_t v = strtoull(q, &end, 16);
    if (end == q || *end != '"') return 0;
    *out = v;
    return 1;
}

/* The set of shards a previous run finished, with what it says their payloads
 * should be.  Open addressing on the shard index; the capacity is fixed and
 * generous (a manifest has one line per shard). */
struct fdlh_done {
    int *idx;                 /* -1 = empty slot */
    long long *bytes;
    uint64_t *digest;
    int *has_digest;
    size_t mask;
    long long n;
};

static inline void fdlh_done_init(struct fdlh_done *d, size_t hint)
{
    size_t hs = 1024;
    while (hs < hint * 4) hs <<= 1;
    d->mask = hs - 1;
    d->n = 0;
    d->idx = (int *)fdlh_alloc(sizeof(int) * hs, "done table");
    d->bytes = (long long *)fdlh_alloc(sizeof(long long) * hs, "done table");
    d->digest = (uint64_t *)fdlh_alloc(sizeof(uint64_t) * hs, "done table");
    d->has_digest = (int *)fdlh_alloc(sizeof(int) * hs, "done table");
    for (size_t i = 0; i <= d->mask; i++) d->idx[i] = -1;
}

static inline size_t fdlh_done_slot(const struct fdlh_done *d, int idx)
{
    size_t s = ((uint64_t)idx * 11400714819323198485ULL) & d->mask;
    while (d->idx[s] >= 0 && d->idx[s] != idx) s = (s + 1) & d->mask;
    return s;
}

static inline void fdlh_done_put(struct fdlh_done *d, int idx, long long bytes,
                          uint64_t digest, int has_digest)
{
    size_t s = fdlh_done_slot(d, idx);
    if (d->idx[s] < 0) { d->idx[s] = idx; d->n++; }
    d->bytes[s] = bytes; d->digest[s] = digest; d->has_digest[s] = has_digest;
}

static inline int fdlh_done_get(const struct fdlh_done *d, int idx, long long *bytes,
                         uint64_t *digest, int *has_digest)
{
    size_t s = fdlh_done_slot(d, idx);
    if (d->idx[s] < 0) return 0;
    *bytes = d->bytes[s]; *digest = d->digest[s]; *has_digest = d->has_digest[s];
    return 1;
}

/* Read every completion record out of a manifest.  Lines that are not exactly
 * a completion record are ignored here (the audit rejects them); lines that are
 * one but lack a length or a digest are recorded as unverifiable, and the
 * caller then redoes that shard rather than trusting it. */
static inline void fdlh_done_read(struct fdlh_done *d, const char *path)
{
    FILE *mf = fopen(path, "r");
    if (!mf) return;
    char line[8192];
    while (fgets(line, sizeof line, mf)) {
        int idx;
        if (!fdlh_manifest_done_line(line, &idx)) continue;
        long long bytes = -1;
        uint64_t dig = 0;
        int has = fdlh_manifest_hex(line, "digest", &dig);
        if (!fdlh_manifest_field(line, "bytes", &bytes)) bytes = -1;
        fdlh_done_put(d, idx, bytes, dig, has && bytes >= 0);
    }
    fclose(mf);
}

/* A completed shard may be skipped only if its payload is still there, still
 * the recorded length, and still hashes to the recorded digest. */
static inline int fdlh_payload_ok(const char *path, long long want_bytes, uint64_t want_digest,
                           const char **why)
{
    FILE *f = fopen(path, "rb");
    if (!f) { *why = "payload missing"; return 0; }
    if (fseek(f, 0, SEEK_END)) { fclose(f); *why = "payload unreadable"; return 0; }
    long long sz = ftell(f);
    if (sz != want_bytes) { fclose(f); *why = "payload is the wrong length"; return 0; }
    if (fseek(f, 0, SEEK_SET)) { fclose(f); *why = "payload unreadable"; return 0; }
    uint64_t h = 1469598103934665603ULL;
    unsigned char buf[1 << 16];
    size_t got;
    while ((got = fread(buf, 1, sizeof buf, f)) > 0)
        for (size_t i = 0; i < got; i++) { h ^= buf[i]; h *= 1099511628211ULL; }
    int err = ferror(f);
    fclose(f);
    if (err) { *why = "payload unreadable"; return 0; }
    if (h != want_digest) { *why = "payload does not match its recorded digest"; return 0; }
    return 1;
}

#endif
