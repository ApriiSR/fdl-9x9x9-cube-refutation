#!/usr/bin/env bash
# verify.sh -- reproduce the order-9 result end to end, and check every artefact
# against checksums.txt.
#
#   ./verify.sh full [options]
#       Recompute everything from nothing, including the 1.1 GiB catalogue of
#       all 14 616 576 supports of [9]^3.  That enumeration is the whole cost:
#       about 35 core-hours.  Everything downstream of it takes minutes.
#       Uses no mathematics beyond Lemmas 1-5.
#
#   ./verify.sh symmetric [options]
#       The same catalogue, enumerated one shard per orbit of the subgroup that
#       fixes the plane x = 0, the rest of the shards produced by applying the
#       group element.  157 shards are enumerated instead of 48 912.  This mode
#       additionally relies on Lemma 6, and it checks that reliance: the n = 8
#       census is reproduced the same way, a sample of mapped n = 9 shards is
#       re-enumerated directly and compared, and the catalogue must still have
#       the canonical SHA-256.
#
#   ./verify.sh fast --catalogue FILE [options]
#       Skip the enumeration.  The catalogue's SHA-256 is checked FIRST, before
#       anything is derived from it, and then every downstream step is redone.
#
#   ./verify.sh test [options]
#       The order-8 controls only (tests/test_n8.py).  A couple of minutes.
#
# Options:
#   --workers N   parallel workers                  (default: all cores)
#   --nice N      run the sweep and the searches at nice -n N   (default: 0)
#   --work DIR    scratch and output directory      (default: ./work)
#   --cap SEC     per-shard wall-clock cap in the sweep          (default: 3600)
#   --sample K    mapped shards re-enumerated directly in `symmetric` (default: 320)
#
# A sweep is resumable per shard: rerunning the same command skips the shards
# already recorded as finished, at any worker count.  Once 200 shards (or a
# quarter of them) have finished it prints throughput and an ETA, extrapolated
# from its own run; the shard order is shuffled, so that prefix is an unbiased
# sample of the universe rather than a structurally poor one.
#
# Nothing is downloaded and nothing is sent anywhere.  Every step writes into
# --work; the repository itself is only read.
set -euo pipefail

MODE=${1:-}
shift 2>/dev/null || true
WORKERS=$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4 )
NICE=0
WORK=work
CAP=3600
SAMPLE=320
CATALOGUE=
ETA_AFTER=200

while [ $# -gt 0 ]; do
    case $1 in
        --workers) WORKERS=$2; shift 2 ;;
        --nice) NICE=$2; shift 2 ;;
        --work) WORK=$2; shift 2 ;;
        --cap) CAP=$2; shift 2 ;;
        --sample) SAMPLE=$2; shift 2 ;;
        --catalogue) CATALOGUE=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

case $MODE in
    full|symmetric|fast|test) ;;
    *) sed -n '2,43p' "$0" | sed 's/^# \{0,1\}//'; exit 2 ;;
esac

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
BIN=$HERE/bin
PY=${PYTHON:-python3}
LOG=$WORK/verify.log
NICER=""
if [ "$NICE" -gt 0 ]; then NICER="nice -n $NICE"; fi

say() { printf '\n=== [%5ds] %s\n' "$SECONDS" "$*" | tee -a "$LOG"; }
note() { printf '%s\n' "$*" | tee -a "$LOG"; }
run() { note "\$ $*"; "$@" 2>&1 | tee -a "$LOG"; }

# Throughput and an ETA, from the run's own first shards.
#   $1 manifest directory (m_*.jsonl)   $2 shards to do   $3 already done
eta_watch() {
    local dir=$1 total=$2 base=$3
    local t0=$SECONDS shown=0 last=0 d n el
    while :; do
        sleep 10
        d=$( { cat "$dir"/m_*.jsonl 2>/dev/null || true; } | grep -c '"done":true' || true)
        n=$((d - base)); el=$((SECONDS - t0))
        [ "$n" -gt 0 ] && [ "$el" -gt 0 ] || continue
        if [ "$n" -ge "$ETA_AFTER" ] || [ $((n * 4)) -ge "$total" ]; then
            if [ "$shown" -eq 0 ] || [ $((el - last)) -ge 300 ]; then
                shown=1; last=$el
                awk -v n="$n" -v el="$el" -v tot="$total" -v w="$WORKERS" -v s="$SECONDS" 'BEGIN {
                    rate = n / el; rem = (tot - n) / rate;
                    printf "    [%5ds] %d/%d shards, %.3f s/shard on %d workers, "\
                           "%.2f core-s/shard, ETA %dh%02dm\n",
                           s, n, tot, el / n, w, el * w / n, int(rem / 3600), int(rem % 3600 / 60);
                }' | tee -a "$LOG"
            fi
        fi
    done
}

# Run a shard sweep with `enum`, resumably, with a progress watcher.
#   $1 shard file   $2 shard dir   $3 manifest dir   $4 label
sweep() {
    local shardfile=$1 sharddir=$2 mdir=$3 label=$4
    local total base k pids=() wpid t0 t1
    mkdir -p "$sharddir" "$mdir/logs"
    total=$(grep -c . "$shardfile")
    { cat "$mdir"/m_*.jsonl 2>/dev/null || true; } > "$mdir/done.jsonl"
    base=$(grep -c '"done":true' "$mdir/done.jsonl" || true)
    note "$label: $total shards, $base already finished, $WORKERS workers, cap ${CAP}s"
    t0=$(date +%s)
    k=0
    while [ "$k" -lt "$WORKERS" ]; do
        $NICER $BIN/enum 9 shards "$shardfile" "$sharddir" "$mdir/m_$k.jsonl" \
            --slice "$k" "$WORKERS" --cap "$CAP" --done "$mdir/done.jsonl" \
            > "$mdir/logs/w$k.log" 2>&1 &
        pids+=($!)
        k=$((k + 1))
    done
    eta_watch "$mdir" "$((total - base))" "$base" &
    wpid=$!
    for p in "${pids[@]}"; do wait "$p"; done
    kill "$wpid" 2>/dev/null || true
    wait "$wpid" 2>/dev/null || true
    t1=$(date +%s)
    { cat "$mdir"/m_*.jsonl 2>/dev/null || true; } > "$mdir/enumerated.jsonl"
    note "$label sweep wall: $((t1 - t0)) s on $WORKERS workers"
    note "BUDGET shards (must be 0): $(grep -c '"status":"BUDGET"' "$mdir/enumerated.jsonl" || true)"
    $PY - "$mdir/enumerated.jsonl" <<'EOF' | tee -a "$LOG"
import json, sys
r = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
w = [x['wall'] for x in r if x.get('done')]
w.sort()
print(json.dumps(dict(shards=len(w), supports=sum(x['count'] for x in r if x.get('done')),
                      nodes=sum(x['nodes'] for x in r if x.get('done')),
                      core_seconds=round(sum(w), 1), core_hours=round(sum(w) / 3600, 3),
                      wall_per_shard=dict(min=round(w[0], 4), median=round(w[len(w) // 2], 4),
                                          max=round(w[-1], 4)))))
EOF
}

: > "$LOG"
say "fdlh9 verify.sh $MODE -- $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
note "host: $(uname -srm)  workers: $WORKERS  nice: $NICE  work: $WORK"
note "python: $($PY -c 'import sys,numpy;print(sys.version.split()[0], "numpy", numpy.__version__)')"

say "0. build"
run make -s all
note "cc: $(${CC:-cc} --version 2>&1 | head -1)"

say "1. the order-8 controls (nothing about order 9 is believed until these pass)"
run $PY tests/test_n8.py --workers "$WORKERS"

if [ "$MODE" = test ]; then
    say "done (test mode)"
    exit 0
fi

N9=$WORK/n9
CAT=$WORK/n9_supports.bin

if [ "$MODE" = full ] || [ "$MODE" = symmetric ]; then
    say "2. the shard universe: every admissible row 0, by brute force over all 9!"
    run $BIN/shards 9 count
    $BIN/shards 9 list "$WORK/n9_shards.txt" 2>&1 | tee -a "$LOG"
    # The sweep runs the shards in a fixed shuffled order so that a partial run
    # is an unbiased sample of the universe rather than a structurally poor
    # prefix.  Indices and rows are unchanged by the shuffle.
    $BIN/shards 9 list "$WORK/n9_shards_shuffled.txt" --shuffle 20260906 2>&1 | tee -a "$LOG"
fi

if [ "$MODE" = full ]; then
    say "3. enumerate T(9): 48 912 shards  (this is the expensive step)"
    sweep "$WORK/n9_shards_shuffled.txt" "$N9/shards" "$N9" "T(9)"
    cp "$N9/enumerated.jsonl" "$N9/manifest.jsonl"

elif [ "$MODE" = symmetric ]; then
    say "2b. the subgroup fixing the plane x = 0, and the shard orbits it induces"
    run $BIN/symmetry 9 group
    run $BIN/symmetry 9 orbits "$WORK/n9_shards.txt" "$WORK/n9_shard_orbits.jsonl" \
        "$WORK/n9_reps.txt"

    say "2c. control (a): the same construction at n = 8 must give the census"
    N8=$WORK/n8sym
    rm -rf "$N8"; mkdir -p "$N8/shards"
    run $BIN/shards 8 list "$WORK/n8_shards.txt"
    run $BIN/symmetry 8 orbits "$WORK/n8_shards.txt" "$WORK/n8_shard_orbits.jsonl" \
        "$WORK/n8_reps.txt"
    run $BIN/enum 8 shards "$WORK/n8_reps.txt" "$N8/shards" "$N8/m_reps.jsonl"
    run $BIN/symmetry 8 expand "$WORK/n8_shard_orbits.jsonl" "$N8/shards" "$N8/m_map.jsonl"
    cat "$N8/m_reps.jsonl" "$N8/m_map.jsonl" > "$N8/manifest.jsonl"
    run $PY src/catalogue.py audit 8 "$N8/manifest.jsonl" "$N8/shards"
    run $PY src/catalogue.py pack 8 "$N8/manifest.jsonl" "$N8/shards" "$WORK/n8_supports_sym.bin"
    run $PY src/check.py setcmp 8 "$WORK/n8_supports_sym.bin" data/n8_supports.bin

    say "3. enumerate one shard per orbit  (157 shards, not 48 912)"
    sweep "$WORK/n9_reps.txt" "$N9/shards" "$N9" "T(9) representatives"

    say "3b. map every other shard onto its orbit representative (Lemma 6)"
    t0=$(date +%s); pids=(); k=0
    while [ "$k" -lt "$WORKERS" ]; do
        $NICER $BIN/symmetry 9 expand "$WORK/n9_shard_orbits.jsonl" "$N9/shards" \
            "$N9/x_$k.jsonl" --slice "$k" "$WORKERS" --done "$N9/x_all.jsonl" \
            > "$N9/logs/x$k.log" 2>&1 &
        pids+=($!); k=$((k + 1))
    done
    for p in "${pids[@]}"; do wait "$p"; done
    t1=$(date +%s)
    cat "$N9"/logs/x*.log | tee -a "$LOG"
    cat "$N9"/x_*.jsonl > "$N9/x_all.jsonl"
    note "mapping wall: $((t1 - t0)) s on $WORKERS workers"
    cat "$N9/enumerated.jsonl" "$N9/x_all.jsonl" > "$N9/manifest.jsonl"

    say "3c. control (b): $SAMPLE mapped shards, re-enumerated directly and compared"
    run $PY src/catalogue.py sample 9 "$WORK/n9_shard_orbits.jsonl" "$SAMPLE" "$WORK/n9_sample.txt"
    rm -rf "$WORK/n9_sample_shards"; mkdir -p "$WORK/n9_sample_shards"
    t0=$(date +%s); pids=(); k=0
    while [ "$k" -lt "$WORKERS" ]; do
        $NICER $BIN/enum 9 shards "$WORK/n9_sample.txt" "$WORK/n9_sample_shards" \
            "$WORK/n9_sample_$k.jsonl" --slice "$k" "$WORKERS" --cap "$CAP" \
            > "$N9/logs/s$k.log" 2>&1 &
        pids+=($!); k=$((k + 1))
    done
    for p in "${pids[@]}"; do wait "$p"; done
    t1=$(date +%s)
    note "sample re-enumeration wall: $((t1 - t0)) s on $WORKERS workers"
    run $PY src/catalogue.py setcmpshards 9 "$WORK/n9_sample.txt" \
        "$WORK/n9_sample_shards" "$N9/shards"
fi

if [ "$MODE" = full ] || [ "$MODE" = symmetric ]; then
    say "4. audit the sweep, then assemble the catalogue"
    run $PY src/catalogue.py audit 9 "$N9/manifest.jsonl" "$N9/shards"
    run $PY src/catalogue.py pack 9 "$N9/manifest.jsonl" "$N9/shards" "$CAT"
    say "4b. the catalogue must have the canonical SHA-256 before anything uses it"
    want=$(awk "/^#/ {next} \$2 == \"n9_supports.bin\" {print \$1}" checksums.txt)
    got=$($PY src/check.py sha256 "$CAT" | awk '{print $1}')
    note "expected $want"
    note "actual   $got"
    [ "$want" = "$got" ] || { echo "CATALOGUE CHECKSUM MISMATCH -- refusing to continue" >&2; exit 1; }
    note "catalogue checksum OK"
else
    [ -n "$CATALOGUE" ] || { echo "fast mode needs --catalogue FILE" >&2; exit 2; }
    say "2-4. use the supplied catalogue, checking its SHA-256 before anything reads it"
    want=$(awk "/^#/ {next} \$2 == \"n9_supports.bin\" {print \$1}" checksums.txt)
    got=$($PY src/check.py sha256 "$CATALOGUE" | awk '{print $1}')
    note "expected $want"
    note "actual   $got"
    [ "$want" = "$got" ] || { echo "CATALOGUE CHECKSUM MISMATCH -- refusing to continue" >&2; exit 1; }
    note "catalogue checksum OK"
    CAT=$CATALOGUE
fi

say "5. every record is a support, checked against the definition in numpy"
k=0
while [ "$k" -lt "$WORKERS" ]; do
    $PY src/check.py verify 9 "$CAT" --part "$k" --parts "$WORKERS" > "$WORK/verify_$k.json" &
    k=$((k + 1))
done
wait
cat "$WORK"/verify_*.json | tee -a "$LOG"
$PY - "$WORK" "$WORKERS" <<'EOF' | tee -a "$LOG"
import json, sys, glob, os
tot = bad = 0
for p in sorted(glob.glob(os.path.join(sys.argv[1], 'verify_*.json'))):
    r = json.load(open(p))
    tot += r['checked']
    bad += r['records_without_n2_distinct_cells'] + r['records_missing_a_line_exactly_once']
print(json.dumps(dict(records_checked=tot, failures=bad)))
sys.exit(1 if bad else 0)
EOF

say "6. the catalogue is closed under the order-9216 cell group"
run $BIN/orbits 9 group
run $BIN/orbits 9 closure "$CAT"

say "7. the roots: supports through the centre cell, and their orbits"
run $BIN/orbits 9 centre "$CAT" "$WORK/n9_roots.bin"
run $BIN/orbits 9 classify "$WORK/n9_roots.bin" "$WORK/n9_orbit_sizes.jsonl" "$WORK/n9_orbit_reps.bin"
run $PY src/catalogue.py orbitstats "$WORK/n9_orbit_sizes.jsonl"
run $PY src/check.py verify 9 "$WORK/n9_orbit_reps.bin"

say "8. companion pools, and an independent re-derivation of every one of them"
run $BIN/pools 9 build "$CAT" "$WORK/n9_orbit_reps.bin" "$WORK/pools" "$WORK/n9_pool_sizes.jsonl"
PW=$WORKERS; if [ "$PW" -gt 6 ]; then PW=6; fi   # each worker holds the catalogue's masks
k=0
while [ "$k" -lt "$PW" ]; do
    $PY src/check.py pools 9 "$CAT" "$WORK/n9_orbit_reps.bin" "$WORK/pools" \
        --part "$k" --parts "$PW" > "$WORK/poolchk_$k.json" &
    k=$((k + 1))
done
wait
tail -n1 -q "$WORK"/poolchk_*.json | tee -a "$LOG"
grep -h '"disagreements"' "$WORK"/poolchk_*.json | grep -qv '"disagreements": 0' \
    && { echo "POOL RE-DERIVATION DISAGREES" >&2; exit 1; } || true

say "9. exhaust every root case, and independently compute its packing ceiling"
t0=$(date +%s)
k=0
while [ "$k" -lt "$WORKERS" ]; do
    $NICER $BIN/pack 9 roots "$WORK/n9_orbit_reps.bin" "$WORK/pools" \
        "$WORK/res_$k.jsonl" --slice "$k" "$WORKERS" --cap 600 &
    k=$((k + 1))
done
wait
t1=$(date +%s)
cat "$WORK"/res_*.jsonl > "$WORK/n9_results_raw.jsonl"
note "root search wall: $((t1 - t0)) s on $WORKERS workers"
run $PY src/catalogue.py ledger "$WORK/n9_results_raw.jsonl" "$WORK/n9_root_results.jsonl"

say "10. the exceptional orbit, written out as an explicit packing and re-checked"
J=$($PY - "$WORK/n9_root_results.jsonl" <<'EOF'
import json, sys
best = max((json.loads(l) for l in open(sys.argv[1])), key=lambda r: r['max_companion_clique'])
print(best['j'])
EOF
)
note "the orbit with the largest packing is j = $J"
run $BIN/pack 9 witness "$WORK/n9_orbit_reps.bin" "$WORK/pools" "$J" "$WORK/exceptional_packing.json"
run $PY src/check.py witness 9 "$WORK/exceptional_packing.json"

say "11. checksums"
PROD=$WORK/checksums.produced
: > "$PROD"
$PY src/check.py sha256 "$CAT" | sed "s| .*|  n9_supports.bin|" >> "$PROD"
( cd "$WORK" && $PY "$HERE/src/check.py" sha256 \
    n9_orbit_reps.bin n9_orbit_sizes.jsonl n9_pool_sizes.jsonl \
    n9_root_results.jsonl exceptional_packing.json ) >> "$PROD"
$PY src/check.py sha256 data/n8_supports.bin | sed "s| .*|  n8_supports.bin|" >> "$PROD"
sort -k2 "$PROD" -o "$PROD"
cat "$PROD" | tee -a "$LOG"

say "12. compare with checksums.txt"
$PY - checksums.txt "$PROD" <<'EOF' | tee -a "$LOG"
import sys
def read(p):
    d = {}
    for line in open(p):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        h, name = line.split()
        d[name.split('/')[-1]] = h
    return d
want, got = read(sys.argv[1]), read(sys.argv[2])
bad = 0
for name in sorted(want):
    if name not in got:
        print(f'MISSING   {name}'); bad += 1
    elif want[name] != got[name]:
        print(f'MISMATCH  {name}\n  expected {want[name]}\n  actual   {got[name]}'); bad += 1
    else:
        print(f'OK        {name}  {got[name]}')
print(f'\n{len(want) - bad} of {len(want)} artefacts match')
sys.exit(1 if bad else 0)
EOF

say "done -- $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
