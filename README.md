# There is no fully diagonalised Latin cube of order 9

This repository reproduces, from nothing, a negative answer to a question of
Walter Taylor.

Call an array `A : [n]^d -> [n]` (with `[n] = {0, 1, ..., n-1}`) a **fully
diagonalised Latin hypercube**, or FDLH, if every *main line* of the cube
receives all `n` symbols.  A main line is obtained by letting a parameter `t`
run over `[n]` and setting each coordinate either to a constant, or to `t`, or
to `n-1-t`, with at least one coordinate non-constant.  So the rows, the
columns, the pillars, all the plane diagonals and all the space diagonals must
each be rainbow.  Taylor introduced these in 1972 under the name *completely
Latin*; the `2^d` corner cells are pairwise joined by main lines, which forces
`n = 1` or `n >= 2^d`, and Taylor conjectured that this is the only
obstruction.  He asked in particular (his Problem 2) for the case `d = 3`,
`n = 9`.

**The answer is no.  There is no FDLH of order 9 in dimension 3.**  Taylor's
conjecture that `M(d) = 2^d` (his Problem 3) is therefore false at `d = 3`, and
the first undecided order in dimension 3 is now 10.

The proof is a finite computation.  This repository separates it into the part
a person can check on paper — five short lemmas below, which is the whole of
the mathematics — and the part a machine has to do, which is stated as a table
of programs, inputs, expected outputs and SHA-256 checksums, so that a reader
can rerun any line of it independently.

```
make                      # five C programs, no libraries
./verify.sh test          # the order-8 controls          (~1 minute)
./verify.sh full          # everything, from nothing    (2 h 41 m on 14 cores)
./verify.sh fast --catalogue n9_supports.bin   # everything downstream of the
                                               # catalogue   (8-12 minutes)
```

---

## Part I. The mathematics

Throughout, `n = 9`, `r(t) = n-1-t`, and cells of `[n]^3` are written `(i,j,k)`
and indexed by `(i*n + j)*n + k`.

**Main lines.**  Choose for each of the three coordinates one of: a constant
`c in [n]`, the pattern `t`, or the pattern `r(t)`; not all three constant.  As
`t` runs over `[n]` this traces out `n` cells.  Replacing `t` by `r(t)`
throughout gives the same set, so each line is counted once by requiring the
first non-constant coordinate to use the pattern `t`.  The number of main lines
of `[n]^d` is

```
sum_{s=1}^{d} C(d,s) * 2^{s-1} * n^{d-s},
```

which for `d = 3` is `3n^2 + 6n + 4`: **301** at `n = 9` (243 axis lines, 54
face diagonals, 4 space diagonals), 244 at `n = 8`.

**Support.**  A set `T` of cells meeting every main line exactly once.  Write
`T(n)` for the set of supports of `[n]^3`.

### Lemma 1 (supports and partitions)

> `A : [n]^3 -> [n]` is an FDLH if and only if each of its `n` level sets
> `A^{-1}(v)` is a support.  Conversely `n` pairwise disjoint supports
> partition `[n]^3`, and labelling them by any bijection to `[n]` gives an
> FDLH.  Hence an FDLH of order `n` exists **iff** `T(n)` contains `n` pairwise
> disjoint members.

*Proof.*  A main line has `n` cells, so it receives all `n` symbols iff each
symbol occurs on it exactly once, i.e. iff `|A^{-1}(v) ∩ L| = 1` for every
symbol `v` and every main line `L`.  That is exactly the statement that every
level set is a support, and the level sets partition the cube.

For the converse, first note that any support has exactly `n^2` cells: the `n^2`
lines `{(i,j,*)}` are main lines (only the third coordinate varies), they are
pairwise disjoint, they cover the cube, and `T` meets each once.  So `n`
pairwise disjoint supports occupy `n * n^2 = n^3` cells and therefore partition
`[n]^3`.  Setting `A(c) = v` for the unique `v` with `c ∈ T_v` gives an array
whose level sets are the `T_v`, which by the first paragraph is an FDLH. ∎

### Lemma 2 (Latin-square form, and what row 0 must be)

> Let `T` be a support of `[n]^3`.  Then
>
> (a) `T = {(i, j, L(i,j)) : i, j ∈ [n]}` for a unique array `L : [n]^2 -> [n]`;
>
> (b) `L` is a Latin square;
>
> (c) the permutation `p = L(0, ·)` has **exactly one fixed point** and
> **exactly one reflected point** (exactly one `t` with `p(t) = t`, and exactly
> one `t` with `p(t) = r(t)`).

*Proof.*  (a)  For fixed `(i,j)` the set `{(i,j,k) : k ∈ [n]}` is a main line
(third coordinate `t`, first two constant).  `T` meets it exactly once, so
there is exactly one `k` with `(i,j,k) ∈ T`; call it `L(i,j)`.

(b)  For fixed `i` and `k`, `{(i,j,k) : j ∈ [n]}` is a main line, met exactly
once, so for each `i` and each `k` there is exactly one `j` with `L(i,j) = k`:
every row of `L` is a permutation.  Symmetrically `{(i,j,k) : i ∈ [n]}` gives
that every column is a permutation.  So `L` is a Latin square.

(c)  The plane `i = 0` contains two main lines with the first coordinate held
constant at `0` and the other two varying: `{(0,t,t)}` and `{(0,t,r(t))}`.  `T`
meets the first exactly once, i.e. there is exactly one `t` with
`L(0,t) = t`; and the second exactly once, i.e. exactly one `t` with
`L(0,t) = r(t)`.  By (b), `p = L(0,·)` is a permutation, so it has exactly one
fixed and exactly one reflected point. ∎

Call such a permutation **admissible**.  Lemma 2(c) is what makes the
enumeration shardable: *every* support has an admissible row 0, so enumerating
the supports with each admissible row 0 in turn misses nothing.  Only that
direction of the lemma is load-bearing — the direction in which an error could
hide a support.  The number of admissible permutations is
[OEIS A007016](https://oeis.org/A007016): `8, 20, 96, 656, 5568, 48912` for
`n = 4..9`.  The package computes it twice, by brute force over all `9! =
362 880` permutations (`shards`) and by inclusion–exclusion in closed form
(`check.py a007016`), and both agree with OEIS.

The enumerator itself does **not** use Lemma 2.  It solves the exact-cover
problem

```
items   = the 301 main lines
options = the 729 cells,  a cell covering exactly the lines through it
```

by Knuth's Algorithm X with dancing links.  The Latin structure comes out of
the cover; it is not put in.  Lemma 2 is used only to justify the sharding, and
the enumerator re-checks each shard's row for admissibility before using it.

### Lemma 3 (every cube has a root)

> `n = 9` is odd, so `[9]^3` has a centre cell `c0 = (4,4,4)`.  In any FDLH of
> order 9 exactly one level set contains `c0`.

*Proof.*  The level sets partition the cube. ∎

Call a support containing `c0` a **root**.  So: if no root is a member of any
partition of `[9]^3` into nine supports, there is no FDLH of order 9.  (An
aside, not needed for the proof: all four space diagonals pass through `c0`, so
a support containing it discharges all four with one cell while a support
avoiding it needs four distinct cells.  That is why roots are not rare —
`11 821 056` of the `14 616 576` supports, 80.9 %, contain the centre.)

### Lemma 4 (orbit reduction)

> Let `rev(t) = r(t) = n-1-t` and let `C(rev) <= S_n` be its centraliser.  For
> `pi ∈ S_3`, `eps ∈ {id, rev}^3` and `tau ∈ C(rev)` define a map on cells by
>
> ```
> g(x)_i = tau( eps_i( x_{pi(i)} ) ).
> ```
>
> These form a group `G` which (i) permutes the main lines, hence maps supports
> to supports; (ii) fixes the centre cell; and (iii) maps partitions to
> partitions.  Consequently a root `T` lies in a partition iff `gT` does, and it
> suffices to test one root per `G`-orbit.  At `n = 8` and `n = 9`,
> `|C(rev)| = 2^{floor(n/2)} * floor(n/2)! = 384` and `|G| = 6 * 8 * 384 / 2 =
> **9216**`.

*Proof.*  (i)  Describe a main line by its triple of patterns (constant `c`,
`t`, or `r(t)`) per coordinate.

*Coordinate permutations.*  `pi` permutes which coordinate carries which
pattern.  The image is again a triple of patterns with at least one
non-constant, so a main line.

*Reversals.*  `eps_i = rev` sends a constant `c` in position `i` to the
constant `r(c)`, the pattern `t` to `r(t)`, and `r(t)` to `t`.  Again a main
line.

*Symbol permutations.*  Applying the same `tau` to all three coordinates sends
a constant `c` to the constant `tau(c)`; and it sends the pattern `t` to
`tau(t)`, which after the substitution `s = tau(t)` (a bijection of `[n]`, so a
legitimate reparametrisation) is again the pattern `s`.  The pattern `r(t)`
goes to `tau(r(t)) = tau(rev(t))`, and this equals `rev(tau(t)) = r(s)`
**exactly because `tau` commutes with `rev`** — which is the whole reason
`C(rev)` and not `S_n` appears.  So the image is again a main line, with the
same pattern triple up to the above substitutions.

A bijection of cells that permutes the main lines carries a set meeting every
line once to a set meeting every line once, so `G` maps `T(n)` to `T(n)`.

(ii)  `rev` has the unique fixed point `4`, and `tau rev = rev tau` forces
`tau` to permute the fixed points of `rev`, so `tau(4) = 4`; also
`eps_i(4) = 4`; and `pi` merely permutes the equal coordinates of `(4,4,4)`.
So `g(c0) = c0` for every `g ∈ G`.

(iii)  `g` is a bijection of the cell set, so `T_1, ..., T_9` are pairwise
disjoint supports covering the cube iff `gT_1, ..., gT_9` are.

Finally the order.  The parametrisation `(pi, eps, tau) -> g` is 2-to-1:
`eps = (rev, rev, rev)` with `tau = id` names the same cell map as `eps = id`
with `tau = rev` (both send `x_i` to `r(x_i)`, and `rev ∈ C(rev)`), and pairing
each triple with the one obtained by flipping all three `eps` and composing
`tau` with `rev` is a fixed-point-free involution on the `6 * 8 * |C(rev)|`
triples that leaves `g` unchanged.  Hence `|G| = 6 * 8 * |C(rev)| / 2`.  A `tau` commuting with
`rev` must permute the `floor(n/2)` pairs `{t, r(t)}` as blocks, possibly
flipping each, and must fix the middle point when `n` is odd, so
`|C(rev)| = 2^{floor(n/2)} * floor(n/2)!` — `384` at `n = 8, 9`.  The program
`orbits` builds all `6 * 8 * 384 = 18 432` maps, deduplicates them, and checks
that `9216` remain and that all of them fix `c0`. ∎

### Lemma 5 (the clique bound: what actually kills order 9)

> For a support `T` let its **companion pool** be
> `P(T) = { S ∈ T(n) : S ∩ T = ∅ }`, and let `Gamma(T)` be the graph on `P(T)`
> joining two companions when they are disjoint.  If `T` belongs to a partition
> of `[n]^3` into `n` supports, then the other `n-1` members lie in `P(T)` and
> are pairwise disjoint — a clique of size `n-1` in `Gamma(T)`.  Hence
>
> ```
> omega(Gamma(T)) < n - 1   ==>   T belongs to no partition.
> ```
>
> In particular, at `n = 9`, `omega(Gamma(T)) <= 4` suffices, and a
> *triangle-free* `Gamma(T)` (`omega <= 2`) suffices very comfortably.

*Proof.*  Immediate: the other `n-1` members of the partition are supports
disjoint from `T`, hence in `P(T)` by completeness of the catalogue, and
pairwise disjoint, hence a clique. ∎

The point of Lemma 5 is that it is not a search.  A depth-first exact cover
that finds nothing can fail to find something for a bad reason; a clique number
is a property of a finite graph, computed by a different program, and it comes
out not merely below 8 but at 2 or 4.

### Theorem

> There is no fully diagonalised Latin cube of order 9.  Equivalently, Taylor's
> `P(9,3)` is false: his Problem 2 has a negative answer, and his conjecture
> that `n >= 2^d` is the only restriction (Problem 3, `M(d) = 2^d`) fails at
> `d = 3`.

*Proof.*  Suppose `A` is an FDLH of order 9.  By Lemma 1 its level sets are
nine pairwise disjoint supports.  By Lemma 3 one of them, `T`, contains the
centre `c0`.  By Lemma 4 we may replace the whole cube by its image under any
`g ∈ G`, so we may assume `T` is the chosen representative of its `G`-orbit.
The computation (Part II) establishes:

* `T(9)` has exactly `14 616 576` members, enumerated exhaustively;
* `11 821 056` of them contain `c0`, falling into exactly `2 049` `G`-orbits
  whose sizes sum to `11 821 056`;
* for each of the `2 049` orbit representatives, the companion pool is complete
  (it is a filter over the complete catalogue, re-derived independently), and
  `omega(Gamma(T)) <= 4` — in fact `2` for `2 048` of them, `4` for one.

By Lemma 5, `T` belongs to no partition.  Contradiction. ∎

Independently of Lemma 5, the same computation runs a depth-first exact cover
over each of the `2 049` pools: all `2 049` are exhausted, with `0` covers
found, and no branch survives past a root plus **one** companion.  Either
argument alone is sufficient; they are produced by different code on the same
data.

### An aside: the best packing there is

Nothing above needs a witness, but one is instructive.  The single exceptional
orbit — the one whose pool graph has clique number 4 rather than 2 — yields a
set of **five** pairwise disjoint supports of `[9]^3`, written out explicitly in
`data/exceptional_packing.json` and re-verified from the definition by
`check.py witness`.  A cube would need nine.  Any future proof of "the packing
dies early at order 9" that does not survive this example is wrong.

---

## Part II. The computation

### The catalogue is a canonical object

The expensive artefact is `n9_supports.bin`: all `14 616 576` supports of
`[9]^3`, `81` bytes each (`byte i*9+j` is `k` for the cell `(i,j,k)`), `1 183
942 656` bytes in all.  It is written in a **canonical order**: shards in
increasing shard index — which is the lexicographic rank of the shard's row 0
among the admissible permutations — and records sorted lexicographically inside
each shard.  Since a record begins with its own row 0, this is simply *all
supports in lexicographic order*.  Its SHA-256 therefore depends on the set
alone, not on which program found it or in what order, which is why comparing
that one hash is a meaningful check between independent implementations.

The same set was enumerated previously by an unrelated program: a constraint
solver applied to the Latin-square model of Lemma 2, sharded the same way.  Put
into the same canonical order, that earlier enumeration has the identical
SHA-256, `f9dd54e4...`.  The two enumerations share no solver, no model and no
code; the exact-cover route in this repository is about three times cheaper.

### What each step does

Times below were measured by the runs recorded in this README (see *Measured
runtimes*), not estimated.  Every artefact's SHA-256 is in `checksums.txt` and
is checked by `verify.sh`.

| # | step | program | input | expected output |
|---|------|---------|-------|-----------------|
| 1 | order-8 controls | `tests/test_n8.py` | — | 23 checks pass (below) |
| 2 | shard universe | `shards 9 list` | — | 48 912 admissible rows = A007016(9) |
| 3 | enumerate `T(9)` | `enum 9 shards` | shard list | 48 912 shards EXHAUSTED, 14 616 576 supports |
| 4 | assemble, audit | `catalogue.py` | shard payloads | `n9_supports.bin`, no duplicates, every record carrying its own shard row |
| 5 | definition check | `check.py verify 9` | catalogue | 14 616 576 records, 0 failures against the 301 main lines |
| 6 | group closure | `orbits 9 closure` | catalogue | 87 699 456 images checked, 0 missing |
| 7 | roots | `orbits 9 centre` | catalogue | 11 821 056 supports through `(4,4,4)` |
| 8 | orbits | `orbits 9 classify` | roots | 2 049 orbits, sizes summing to 11 821 056 |
| 9 | companion pools | `pools 9 build` | catalogue, 2 049 reps | pool sizes 764 / 980 / 988.3 / 1 552 (min/median/mean/max) |
| 10 | pools re-derived | `check.py pools 9` | catalogue, pools | 0 disagreements |
| 11 | exhaust + cliques | `pack 9 roots` | 2 049 reps, pools | 2 049 EXHAUSTED, 0 covers, max depth 2, 183 625 nodes; max packing 3 (2 048 orbits) and 5 (1 orbit) |
| 12 | witness | `pack 9 witness` | exceptional orbit | 5 pairwise disjoint supports, re-verified |
| 13 | checksums | `check.py sha256` | all artefacts | match `checksums.txt` |

### The artefacts and their checksums

| artefact | size | SHA-256 |
|---|---:|---|
| `n9_supports.bin` (not in the repository) | 1 183 942 656 B | `f9dd54e401c69327d9383567e0050ca65ffd9bfeca7d8b7e8609059508fe3eb8` |
| `data/n8_supports.bin` | 835 584 B | `f780d243569b1e54d1d7482ffad7c08aa70e2c7caf9849d9c553fac640f15abc` |
| `data/n9_orbit_reps.bin` | 165 969 B | `36fcd900c16adb496c5132b8205826e1ad99a90075058419ac379881008bdb80` |
| `data/n9_orbit_sizes.jsonl` | 58 300 B | `f5a98d97d93fa881018058822afe2088a78725791081b93a0ee3ac20b7f99b43` |
| `data/n9_pool_sizes.jsonl` | 44 769 B | `1528986d3c6391065114d9eceb5c8fd5a00fbe4b16d5351364d75bd63442d39d` |
| `data/n9_root_results.jsonl` | 362 770 B | `55f3730690d06a0fc2ae37f33e2e7b1df9a49cdb446d5e9a92f5fc2e38cce386` |
| `data/exceptional_packing.json` | 940 B | `1c4fd1c21289ac11f1d569d929bead3a13b155e17cd07b4a0a97bbd30f028549` |

`n9_roots.bin` (957 505 536 B, the 11 821 056 supports through the centre) is a
one-second filter over the catalogue and is not checksummed separately.
`n9_root_results.jsonl` carries no timings, so that it is reproducible byte for
byte; the raw per-query output with wall times is left beside it as
`n9_results_raw.jsonl`.

### The order-8 controls

Order 8 is where a cube is known to exist and its census is known
independently, so the whole pipeline is exercised at an order where the answers
were not produced by this code.  `./verify.sh test` runs 23 checks; the
substantive ones are

* `T(8)` enumerated from nothing is **13 056** supports, set-equal *and*
  byte-identical to the shipped `data/n8_supports.bin`;
* all 13 056 pass the definition-level check against the 244 main lines;
* `T(8)` is closed under the order-9216 cell group and falls into **6** orbits,
  of sizes 768, 768, 2 304, 2 304, 2 304, 4 608;
* the exact cover of `[8]^3` by supports has exactly **198 624** solutions, with
  node counts by depth `1, 1 632, 26 016, 60 672, 69 513, 147 292, 154 660,
  198 624, 198 624`.  Relabelling the
  eight classes of a partition by two different permutations gives two
  different arrays (a permutation moving `v` to `w != v` changes the array,
  since the classes `T_v` and `T_w` are different sets), so this also counts
  the labelled cubes: there are `8! * 198 624 = 8 008 519 680` FDLHs of
  order 8;
* the **root search agrees with the census on every one of the 13 056 supports**
  — for each support, the number of covers containing it, computed once by the
  census and once by the same root machinery that is run at order 9.  Eight
  spot checks would make the point; running all 13 056 costs two seconds, so
  all of them are run, and eight sampled ones are also printed individually;
* the packing ceiling and the search agree at order 8 too: a support lies in a
  cover **iff** its pool graph has a 7-clique;
* all four wall-clock caps are tested to fire, with the shard-level cap shown to
  leave no payload behind, not to mark the shard done, and to have it redone on
  the next pass.  A stopping rule that never fires looks exactly like an
  exhaustion, which is the failure mode a negative result most needs excluded.

### Measured runtimes

Two machines, both arm64.  **A** is a 20-core Linux box (Cortex-X925 +
Cortex-A725, 121 GB), running `full` with 14 workers at `nice -n 10` while
other jobs of its owner's were also running.  **B** is a 16-core macOS laptop
(64 GB), running `test` and `fast` with 14 workers.

`verify.sh full` on **A**, start to finish: **2 h 41 m 28 s**.

| step | measurement |
|---|---|
| order-8 controls (23 checks), run on their own | 12 s on A, 21 s on B |
| enumerate `T(9)`, 48 912 shards | **9 271 s wall** on 14 workers; **124 188 core-seconds = 34.5 core-hours**; 7.83 x 10^11 search nodes; wall per shard min 0.000 s, median 2.486 s, max 11.279 s; 0 shards hit the 3 600-second cap |
| assemble the catalogue and audit it | part of the 417 s that everything other than the sweep took on A |
| check all 14 616 576 records against the definition | 18 s on B (14 numpy workers) |
| closure under the cell group | 50.2 s on A, 43.5 s on B |
| extract the 11 821 056 roots | 1.4 s |
| classify them into 2 049 orbits | 48.9 s on A, 39.9 s on B |
| build 2 049 companion pools | 124.9 s on A, 91.8 s on B |
| re-derive all 2 049 pools in numpy | 320 s on B (6 workers) |
| exhaust 2 049 root cases and compute their cliques | **1 s** wall on 14 workers: 1.70 CPU-seconds of exact cover and 4.5 CPU-seconds of clique computation for all 2 049 |
| extract and re-check the 5-packing | under a second |

`verify.sh fast` on **B**, start to finish: **448 s** and **740 s** on two
runs — the spread is almost all in the numpy re-derivation of every companion
pool (320 s against 493 s), which is bound by how much of the 1.1 GiB
catalogue is in page cache.  Of the 448 s, 21 s is the order-8 suite and 320 s
is that re-derivation: everything the order-9 argument actually needs, given
the catalogue, is about two minutes, and the rest is checking.

The shape of the cost is worth stating plainly: **the enumeration is the whole
expense and the refutation is free.**  Exhausting all 2 049 root cases costs
1.70 CPU-seconds; the independent clique computation costs another 4.5.
Building the object they run over costs 34.5 core-hours.

### Reproducing it

```
make                                     # bin/{enum,shards,orbits,pools,pack}
./verify.sh full  --workers 14           # everything
```

Requirements: a C99 compiler, GNU make, `bash`, and Python 3 with **numpy**
(the only Python dependency).  If numpy lives in a virtualenv, set
`PYTHON=/path/to/python`.  Disk: about 3.5 GB under `--work`.  Memory: the
largest single process holds the catalogue plus one bitmask per record, about
3 GB; the parallel pool re-derivation is capped at 6 workers for that reason.

`./verify.sh fast --catalogue FILE` skips step 3 and starts from a catalogue you
already have.  It checks that file's SHA-256 **before** anything reads it and
refuses to continue on a mismatch, so a wrong or truncated catalogue cannot
silently propagate into the downstream results.

Every binary takes `--help`.  Nothing in this repository downloads anything or
sends anything anywhere; all output goes under `--work` (default `./work`).

---

## What a referee still has to take on trust

The point of Part I is that the reduction from "is there a cube of order 9" to
"here are 2 049 finite graphs, look at their clique numbers" is checkable by
hand.  What remains is:

1. **That the programs implement the definitions.**  The main lines are built
   twice — in C from the pattern description, and independently in numpy — and
   the two counts and incidence structures are compared; every record ever
   produced is checked against the numpy version of the definition.  But the
   definition itself is written down twice by the same author, and a reader who
   disagrees with the definition of a main line will disagree with everything
   downstream.  Read `src/lines.h` and `main_lines()` in `src/check.py`; they
   are about twenty lines each.

2. **That the exhaustive enumeration is exhaustive.**  Two things support this
   beyond the code: at order 8 the identical pipeline reproduces a census that
   is known independently; and the enumerated `T(9)` is *closed under the
   order-9216 group* (87 699 456 images checked, none missing), so a missing
   support would have had to be missed together with its entire orbit.  A
   second, unrelated enumeration by a constraint solver on a different model
   agrees exactly (same canonical SHA-256) — but that agreement is history this
   repository does not itself recompute; what it recomputes is the enumeration
   and the checks.

3. **Common-mode error in the model.**  Both the exact-cover enumerator and the
   definition-level checker are built on the same reading of "main line" and of
   the record format.  The order-8 controls exclude a good deal of this: they
   agree with a census, a node profile, an orbit decomposition and an
   automorphism-free count that were established elsewhere and are quoted here
   as targets, not derived.

4. **Floating point plays no role**; there is none in the argument.  Wall-clock
   caps are the only timing-dependent behaviour, and they are configured to
   values the sweep never approaches (the slowest shard takes about 16 seconds
   against a 3 600-second cap), and are separately tested to fire.

5. **The scope.**  This says nothing about order 10 or order 12 in dimension 3,
   and nothing about dimension 4 or above.  Order 10 is now the case that
   decides whether the threshold in dimension 3 is 10 or 11.

---

## Files

```
README.md            this file
Makefile             builds the five C programs into bin/
verify.sh            full / fast / test
checksums.txt        SHA-256 of every data artefact
src/lines.h          the main lines of [n]^3, from the definition
src/enum.c           the enumerator: exact cover by dancing links
src/shards.c         the admissible row-0 permutations (A007016)
src/orbits.c         the order-9216 cell group; closure and orbit classification
src/pools.c          companion pools as a filter over a complete catalogue
src/pack.c           exhaustive exact cover, and the packing ceiling
src/check.py         definition-level checks in numpy, independent of the C
src/catalogue.py     assembly, audit and ledger canonicalisation
tests/test_n8.py     the order-8 controls
data/n8_supports.bin          the 13 056 supports of [8]^3            (835 KB)
data/n9_orbit_reps.bin        the 2 049 root orbit representatives    (166 KB)
data/n9_orbit_sizes.jsonl     their orbit sizes
data/n9_pool_sizes.jsonl      their companion pool sizes
data/n9_root_results.jsonl    per-root: pool, covers, nodes, depth, edges,
                              triangles, clique number
data/exceptional_packing.json the five pairwise disjoint supports
```

The 1.1 GiB catalogue and the 957 MB root file are not in the repository.
Recompute them with `./verify.sh full`, or obtain the catalogue elsewhere and
use `./verify.sh fast`, which checks its SHA-256 first.

## References

* W. Taylor, *On the coloration of cubes*, Discrete Mathematics **2** (1972)
  187–190.  The objects, the `n >= 2^d` corner bound, the conjecture, and
  Problem 2 (`d = 3`, `n = 9`).
* D. E. Knuth, *Dancing links*, in J. Davies, B. Roscoe, J. Woodcock (eds.),
  *Millennial Perspectives in Computer Science*, Palgrave (2000), 187–214.
  Algorithm X and the doubly linked cover matrix used by `src/enum.c`.
* [OEIS A007016](https://oeis.org/A007016) — permutations with exactly one
  fixed point and exactly one reflected point; the shard count.
