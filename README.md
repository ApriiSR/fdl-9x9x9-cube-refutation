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
conjecture that `n >= 2^d` is the only obstruction is therefore false at
`d = 3`.  This computation says nothing about order 10 or order 12, and
establishes no general threshold; see *The scope* at the end.

The proof is a finite computation.  This repository separates it into the part
a person can check on paper — five short lemmas below, which is the whole of
the mathematics — and the part a machine has to do, which is stated as a table
of programs, inputs, expected outputs and SHA-256 checksums, so that a reader
can rerun any line of it independently.  A sixth lemma is stated and proved
afterwards; it proves nothing about order 9 and is used only to make one of the
three verification modes cheap.

```
make                      # six C programs, no libraries
./verify.sh test          # the order-8 controls          (12-97 s)
./verify.sh full          # everything, from nothing    (2 h 41 m on 14 cores)
./verify.sh symmetric     # the same, one shard per symmetry orbit (9-24 min)
./verify.sh fast --catalogue n9_supports.bin   # everything downstream of the
                                               # catalogue     (6-13 minutes)
```

Three modes, because there are three things a reader might want.  `full`
recomputes the catalogue from nothing and uses no mathematics beyond Lemmas 1-5.
`symmetric` recomputes it too, but enumerates only 157 of the 48 912 shards and
obtains the rest by symmetry — three hundred times less search, at the price of
one extra lemma, which it checks three ways.  `fast` takes the catalogue as
given, verifies its SHA-256 before reading it, and redoes everything downstream.

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
enumeration shardable.  Admissibility is a *necessary* condition on row 0, not
a sufficient condition for extension to a support: exhaustive enumeration
therefore requires every admissible row to be considered, including rows whose
shards turn out to be empty.  That necessary direction is all the sharding
argument uses, and no converse is asserted or needed.  The number of admissible
permutations is
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

Finally, that these maps form a group, and how many of them there are.

*Closure.*  If `g = (pi, eps, tau)` and `h = (rho, delta, sigma)`, then
`g ∘ h` is the map with parameters

```
pi'   = rho ∘ pi,
eps'_i = eps_i ∘ delta_{pi(i)},
tau'  = tau ∘ sigma,
```

which is again of the displayed form: `tau ∘ sigma` lies in `C(rev)` because
`C(rev)` is a subgroup, and the reversals can be gathered on the left of the
symbol permutation precisely because `tau` and `sigma` each commute with `rev`.
The identity is `(id, id, id)`, and a finite composition-closed family of
bijections of a finite set contains inverses.  So `G` is a group.

*The count.*  Suppose two parameter triples define the same cell map.  Each
output coordinate depends bijectively on exactly one input coordinate, so the
two coordinate permutations agree; write `pi` for both.  Equating the `i`-th
output coordinates gives `tau(eps_i(x)) = tau'(eps'_i(x))` for all `x`, i.e.
`tau'^{-1} tau = eps'_i eps_i^{-1}` for every `i`.  The right-hand side is
`id` or `rev`, and it is the *same* element for every `i` because the left-hand
side does not depend on `i`.  So either all three reversal bits agree and
`tau = tau'` — the triples are identical — or all three are flipped and
`tau' = tau ∘ rev`.  Every map therefore has exactly two descriptions, and
`|G| = 6 * 8 * |C(rev)| / 2`.

A `tau` commuting with `rev` permutes the `floor(n/2)` pairs `{t, r(t)}` as
blocks and may flip each, and fixes the middle point when `n` is odd, so
`|C(rev)| = 2^{floor(n/2)} * floor(n/2)!` — `384` at `n = 8` and `n = 9`, and
`|G| = 6 * 8 * 384 / 2 = 9216`.  The program `orbits` builds all
`6 * 8 * 384 = 18 432` maps, deduplicates them, and checks that `9216` remain,
that all of them fix `c0`, and that six explicitly named generators generate
all of them. ∎

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

The recorded counts already suffice, without any clique computation at all.  An
8-clique contains `C(8,3) = 56` triangles, whereas every one of the 2 049
computed companion graphs has either `0` or `8` triangles.  None of them can
contain an 8-clique, so by Lemma 5 no root lies in a partition.  A triangle
count is an elementary property of a finite graph; it is what `catalogue.py
validate` checks, and it is the version of the argument a reader can audit
fastest.

`pack` also computes each maximum clique outright, by the recursive branch and
bound `bk()`, and gets 2 for 2 048 of the orbits and 4 for one.  That is a
different algorithm from the exact cover — but the two live in the same program
and share its record loader, its cell masks and its `disjoint()`, so they are
independent as *algorithms*, not as implementations.

### Theorem

> There is no fully diagonalised Latin cube of order 9.  Taylor's Problem 2
> therefore has a negative answer, and his conjecture that `n >= 2^d` is the
> only restriction fails at `d = 3`.

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
argument alone is sufficient.  They are different algorithms run over the same
pools by the same program; `catalogue.py validate` requires both of them, and
the triangle bound, to hold before the run is called a success.

### An aside: a maximum packing containing a root

Nothing above needs a witness, but one is instructive.  The single exceptional
orbit — the one whose pool graph has clique number 4 rather than 2 — yields a
set of **five** pairwise disjoint supports of `[9]^3`, written out explicitly in
`data/exceptional_packing.json` and re-verified from the definition by
`check.py witness`.  A cube would need nine.  Any future proof of "the packing
dies early at order 9" that does not survive this example is wrong.

This is a maximum packing *containing a root*, which is all the computation
bounds.  It is not a claim about incomplete packings that avoid the centre
cell: Lemma 3 applies to partitions, not to arbitrary packings, so nothing here
says five is the global maximum.

### Lemma 6 (the plane-fixing subgroup, and the shard orbits)

**The theorem above does not use this lemma**, and neither do `verify.sh full`
and `verify.sh fast`.  It is used by `verify.sh symmetric` alone — the mode that
enumerates 157 shards instead of 48 912 — and is stated and proved here so that
a reader who wants the cheap mode can see exactly what extra is being asked of
them, and one who does not can run `full` and lose nothing.

Write `P = {x_0 = 0}` for the plane whose contents define a shard: by Lemma 2 a
support's intersection with `P` is `{(0, j, p(j))}` for an admissible
permutation `p`, and the shard `S_p` is the set of supports with that row 0.

> Let `H = { g ∈ G : g(P) = P }` be the setwise stabiliser of `P` in the group
> `G` of Lemma 4.  Then
>
> (a) `H` consists of exactly those `g = g(pi, eps, tau)` with `pi(0) = 0` and
> `tau(eps_0(0)) = 0`, and
>
> ```
> |H| = |G| / (3 * 2*floor(n/2)) = 9216 / 24 = 384      at n = 8 and n = 9;
> ```
>
> (b) writing `u = tau . eps_1` and `v = tau . eps_2`, an `h ∈ H` acts on `P`
> by `(0,j,k) -> (0, u(j), v(k))` if `pi` fixes the last two coordinates and by
> `(0,j,k) -> (0, u(k), v(j))` if `pi` swaps them; so it carries the row `p` to
> `p^h = v.p.u^{-1}` or `v.p^{-1}.u^{-1}`, which is again admissible.  Hence the
> shard universe is closed under `H`, which permutes it;
>
> (c) for every `h ∈ H` and every `p`, `h(S_p) = S_{p^h}` — the image of a shard
> is the whole of the image shard, record by record and cell by cell.  In
> particular `|S_p| = |S_{p^h}|`.

*Proof.*  (a)  The image's first coordinate is `tau(eps_0(x_{pi(0)}))`.  On `P`
the coordinates `x_1, x_2` are free, so if `pi(0) != 0` that expression takes
all `n` values on `P` and the image cannot be contained in the plane `x_0 = 0`;
hence `pi(0) = 0`,
and then the image's first coordinate is the constant `tau(eps_0(0))`, which
must be `0`.  Conversely every such `g` maps `P` into `P`, and an injection of a
finite set into itself is onto it.  `H` is the stabiliser of a subset, hence a
subgroup, and `[G:H]` is the size of the `G`-orbit of `P`.  That orbit is
`{ {x_a = c} }` with `a` any of the three axes and `c` any value of `tau(0)` or
`tau(n-1)`: since `tau` permutes the pairs `{t, n-1-t}` as blocks, `c` ranges
over every symbol except the middle one of an odd `n`, so over `2*floor(n/2)`
values.  Hence `[G:H] = 3 * 2*floor(n/2) = 24` and `|H| = 384`.  (`symmetry group` recomputes both by construction: it selects the elements of
`G` that fix `P`, exhibits the 24 planes, and forms all `384^2 = 147 456`
products of pairs of selected elements, requiring each to be **an element of the
selected set**, found by lookup.  Checking instead that a product preserves `P`
would check nothing, since a composition of two plane-preserving maps preserves
the plane whatever else it does.)

(b)  With `pi(0) = 0`, `pi` restricts to a permutation of `{1, 2}`, which gives
the two displayed forms.  The graph `{(j, p(j))}` is carried to
`{(u(j), v(p(j)))}`, the graph of `v.p.u^{-1}`, or to `{(u(p(j)), v(j))}`, the
graph of `v.p^{-1}.u^{-1}`.  For admissibility, note that `u` and `v` are each
`tau` or `tau.rev`, and that `tau` commutes with `rev`.  Take the first form and
put `s = u^{-1}(t)`.  Then `t` is a fixed point of `p^h` iff `v(p(s)) = u(s)`,
and a reflected point iff `v(p(s)) = rev(u(s))`.  Cancelling `tau` from both
sides — legitimate because `tau` is a bijection commuting with `rev` — turns
those two conditions into `p(s) = s` and `p(s) = rev(s)`, in that order when
`eps_1 = eps_2` and in the opposite order when not.  Each has exactly one
solution because `p` is admissible, so `p^h` has exactly one fixed and exactly
one reflected point.  The second form is the same computation after the
substitution `s = p^{-1}(u^{-1}(t))`.  So `p^h` is admissible, and
`h -> (p -> p^h)` is an action of `H` on the `48 912` admissible permutations.
(`symmetry orbits` checks all `48 912 × 384 = 18 782 208` images one by one —
every row against every subgroup element, not only the 157 representatives —
requiring each to be admissible and to be present in the brute-force shard
list.  It reports the number of images it checked, and stops if that is not
`shards × |H|`.)

(c)  `h` is a bijection of the cells with `h(P) = P`, so for any set `T`,
`h(T) ∩ P = h(T ∩ P)`.  Let `T` be a support whose row 0 is `p`.  By Lemma 4,
`hT` is a support, and its row 0 is `h(T ∩ P) = p^h` by (b).  Hence
`h(S_p) ⊆ S_{p^h}`.  Applying the same to `h^{-1}`, which is in `H` because `H`
is a group, gives `h^{-1}(S_{p^h}) ⊆ S_p`, i.e. `S_{p^h} ⊆ h(S_p)`.  The two
are therefore equal. ∎

At `n = 9` the `48 912` shards fall into **157** `H`-orbits, of sizes 384 (107
orbits), 192 (36), 96 (8), 48 (2) and 12 (4); at `n = 8` the `5 568` shards fall
into 25.  By (c) the number of supports in a shard is constant on an orbit, and
that is visible in the finished catalogue: across all 157 orbits, no orbit
contains two shards with different support counts.

`verify.sh symmetric` does exactly what (c) licenses: it enumerates the 157
representatives with the same exact-cover search `full` uses, writes every other
shard as the image of its representative's payload under a recorded element of
`H`, and then checks three things — that at `n = 8` the same construction from
25 enumerated shards reproduces `data/n8_supports.bin` set-equal *and*
byte-identical; that at `n = 9` a uniform sample of the *mapped* shards (320 by
default, drawn with a fixed seed from the 48 755 never enumerated) agrees
byte-identically with a direct re-enumeration; and that the assembled catalogue
has the canonical SHA-256, which is a statement about the whole set of
`14 616 576` supports and is where an error anywhere in the mapping would
surface.  A referee who distrusts Lemma 6 need not argue with it: `full` and
`symmetric` produce the same `1 183 942 656` bytes, and `full` never mentions
it.

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

This repository contains one support enumerator, `enum`.  Both `full` and
`symmetric` use it; their agreement checks the symmetry shortcut, not the
enumerator.  A catalogue produced by another implementation can be compared
against this one using the canonical byte order above and SHA-256, and the
author has such an agreement on record from an earlier constraint-solver
enumeration — but that implementation and its run records are not in this
repository, so nothing here recomputes it, and any claim of independent
historical agreement would need that separate provenance to be worth
anything.

### What each step does

Steps 3a-3d run only in `symmetric` mode; everything else runs in every mode
that reaches it.  Times below were measured by the runs recorded in this README
(see *How long it takes*), not estimated.  Every artefact's SHA-256 is in `checksums.txt` and
is checked by `verify.sh`.

| # | step | program | input | expected output |
|---|------|---------|-------|-----------------|
| 1 | order-8 controls | `tests/test_n8.py` | — | 37 checks pass (below) |
| 2 | shard universe | `shards 9 list` | — | 48 912 admissible rows = A007016(9) |
| 3 | enumerate `T(9)` | `enum 9 shards` | shard list | 48 912 shards EXHAUSTED, 14 616 576 supports |
| 3a | *(`symmetric` only)* the plane-fixing subgroup | `symmetry 9 group` | — | `\|H\| = 384`, index 24, all 147 456 products elements of `H` |
| 3b | *(`symmetric` only)* shard orbits | `symmetry 9 orbits` | shard list | 157 orbits, all 48 912 x 384 = 18 782 208 images admissible and present |
| 3c | *(`symmetric` only)* enumerate, then map | `enum 9 shards`, `symmetry 9 expand` | 157 representatives | the same 48 912 shards, 14 616 576 supports |
| 3d | *(`symmetric` only)* controls | `catalogue.py sample`, `enum 9 shards`, `catalogue.py setcmpshards` | 320 mapped shards | set-equal and byte-identical to direct enumeration |
| 4 | audit, assemble | `catalogue.py audit`, `pack` | manifest, payloads, shard universe | the universe covered exactly, records in canonical order, then `n9_supports.bin` |
| 5 | definition check | `check.py verify 9` | catalogue | 14 616 576 records, 0 failures against the 301 main lines |
| 6 | group closure | `orbits 9 closure` | catalogue | 87 699 456 images checked, 0 missing |
| 7 | roots | `orbits 9 centre` | catalogue | 11 821 056 supports through `(4,4,4)` |
| 8 | orbits | `orbits 9 classify` | roots | 2 049 orbits, sizes summing to 11 821 056 |
| 9 | companion pools | `pools 9 build` | catalogue, 2 049 reps | pool sizes 764 / 980 / 988.3 / 1 552 (min/median/mean/max) |
| 10 | pools re-derived | `check.py pools 9`, `catalogue.py reports` | catalogue, pools | 0 disagreements, and one complete report per worker covering all 2 049 queries |
| 11 | exhaust + cliques | `pack 9 roots` | 2 049 reps, pools | 2 049 EXHAUSTED, 0 covers, max depth 2, 183 625 nodes; max packing 3 (2 048 orbits) and 5 (1 orbit) |
| 11b | the results say what the theorem needs | `catalogue.py validate` | ledger | 2 049 queries, all EXHAUSTED, 0 covers, cliques ≤ 4, triangles < 56 |
| 12 | witness | `pack 9 witness` | exceptional orbit | 5 pairwise disjoint supports, re-verified |
| 13 | checksums | `check.py sha256`, `catalogue.py checksums` | all artefacts | the exact seven names, all matching `checksums.txt` |

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
byte; the raw output with wall times is left beside it as `n9_results_raw.jsonl`.

### The order-8 controls

Order 8 is where a cube exists and the census is known, so the whole pipeline is
exercised at an order where the answers were not produced by this run.
`./verify.sh test` runs 37 checks; the substantive ones are

* `T(8)` enumerated from nothing is **13 056** supports, set-equal *and*
  byte-identical to the shipped `data/n8_supports.bin`;
* all 13 056 pass the definition-level check against the 244 main lines, and the
  C and numpy constructions of the main lines agree as *sets of lines*, in
  canonical form — at `n = 8` and at `n = 9`;
* `T(8)` is closed under the order-9216 cell group and falls into **6** orbits,
  of sizes 768, 768, 2 304, 2 304, 2 304, 4 608;
* the exact cover of `[8]^3` by supports has exactly **198 624** solutions, with
  node counts by depth `1, 1 632, 26 016, 60 672, 69 513, 147 292, 154 660,
  198 624, 198 624`.  Relabelling the eight classes of a partition by two
  different permutations gives two different arrays, so this also counts the
  labelled cubes: `8! * 198 624 = 8 008 519 680` FDLHs of order 8;
* the **root search agrees with the census on every one of the 13 056 supports**
  — for each support, the number of covers containing it, computed once by the
  census and once by the same root machinery that is run at order 9;
* the packing ceiling and the search agree at order 8 too: a support lies in a
  cover **iff** its pool graph has a 7-clique;
* the whole of `verify.sh symmetric` is run at order 8: the plane-fixing
  subgroup has order 384 and index 24 and all `384^2` of its products are
  elements of it, the 5 568 shards fall into **25** orbits with every one of the
  `5 568 x 384` images admissible and present in the shard list, and mapping the
  5 543 non-representative shards from the 25 enumerated reproduces the census
  **byte for byte**;
* the three capped command paths — the sweep, the census and the root search —
  are each tested to fire, driven by an injected clock so that a cap fires after
  an exact number of readings rather than eventually; the shard-level cap is
  shown to leave no payload, not to mark the shard done, and to have it redone
  on the next pass; and `pack`'s cap is shown to fire in the clique stage as
  well as in the search, since it is a deadline on the whole query.  A stopping
  rule that never fires looks exactly like an exhaustion — the failure mode a
  negative result most needs excluded;
* and, in the same spirit, the mapping's own guard is tested to fire: one shard
  is pointed at a different element of the subgroup — in a full-size orbit, so
  the stabiliser is trivial and any other element does land the records
  elsewhere — and `symmetry expand` must refuse it.  A check that never rejects
  looks exactly like agreement.

The rest of the suite is failure paths, for the same reason.  A record byte that
is not a coordinate in `[n]` must be refused by all five programs that load
records, before it is used as an index; an order beyond the compiled capacities
must be refused before any construction; a pool whose members meet their own
query must be refused rather than searched; a killed and damaged sweep must
resume correctly (above).  The audit must reject an empty manifest, a shard
outside the universe, a malformed line, a missing payload digest, a missing
universe and records out of order; the report validator must reject an empty
report, a missing one and one covering too few queries; the checksum comparison
must reject a duplicated name, a malformed digest and a required artefact that
is not listed.  And the result validator, handed the order-8 ledger — the case
where a cube *does* exist — must refuse it.

### How long it takes

Two machines, both arm64.

* **A** — a 20-core Linux box (Cortex-X925 + Cortex-A725, 121 GB), 14 workers
  at `nice -n 10`, with other jobs of its owner's also running.
* **B** — a 16-core macOS laptop (64 GB), 6 workers at `nice -n 10`, kept
  deliberately below its core count so the machine stays usable.

| mode | A (20-core Linux, 14 workers) | B (16-core laptop, 6 workers) |
|---|---|---|
| `test` | 12 s | 97 s |
| `full` | **2 h 41 m 28 s** | **7 h 30 m to 13 h** *(projected — see below)* |
| `symmetric` | ≈ 9 m *(projected from A's own steps)* | **23 m 46 s** |
| `fast` | ≈ 7 m *(projected from A's own steps)* | **12 m 49 s** |

Every figure is a wall-clock measurement except the two marked as projected.
B's column is the run recorded here; two earlier sessions on the same laptop
gave 11 m 32 s and 15 m 42 s for `symmetric` and 6 m 21 s and 8 m 17 s for
`fast`, and A's row was measured in one of them.  Those earlier numbers predate
the checks this package now runs, which add about half a minute of order-8
controls and three seconds at order 9; the rest of that spread is the
instrument, not the work.

**A caveat on all of B's absolute numbers.**  A laptop under sustained all-core
load is not a stable instrument.  The 320-shard control below is the same 320
shards and, as its identical node count confirms, exactly the same search every
time; it cost 3.241 core-seconds per shard in one measurement session, 4.960 in
another and 5.665 in the one tabulated here, and every other step moved with
it.  A pair of single-worker runs of the same 24 shards, one with this code and
one with the code as it stood before the checks were strengthened, came out at
115.6 and 119.8 core-seconds, so the drift is the machine.  Ratios measured
*within* one run are unaffected, which is why the speed-up below is quoted that
way.

**How B's `full` row is projected.**  A uniform sample of **1 024 of the 48 912
shards** was swept on B (the shard order is shuffled by a fixed seed, so a
prefix is a uniform sample, not a structurally poor one): **562 s wall on 6
workers**, 3 258.8 core-seconds, 16.29 x 10^9 search nodes, 301 942 supports,
nothing near the cap.  That is **3.182 core-seconds per shard**, so the whole
sweep is 155 700 core-seconds — 43.2 core-hours, or 7 h 27 m wall on 6 workers
— plus the few minutes everything downstream takes.  Two checks on the
extrapolation: the sample's node count scales to 7.78 x 10^11 against A's
measured 7.83 x 10^11 (0.7 % apart), and its mean supports per shard is 294.9
against the true 298.8.  That sample was taken in the fastest of the three
sessions, so **7 h 30 m is a floor**; at the per-shard cost of the slowest it
would be 13 h.

**Where the time goes** (B, the two runs recorded here):

| step | `symmetric` | `fast` |
|---|---:|---:|
| order-8 controls, 37 checks | 157 s | 102 s |
| shard universe (brute force over 9!), twice, and the shard orbits | 6 s | — |
| control (a): the n = 8 census, symmetrically | 8 s | — |
| enumerate the **157** orbit representatives | 236 s (**1 295 core-s**) | — |
| map the other **48 755** shards | 16 s | — |
| control (b): 320 mapped shards re-enumerated directly | 326 s (1 813 core-s) | — |
| audit 48 912 shards, assemble and hash the catalogue | 46 s | — |
| SHA-256 of the supplied 1.1 GiB catalogue, before anything reads it | — | 2 s |
| all 14 616 576 records checked against the definition in numpy | 34 s | 35 s |
| closure under the order-9216 cell group, 87 699 456 images | 66 s | 62 s |
| roots, and their classification into 2 049 orbits | 64 s | 52 s |
| 2 049 companion pools, and their re-derivation in numpy | 465 s | 511 s |
| exhaust 2 049 root cases and compute their cliques | < 1 s | 3 s |
| the 5-packing, the validators and the checksums | 1 s | 2 s |

**The speed-up.**  The exact part of it is a count, not a time: the 157
representatives cost **2 512 814 852** search nodes against A's measured
7.83 x 10^11 for all 48 912, a factor of **312** — the reduction factor
48 912 / 157 almost exactly, as Lemma 6 leads one to expect.  Lemma 6 equates
the number of supports in symmetry-related shards, not their search effort, so
that agreement is a measurement rather than a consequence.  In wall-clock terms,
comparing within a single run — the representatives against the per-shard cost
of the 320-shard control drawn from the same sweep, at the same moment — the
three sessions give **326**, **313** and **214**.  Counting both controls as
part of the price, and they should be counted since they are what makes the mode
believable, the factor is **89** to **103**.  End to end, `symmetric` took 1 426 s
against `full`'s projection of 46 900 s on the same laptop in the same session:
**33x**.

Two thirds of `fast` is the numpy re-derivation of every companion pool, which
is bound by memory bandwidth and page cache, and it is where most of the
run-to-run spread lives (258 s in the fastest session against 511 s here).  The
shape of the cost is worth stating plainly: **the enumeration is the whole
expense and the refutation is free.**  Exhausting all 2 049 root cases costs
1.70 CPU-seconds and the independent clique computation another 4.5; building
the object they run over costs 35-43 core-hours.

### Memory and disk

Peak resident set per process, measured on B with `/usr/bin/time -l`:

| process | peak RSS |
|---|---:|
| `shards`, `pack`, `symmetry orbits`, `symmetry expand` | 3-31 MB |
| `enum` — one sweep worker, in any mode | **18 MB** |
| `catalogue.py pack` / `audit` | 79 MB |
| `check.py verify` — one worker | 1.30 GB |
| `orbits closure` / `orbits classify` | 1.26 GB / 1.14 GB |
| `pools build` | 2.41 GB |
| `check.py pools` — one worker | **4.75 GB** |

So the enumeration, the long part, is free of memory pressure: 18 MB per worker
whatever the worker count.  The ceiling is the numpy pool re-derivation.  Of its
4.75 GB the 1.1 GiB catalogue is memory-mapped and shared between workers; the
remaining ~3.6 GB — one 12-word cell bitmask per support (729 cells need 12
64-bit words) and the block temporaries that build them — is private to each.  A
run wants roughly `3.6 GB x pool workers + 1.2 GB`: about 23 GB at 6, which is
why the laptop measurements were taken on a 64 GB machine.

**`--pool-workers N` sizes that stage on its own** (default `min(--workers, 6)`).
Two pool workers want about 8.4 GB before the operating system is allowed
anything and three about 12 GB, so **on an 8 GB machine pass
`--pool-workers 1`**.  Turning `--workers` down instead also works, but it
controls the enumeration too, so it buys the memory back at a cost measured in
hours rather than minutes.

Disk under `--work`: about 3.4 GB for `full` or `symmetric` (1.2 GB of shard
payloads, the 1.1 GiB catalogue, the 957 MB root file, 160 MB of pools).  `fast`
generates about 1.1 GB (the root file and the pools) **excluding** the supplied
catalogue; counting the catalogue it needs about 2.3 GB, plus filesystem
overhead.

### Reproducing it

The four commands are at the top of this file.  `--workers N` defaults to the machine's core count; `--pool-workers N` sizes the
one memory-hungry stage separately (see *Memory and disk*); `--nice N` keeps the
machine usable; `--work DIR` puts the scratch somewhere else; `--cap SEC` is the
per-shard wall-clock cap in the sweep (default 3 600 s; the slowest shard
observed anywhere is 26.0 s, and that on a throttled laptop).  In `enum` the cap
bounds the search, which is essentially the whole of a shard's work; in `pack`
it is a deadline on the whole query, consulted by the graph construction and the
clique bound as well as the exact cover.  Setting `FDLH_CLOCK_STEP` to a
positive number of seconds replaces the monotonic clock with a virtual one that
advances by exactly that much per reading, which is how the caps are tested.

Both sweeping modes are **resumable per shard**, at any worker count: each
worker appends to its own manifest, and the manifests are merged and handed to
every worker at the start of the next run.  A payload is written to a `.part`
name, flushed and renamed before its manifest line is appended, so an interrupt
leaves no half-written shard, and a stray `.part` is removed before its shard is
redone.  A shard is skipped only if **both** its manifest record has the exact
shape a writer produces **and** its payload is still on disk at the recorded
length and FNV-1a digest; otherwise it is redone and the reason is printed.  A
manifest line torn by a kill is discarded before the next append, so two records
can never be joined; a malformed line anywhere *else* is an error.  A shard that
hits `--cap` is reported, is not marked done, and leaves no payload; if any
shard is unfinished when the sweep ends, the run stops with the list rather than
assembling a catalogue quietly short of a shard.  Once 200 shards — or a quarter
of them — have finished, the sweep prints its throughput and an ETA.

Every parallel stage keeps its workers' process ids and requires each to exit 0,
and every stage that writes one report per worker requires exactly one complete
report per worker and requires the reports between them to cover the whole
universe of work — all 2 049 query ids for the pool re-derivation, all
14 616 576 records for the definition check.  A killed or out-of-memory worker
is a failure, not a silence.

Requirements: a 64-bit POSIX platform, and GCC or Clang.  The C is C99 plus
POSIX 2001 (`clock_gettime(CLOCK_MONOTONIC)`, `mkdir`, `ftruncate`, `unlink`)
and the GCC/Clang builtins `__builtin_popcountll` and `__builtin_ctzll`; the
Makefile asks for `-std=c99 -D_POSIX_C_SOURCE=200809L`.  A platform without
`CLOCK_MONOTONIC` will not build.  Also GNU make, `bash`, and **Python 3.8 or
later** (`math.comb`) with **numpy**, the only Python dependency.  If numpy
lives in a virtualenv, set `PYTHON=/path/to/python`.  Tested on macOS 15 arm64
with Apple Clang 21 and Python 3.14.5 / numpy 2.5.1, and on Linux arm64 with
GCC.

`fast` checks the supplied catalogue's SHA-256 **before** anything reads it and
refuses to continue on a mismatch, so a wrong or truncated catalogue cannot
propagate into the downstream results; `full` and `symmetric` check the same
hash on the catalogue they have just assembled, before anything downstream looks
at it.

Every binary takes `--help`.  Nothing in this repository downloads anything or
sends anything anywhere.  All output goes under `--work` (default `./work`),
with two exceptions: `make` writes the six binaries into `bin/`, and the order-8
suite makes its scratch directory under `--work` when `verify.sh` runs it and
under the system temporary directory when run on its own.

---

## What a referee still has to take on trust

The point of Part I is that the reduction from "is there a cube of order 9" to
"here are 2 049 finite graphs, count their triangles" is checkable by hand.
What remains is:

The computational premise is that the catalogue is exhaustive, that the
representative list covers every root orbit, and that each companion graph is
constructed and bounded correctly.  Checking every stored record establishes
validity, not completeness; closure under the group does not exclude a whole
missing orbit; and a matching checksum establishes agreement with reference
bytes rather than exhaustion.  What each mode does about that:

* `full` establishes completeness by exhaustive search over every admissible
  row.  `symmetric` uses exhaustive representative searches together with
  Lemma 6.  `fast` **assumes** the supplied catalogue is complete, checking its
  reference digest before performing anything downstream.
* In every mode the sweep is audited against the shard universe: the manifest
  must cover it exactly, every record must carry its own shard's row 0, and the
  records must be in strict lexicographic order within each shard and the
  shards in strict row order between them.  That is an exhaustion audit rather
  than a digest comparison, and it explains a discrepancy instead of merely
  detecting one.

1. **That the programs implement the definitions.**  The main lines are built
   twice — in C from the pattern description, and independently in numpy — and
   the two are compared as *sets of lines*, in canonical form, not by their
   counts; every record ever produced is checked against the numpy version of
   the definition.  But the definition itself is written down twice by the same
   author, and a reader who disagrees with the definition of a main line will
   disagree with everything downstream.  Read `src/lines.h` and `main_lines()`
   in `src/check.py`; they are about twenty lines each.

2. **That the exhaustive enumeration is exhaustive.**  At order 8 the identical
   pipeline reproduces a census whose values — 13 056 supports, 198 624 covers,
   the node profile, the six-orbit decomposition — are supplied to the suite as
   regression targets rather than derived by it; they are the author's, quoted
   from an earlier run, so they check that the pipeline still computes what it
   computed and are not an independent authority.  And the enumerated `T(9)` is
   *closed under the order-9216 group* (87 699 456 images checked, none
   missing), so a missing support would have had to be missed together with its
   whole orbit.

3. **Common-mode error in the model.**  Both the exact-cover enumerator and the
   definition-level checker are built on the same reading of "main line" and of
   the record format.  The Python checker is a separate implementation of the
   same definitions, checking data the C programs produced; the two packing
   algorithms share their representation code and their input.  These checks
   corroborate; they do not eliminate a common error, and the argument still
   rests on correct compilation and execution.

4. **Floating point plays no role**; there is none in the argument.  Wall-clock
   caps are the only timing-dependent behaviour, and they are configured to
   values the sweep never approaches (the slowest shard measured anywhere takes
   26.0 seconds against a 3 600-second cap), and are separately tested to fire —
   against an injected clock, so the test does not itself depend on wall time.

5. **Which mode was run.**  `full` and `fast` rest on Lemmas 1-5 only.
   `symmetric` additionally rests on Lemma 6 and on `src/symmetry.c`
   implementing it: a wrong subgroup, element index or record image would mean
   48 755 of the 48 912 shards came from an untrusted map rather than a search.
   Three things stand against that — the n = 8 census reproduced the same way,
   320 mapped shards re-enumerated and found byte-identical, and the canonical
   SHA-256 of the whole catalogue — of which the first two are samples and the
   third is not.  A reader who wants no lemma beyond 1-5 should run `full`.

6. **The scope.**  This computation establishes nonexistence at order 9 in
   dimension 3.  It does not determine existence at order 10 or order 12, says
   nothing about dimension 4 or above, and establishes no general existence
   threshold.

---

## Files

```
README.md            this file
Makefile             builds the six C programs into bin/
verify.sh            full / symmetric / fast / test
checksums.txt        SHA-256 of every data artefact
src/lines.h          the main lines of [n]^3, from the definition
src/enum.c           the enumerator: exact cover by dancing links
src/shards.c         the admissible row-0 permutations (A007016)
src/orbits.c         the order-9216 cell group; closure and orbit classification
src/pools.c          companion pools as a filter over a complete catalogue
src/pack.c           exhaustive exact cover, and the packing ceiling
src/symmetry.c       the plane-fixing subgroup, and the shard orbits (Lemma 6)
src/util.h           record loading and byte-range checks, the clock (real or
                     injected), payload digests and manifest recovery
src/check.py         definition-level checks in numpy, independent of the C
src/catalogue.py     assembly, the exhaustion audit, ledger canonicalisation,
                     and the strict validators for results, worker reports and
                     the checksum manifest
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

## Provenance

The programs, the verification script and this README were written with
Claude Code (Anthropic), working from the mathematics above; the lemmas, the
checks and the checksums were reviewed by a human before being published, and
the point of the three verification modes is that nothing here need be taken
on trust from either.

## References

* W. Taylor, *On the coloration of cubes*, Discrete Mathematics **2** (1972)
  187–190.  The objects, the `n >= 2^d` corner bound, the conjecture, and
  Problem 2 (`d = 3`, `n = 9`).
* D. E. Knuth, *Dancing links*, in J. Davies, B. Roscoe, J. Woodcock (eds.),
  *Millennial Perspectives in Computer Science*, Palgrave (2000), 187–214.
  Algorithm X and the doubly linked cover matrix used by `src/enum.c`.
* [OEIS A007016](https://oeis.org/A007016) — permutations with exactly one
  fixed point and exactly one reflected point; the shard count.
