# There is no fully diagonalised Latin cube of order 9

This repository aims to cleanly demonstrate a negative answer to a question which was (to my knowledge) first posed by Walter Taylor in a 1972 paper.

A d-dimensional **fully diagonalised Latin hypercube** (FDLH) (called *completely Latin* by Arkin, Hoggatt and Straus) is a coloring $A : [n]^d \to [n]$, where $[n] = \lbrace0, 1, \ldots, n-1\rbrace$, in which each color occurs exactly once on every *line*.  A line is obtained by letting $t$ run over $[n]$ and taking each coordinate to be a constant, $t$, or $n-1-t$, with at least one coordinate varying.  In dimension 3 the lines are the rows, the columns and the pillars, the diagonals of every planar cross-section of the cube, and the four space diagonals; each must therefore contain every color.

Every pair of the $2^d$ corner cells lies on a line, which forces $n \le 1$ or $n \ge 2^d$ (Taylor's Proposition 4).  Taylor's 1972 Problem 3 asks two questions — in his notation, where P(m, n) asks whether there exists an n-cube of order m: "For every $n$ does there exist $M$ such that $P(m, n)$ whenever $m \ge M$?  May one take $M = 2^n$?"  His Problem 2 asks in particular about order 9 in dimension 3. 


(His Problem 1, the 12x12 square, was settled almost immediately: Hilton (1973) and, independently, Faber constructed doubly diagonalized Latin squares of every order at least 4, as Taylor notes in proof, and Gergely (1974) gave a simpler construction.)


While there is a fully diagonalized Latin cube of order 8 (as shown by Taylor), **there is none of order 9.**  This answers Problem 2 negatively and rules out $M = 2^d$ at $d = 3$; it leaves open the first question of Problem 3, whether *some* $M$ exists in each dimension.  The method does not extend to order 10: the catalogue of candidate color classes would be thousands of times larger (a pilot suggests at least 6 x 10^10 entries and a core-year of enumeration), and the case analysis here relies on the cube having a center cell, which an even cube does not (though the order 8 case is sufficiently small to use as a sanity check anyways).  Settling order 10 seems to need a new idea, though if an order-10 cube exists a search might of course find it long before exhausting anything.


The rough outline of this algorithm is as follows:
1. Enumerate every *support*, a set of cells that could be the locations of a single color within an FDLH (that is, sets containing precisely one cell per line).
2. Track those supports containing the center cell, one per symmetry class. 
3. Form the graph of which of them can occur in the same FDLH — two vertices share an edge if they represent disjoint supports.  A 9x9x9 would need eight mutually compatible companions to some fixed center support — that is, an 8-clique in that support's graph.  Any 8-clique contains 56 triangles (mutually compatible triples), so counting triangles is enough: for 2 048 of the 2 049 center supports the graph has no triangles at all, and for the last one it has 8.  No graph can contain an 8-clique, so no center support extends to a 9x9x9, and every 9x9x9 would have to contain one.


{"Five short lemmata justify that reduction" the reducion seems self-evidently justified to me. Whatever the lemmata justify, it isn't the rough outline above, which doesn't really require justification. Perhaps they justify that our operationalization of "one per symmetry class" is valid? That seems like the one step that actually needs much elaboration.}


{"Part II states the machine's part as a table of programs, inputs,
expected outputs and SHA-256 checksums, so a reader can rerun any line
independently.  *From nothing* means regenerating the order-9 catalogue rather
than supplying a precomputed one.  A fifth lemma, proved afterwards, proves
nothing about order 9: it justifies a symmetry shortcut used by one
verification mode." i probably want to rewrite this but i need to go see what's actually in part 2 first.}

Start with `make` and `./verify.sh test`; then pick `full`, `symmetric` or
`fast` according to how much you want to recompute.

```
make                      # six C programs, no libraries
./verify.sh test          # the order-8 controls          (12-97 s)
./verify.sh full          # everything, from nothing    (2 h 41 m on 14 cores)
./verify.sh symmetric     # the same, one shard per symmetry orbit (9-24 min)
./verify.sh fast --catalogue n9_supports.bin   # everything downstream of the
                                               # catalogue     (7-13 minutes)
```

You need a 64-bit POSIX platform with GCC or Clang, GNU make, `bash`, and
Python 3.8+ with numpy (details under *Reproducing it*).  The pool re-derivation
is the one memory-hungry stage: on an 8 GB machine pass `--pool-workers 1`.

There are three modes.  `full` recomputes the catalogue from nothing and uses no mathematics beyond the clique bound (Lemmata 1-4).
`symmetric` recomputes it too, but enumerates only 157 of the 48,912 shards (sets of supports sharing the same first face) and obtains the rest by symmetry — three hundred times less search, at the cost of also relying on the shard symmetry (Lemma 5).  `fast` takes the catalogue as
given, verifies its SHA-256 before reading it, and redoes everything downstream.

A failed check aborts the run, so reaching the final `done --` line means every
step passed, the last of them comparing all seven SHA-256s against
`checksums.txt`.  The result itself is then `work/n9_root_results.jsonl`, one
line per root orbit: `status` is `EXHAUSTED` and `covers` is `0` on all 2 049.

---

## Part I. The mathematics

Throughout, $n$ is the order — $9$ in every concrete example, $8$ in the
controls — $r(t) = n-1-t$, and cells of $[n]^3$ are written $(i,j,k)$ and
indexed (in C) by $(i \cdot n + j) \cdot n + k$, which at $n = 9$ is $81i + 9j + k$.

### Lemma 1 (supports and partitions)

A line has $n$ cells, so it carries all $n$ colors exactly when it
carries each of them once.  Call a set of cells meeting every line exactly
once a **support**, and write $T(n)$ for the set of supports of $[n]^3$.  The
*color classes* $A^{-1}(\text{color})$ of an FDLH are therefore supports, and an FDLH of
order $n$ exists **iff** $T(n)$ contains $n$ pairwise disjoint members.

*Proof.*  Every support has exactly $n^2$ cells: the $n^2$ pillars $\lbrace(i,j,\ast)\rbrace$
are lines (only the third coordinate varies), they are pairwise disjoint,
they cover the cube, and a support meets each once.  So $n$ pairwise disjoint
supports occupy $n \cdot n^2 = n^3$ cells and therefore partition $[n]^3$.  Give
each a different color; every line then carries each color exactly once,
so the array is an FDLH.  Conversely the color classes of an FDLH partition
the cube and meet every line once, so each is a support. ∎

**Listing and counting lines.**  Replacing $t$ by $r(t)$ throughout traces
the same line in the opposite direction, so to list each line once, require its
first varying coordinate to use the pattern $t$.  In dimension $d$, choose the
$s$ varying coordinates, their $2^{s-1}$ sets of directions ($t$ or $r(t)$ for each varying coordinate besides the first), and the fixed values of the
$d-s$ others:

$$\sum_{s=1}^{d} \binom{d}{s}\, 2^{s-1}\, n^{d-s},$$

which for $d = 3$ is $3n^2 + 6n + 4$: **301** at $n = 9$ (243 axis lines, 54
plane diagonals, 4 space diagonals), 244 at $n = 8$.

### Lemma 2 (Latin-square form, and what row 0 must be)

Every support $T$ of $[n]^3$ must be of the form

$$T = \lbrace(i, j, L(i,j)) : i, j \in [n]\rbrace$$

for a unique Latin square $L : [n]^2 \to [n]$, and its first row $p = L(0, \cdot)$
must have **exactly one fixed point** and **exactly one reflected point**: exactly one
solution each to $p(t) = t$ and $p(t) = r(t)$.


Supports have the form $\lbrace(i, j, L(i,j)) : i, j \in [n]\rbrace$ merely due to the requirement that the support contain one cell from each *axis* line — that part of the statement would be true even for supports of merely Latin cubes (i.e. the more general notion of support that permits choosing zero or multiple cells from a diagonal).  The fixed point comes from the requirement that the support hit the positive diagonal of the plane $i = 0$, while the reverse point comes from the requirement to hit the negative diagonal:


*Proof.*  The pillar $\lbrace(i,j,\ast)\rbrace$ is a line (third coordinate $t$, the other
two constant), and $T$ meets it exactly once, which picks out the single value
$L(i,j)$.  Meeting each of the lines $\lbrace(i,\ast,k)\rbrace$ and $\lbrace(\ast,j,k)\rbrace$ exactly
once says that every row and every column of $L$ contains every value once, so
$L$ is a Latin square and its first row $p$ is a permutation.  Finally the plane $i = 0$
contains the two lines $\lbrace(0,t,t)\rbrace$ and $\lbrace(0,t,r(t))\rbrace$ — first coordinate
constant at $0$, the other two varying — and meeting each exactly once gives
the two conditions on $p$. ∎

Call a permutation with one fixed point and one reflected point **admissible**.  A **shard** is the set of all supports whose Latin squares $L$ have one specified first row, or equivalently the set of all supports that share one specified $i=0$ plane. Since every support's first row is admissible (Lemma 2), an exhaustive enumeration need only consider the shards corresponding to every admissible row (some of which may turn out to be empty).


The numbers of admissible permutations are [OEIS A007016](https://oeis.org/A007016): $8, 20, 96, 656, 5568, 48912$ for $n = 4..9$.  This repository computes the order-9 count in two different ways, by brute force over all $9! = 362\,880$ permutations (`shards.c`) and by inclusion–exclusion in closed form (`check.py a007016 N`), and both agree with OEIS.

Within each shard the enumerator solves an exact-cover problem: the 301 lines
must each be covered exactly once, and choosing a cell covers the lines through
it.  Nine cells are fixed by $L$'s first row / the support's $i=0$ plane.  The search is [Knuth's Algorithm X with dancing links](https://en.wikipedia.org/wiki/Knuth's_Algorithm_X).

### Roots

When $n$ is odd, $[n]^3$ has a center cell $c_0 = (m,m,m)$ with $m = (n-1)/2$
— $(4,4,4)$ at $n = 9$ — and exactly one color class of an FDLH of order $n$
contains it, because the color classes partition the cube.  Call a support
containing $c_0$ a **root**.  To rule out an FDLH of odd order $n$ it therefore
suffices to show that no root belongs to a partition of $[n]^3$ into $n$
supports.  (An even cube has no center cell and no roots; the order-8 control
runs its root search on every support instead, see below.)

An aside, not needed for the proof: all four space diagonals pass through $c_0$. A root meets all four with that one cell, while a support avoiding the center must meet them at four distinct cells.  That may be why roots are more common than one might naively expect — $11\,821\,056$ of the $14\,616\,576$ supports, 80.9 %, contain the center.

### Lemma 3 (orbit reduction)

We need test only one root from each symmetry class, for some suitable notion of symmetry classes. For example, if two roots are reflections of each other, ruling out that one of them occurs in an FDLH also implies that the other does not. To maximize the efficiency of our exhaustive enumeration, we wish to use the largest valid choice of symmetry group we can, i.e. the group which maximizes the size of each root's symmetry class while still ensuring that checking one representative from each class is sufficient.

The symmetries used here permute the axes, reverse individual axes, and relabel the coordinate *values* by one permutation applied on every axis at once.  Axis lines survive any
relabelling, but a diagonal such as $(t, 8-t, 3)$ survives only if the
relabelling respects the pairing $u \leftrightarrow 8-u$: relabel $0 \leftrightarrow 1$ alone and the
cells $(0,8,3), (1,7,3), \ldots$ become $(1,8,3), (0,7,3), \ldots$, which lie on no line.  So
the relabellings used are those that shuffle the pairs $\lbrace0,8\rbrace, \lbrace1,7\rbrace, \lbrace2,6\rbrace, \lbrace3,5\rbrace$ as blocks and optionally flip each, leaving $4$ fixed — the permutations
that commute with reversal.

Precisely, let $C(r)$ be the centraliser of $r$ in $S_n$, the symmetric group
on $[n]$: that is, the set of elements of $S_n$ that communte with $r$. $C(r)$ consists of the permutations $\tau$ of the coordinate values such that $\tau(r(t)) = r(\tau(t))$ for all $t$.  ($r$ is written `rev` in `src/orbits.c`.)  

For $\pi \in S_3$, $\varepsilon \in \lbrace\mathrm{id}, r\rbrace^3$ and $\tau \in C(r)$ define a map $g \colon [n]^3 \to [n]^3$ on cells $x = (x_0, x_1, x_2)$ by

$$g(x)_i = \tau\bigl( \varepsilon_i\bigl( x_{\pi(i)} \bigr) \bigr).$$

These form a group $G$ (a subgroup of $`S_{n^3}`$).  The claim is that $G$ 
(i) permutes the lines, hence maps supports to
supports; (ii) fixes the center cell (when n is odd); and (iii) maps partitions to partitions.

Consequently a root $T$ lies in a partition iff $gT$ does for some $g$, and it suffices to
test one root per $G$-orbit.  At $n = 8$ and $n = 9$,
$\lvert C(r)\rvert = 2^{\lfloor n/2\rfloor} \cdot \lfloor n/2\rfloor! = 384$ and
$\lvert G\rvert = 6 \cdot 8 \cdot 384 / 2 = \mathbf{9216}$, the $/2$ being proved with the closure below.

*Proof.*  (i)  Describe a line by its triple of patterns (constant $c$,
$t$, or $`r(t)`$) per coordinate.  $\pi$ permutes which coordinate carries which
pattern, leaving at least one non-constant.  $\varepsilon_i = r$ sends a constant $c$
in position $i$ to the constant $r(c)$, the pattern $t$ to $r(t)$, and $r(t)$
to $t$.  Applying the same $\tau$ to all three coordinates: take a line with, say, one
constant coordinate and both varying patterns, and follow it through:

$$\begin{aligned}
L   &= \lbrace (c,\ t,\ r(t)) : t \in [n] \rbrace \\
\tau L  &= \lbrace (\tau(c),\ \tau(t),\ \tau(r(t))) : t \in [n] \rbrace \\
    &= \lbrace (\tau(c),\ \tau(t),\ r(\tau(t))) : t \in [n] \rbrace \quad\text{since } \tau \text{ commutes with } r \\
    &= \lbrace (\tau(c),\ s,\ r(s)) : s \in [n] \rbrace \quad\text{writing } s = \tau(t);\ \tau \text{ is a bijection, so } s \text{ runs over all of } [n]
\end{aligned}$$

The result is again a line: the constant went to a constant and the varying
coordinates still carry the patterns $s$ and $r(s)$.  The middle step is the
only place the commuting condition is used — without it the third coordinate
would be $\tau(r(\tau^{-1}(s)))$, which is neither $s$ nor $r(s)$, and $\tau L$ would
not be a line.  Any other line is the same computation with a different choice
of pattern per coordinate.  In each case the image is again a line.

A bijection of cells that permutes the lines carries a set meeting every
line once to a set meeting every line once, so $G$ maps $T(n)$ to $T(n)$.

(ii)  For odd $n$, $r$ has the unique fixed point $m = (n-1)/2$, and
$\tau r = r \tau$ forces $\tau$ to permute the fixed points of $r$, so $\tau(m) = m$; also
$\varepsilon_i(m) = m$; and $\pi$ merely permutes the equal coordinates of $(m,m,m)$.
So $g(c_0) = c_0$ for every $g \in G$.  (For even $n$ there is no fixed point and
no center; (i) and (iii) are unaffected.)

(iii)  $g$ is a bijection of the cell set, so $T_1, \ldots, T_n$ are pairwise
disjoint supports covering the cube iff $gT_1, \ldots, gT_n$ are.

Finally, that these maps form a group, and how many of them there are.

*Closure.*  If $g = (\pi, \varepsilon, \tau)$ and $h = (\rho, \delta, \sigma)$, then
$g \circ h$ (apply $h$ first) is the map with parameters

$$\begin{aligned}
\pi'   &= \rho \circ \pi, \\
\varepsilon'_i &= \varepsilon_i \circ \delta_{\pi(i)}, \\
\tau'  &= \tau \circ \sigma.
\end{aligned}$$

To see this, substitute $`h(x)_j = \sigma(\delta_j(x_{\rho(j)}))`$ into the
definition of $g$, with $j = \pi(i)$:

$$\begin{aligned}
(g \circ h)(x)_i &= \tau\bigl(\varepsilon_i\bigl(h(x)_{\pi(i)}\bigr)\bigr) \\
  &= \tau\bigl(\varepsilon_i\bigl(\sigma\bigl(\delta_{\pi(i)}\bigl(x_{\rho(\pi(i))}\bigr)\bigr)\bigr)\bigr) \\
  &= \tau\bigl(\sigma\bigl(\varepsilon_i\bigl(\delta_{\pi(i)}\bigl(x_{\rho(\pi(i))}\bigr)\bigr)\bigr)\bigr)
     \quad\text{since } \varepsilon_i \in \lbrace\mathrm{id}, r\rbrace \text{ commutes with } \sigma \in C(r) \\
  &= (\tau \circ \sigma)\Bigl((\varepsilon_i \circ \delta_{\pi(i)})\bigl(x_{(\rho \circ \pi)(i)}\bigr)\Bigr).
\end{aligned}$$

The third line is the one step that needs anything: a reversal has to be moved
past a value permutation, and that is legitimate exactly because $\sigma$
commutes with $r$.  The result is of the displayed form, with $\tau \circ \sigma$
in $C(r)$ because $C(r)$ is a subgroup and each $\varepsilon_i \circ \delta_{\pi(i)}$
in $\lbrace\mathrm{id}, r\rbrace$ because $r \circ r = \mathrm{id}$.
The identity is $(\mathrm{id}, \mathrm{id}, \mathrm{id})$, and a finite composition-closed family of
bijections of a finite set contains inverses.  So $G$ is a group.

*The count.*  Suppose two parameter triples define the same cell map.  Each
output coordinate depends bijectively on exactly one input coordinate, so the
two coordinate permutations agree; write $\pi$ for both.  Equating the $i$-th
output coordinates gives $\tau(\varepsilon_i(x)) = \tau'(\varepsilon'_i(x))$ for all $x$, i.e.
$\tau'^{-1} \tau = \varepsilon'_i \varepsilon_i^{-1}$ for every $i$.  The right-hand side is
$\mathrm{id}$ or $r$, and it is the *same* element for every $i$ because the left-hand
side does not depend on $i$.  So either all three reversal bits agree and
$\tau = \tau'$ — the triples are identical — or all three are flipped and
$\tau' = \tau \circ r$.  Every map therefore has exactly two descriptions, and
$\lvert G\rvert = 6 \cdot 8 \cdot \lvert C(r)\rvert / 2$.

A $\tau$ commuting with $r$ permutes the $\lfloor n/2 \rfloor$ pairs $\lbrace t, r(t)\rbrace$ as
blocks and may flip each, and fixes the middle point when $n$ is odd, so
$\lvert C(r)\rvert = 2^{\lfloor n/2\rfloor} \cdot \lfloor n/2\rfloor!$ — $384$ at $n = 8$ and $n = 9$, so
$\lvert G\rvert = 9216$.  The program `orbits` builds all
$6 \cdot 8 \cdot 384 = 18\,432$ maps, deduplicates them, and checks that $9216$ remain,
that all of them fix $c_0$, and that six explicitly named generators generate
all of them. ∎

### Lemma 4 (the clique bound: what actually kills order 9)

Completing a support $T$ to an FDLH requires $n-1$ further supports, disjoint
from $T$ and from one another.  Let $P(T) = \lbrace S \in T(n) : S \cap T = \emptyset \rbrace$ be the
**companion pool** of $T$, and let $\Gamma(T)$ be the graph on $P(T)$ joining
two companions when they are disjoint.  Because the catalogue is complete,
every possible companion is in the pool, so $n-1$ such supports are exactly a
clique of size $n-1$ in $\Gamma(T)$, and any such clique together with $T$ is
$n$ pairwise disjoint supports, hence a partition (Lemma 1).  Hence, for every $n$,

$$T \text{ belongs to a partition} \iff \Gamma(T) \text{ contains a clique of size } n - 1,$$

and in particular $\omega(\Gamma(T)) \lt n - 1$ rules $T$ out (here $\omega$ is the
*clique number*, the size of the largest clique).  At $n = 9$ a bound of
$\omega(\Gamma(T)) \le 4$ already suffices, and a *triangle-free* $\Gamma(T)$ ($`\omega \le 2`$)
suffices very comfortably.  At $n = 8$ the control checks the equivalence
itself: a support lies in a cover iff its pool graph has a 7-clique.

The recorded triangle counts alone rule out every root case, with no clique
computation at all.  An 8-clique contains $C(8,3) = 56$ triangles, whereas
every one of the 2 049 computed companion graphs has either $0$ or $8$.  None
of them can contain an 8-clique, so no root lies in a partition, which would
require one (Lemma 4).
Counting triangles is a cheap way to rule out 8-cliques: `pack` builds the
adjacency bitmap of each companion graph, and for every edge $ab$ counts the
common neighbours of $a$ and $b$ with a popcount over the two rows' AND;
summing over edges counts each triangle three times.  `catalogue.py validate`
checks the recorded counts against the $56$ threshold.

`pack` also computes each maximum clique outright, by the recursive branch and
bound `bk()`, and gets 2 for 2 048 of the orbits and 4 for one.  That is a
different algorithm from the exact cover — but the two live in the same program
and share its record loader, its cell masks and its `disjoint()`, so they are
independent as *algorithms*, not as implementations.

### What the order-8 control checks

Order 8 is where a cube exists and the census (Part II) is known
independently, so every lemma above is exercised at an order whose answers
this code did not produce.  Each check is an instance of a lemma at $n = 8$:

* the catalogue is built from the $5\,568$ admissible first rows — Lemma 2's
  shard universe at $n = 8$ — and every record passes the definition of a
  support against the $244$ lines (Lemma 1's definition);
* $T(8)$ is closed under $G$ and splits into six orbits — Lemma 3 (i);
  clause (ii) is not used, since there is no center;
* the exact cover of $[8]^3$ by supports has $198\,624$ solutions, i.e.
  $8! \cdot 198\,624$ labelled cubes — the correspondence of Lemma 1;
* for every one of the $13\,056$ supports, the root search finds exactly the
  number of covers through it that the census predicts, and a support lies in
  a cover iff its pool graph has a $7$-clique — Lemma 4 at $n = 8$, checked
  exhaustively rather than one root per orbit, because without a center there
  is no root;
* the symmetric mode reproduces the catalogue byte for byte from $25$ of the
  $5\,568$ shards — Lemma 5 at $n = 8$.

### Theorem

> There is no fully diagonalised Latin cube of order 9.

*Proof.*  Suppose $A$ is an FDLH of order 9.  Its color classes are nine
pairwise disjoint supports (Lemma 1).  One of them, $T$, contains the
center $c_0$.  The symmetries in $G$ carry partitions to partitions and fix
the center (Lemma 3), so we may replace the whole cube by its image under any
$g \in G$, and assume $T$ is the chosen representative of its $G$-orbit.
The computation (Part II) establishes:

* $T(9)$ has exactly $14\,616\,576$ members, enumerated exhaustively;
* $11\,821\,056$ of them contain $c_0$, falling into exactly $2\,049$ $G$-orbits
  whose sizes sum to $11\,821\,056$;
* for each of the $2\,049$ orbit representatives, the companion pool is complete
  (it is a filter over the complete catalogue, re-derived independently), and
  $\omega(\Gamma(T)) \le 4$ — in fact $2$ for $2\,048$ of them, $4$ for one.

A partition containing $T$ would give a clique of size 8 in $\Gamma(T)$
(Lemma 4), so $T$ belongs to no partition.  Contradiction. ∎

A second argument, independent of the clique bound (Lemma 4), runs a depth-first exact cover over
each of the $2\,049$ pools: all $2\,049$ are exhausted, with $0$ covers found,
and no branch survives past a root plus **one** companion.  That is the depth
the exact-cover search reaches, not a bound on the size of a disjoint packing:
the search must cover every cell of the cube and kills a branch as soon as some
uncovered cell has no surviving candidate, which a mere packing is never asked
to do.  The five-support packing below is compatible with it.

Either argument alone is sufficient.  They are different algorithms run over
the same pools by the same program; `catalogue.py validate` requires both of
them, and the triangle bound, to hold before the run is called a success.

### An aside: a maximum packing containing a root

Nothing above needs a witness, but one is instructive.  The single exceptional
orbit — the one whose pool graph has clique number 4 rather than 2 — yields a
set of **five** pairwise disjoint supports of $[9]^3$, written out explicitly in
`data/exceptional_packing.json` and re-verified from the definition by
`check.py witness`.  A cube would need nine.  Any future proof of "the packing
dies early at order 9" that does not survive this example is wrong.

This is a maximum packing *containing a root*, which is all the computation
bounds.  It is not a claim about incomplete packings that avoid the center
cell: the root observation applies to partitions, not to arbitrary packings, so nothing here
says five is the global maximum.

### Lemma 5 (the plane-fixing subgroup, and the shard orbits)

This optional lemma is what `verify.sh symmetric` rests on.  It groups the
shards (the sets of supports sharing a given first face, one shard per
admissible row) into orbits under a symmetry group, such that once one shard
in an orbit has been enumerated, the contents of every other shard in that
orbit can be written down cheaply by transforming it, with no further search.
Precisely: let $H$ be the symmetries in $G$ that map the plane $x_0 = 0$ to
itself.  Each $h \in H$ sends a first-face row $p$ to another admissible row
$p^h$, and the claim is that $h$ maps the whole shard $S_p$ onto the whole
shard $S_{p^h}$ — $S_{p^h} = h(S_p)$, with equality, not just containment.
The shards therefore fall into $H$-orbits, 157 of them at order 9, and the
full catalogue is recovered from one enumerated shard per orbit by applying
every $h \in H$.  Nothing downstream changes: roots, pools and cliques are
computed from the reconstructed catalogue exactly as `full` computes them from
the searched one, and the two catalogues are compared byte for byte.  **The
theorem above does not depend on this lemma**, and neither do `verify.sh full`
and `verify.sh fast`.

Write $P = \lbrace x_0 = 0\rbrace$ for the plane whose contents define a shard: a
support's intersection with $P$ is $\lbrace(0, j, p(j))\rbrace$ for an admissible
permutation $p$ (Lemma 2), and the shard $S_p$ is the set of supports with that row 0.

> Let $H = \lbrace g \in G : g(P) = P \rbrace$ be the setwise stabiliser of $P$ in the group
> $G$ of the orbit reduction (Lemma 3).  Then
>
> (a) $H$ consists of exactly those $g = g(\pi, \varepsilon, \tau)$ with $\pi(0) = 0$ and
> $\tau(\varepsilon_0(0)) = 0$, and
>
> $$\lvert H \rvert = \lvert G \rvert / (3 \cdot 2\lfloor n/2 \rfloor) = 9216 / 24 = 384 \quad \text{at } n = 8 \text{ and } n = 9;$$
>
> (b) writing $u = \tau \circ \varepsilon_1$ and $v = \tau \circ \varepsilon_2$, an $h \in H$ acts on $P$
> by $(0,j,k) \to (0, u(j), v(k))$ if $\pi$ fixes the last two coordinates and by
> $(0,j,k) \to (0, u(k), v(j))$ if $\pi$ swaps them; so it carries the row $p$ to
> $p^h = v \circ p \circ u^{-1}$ or $v \circ p^{-1} \circ u^{-1}$, which is again admissible.  Hence the
> shard universe is closed under $H$, which permutes it;
>
> (c) for every $h \in H$ and every $p$, $h(S_p) = S_{p^h}$ — the image of a shard
> is the whole of the image shard, record by record and cell by cell.  In
> particular $\lvert S_p\rvert = \lvert S_{p^h}\rvert$.

*Proof.*  (a)  The image's first coordinate is $\tau(\varepsilon_0(x_{\pi(0)}))$.  On $P$
the coordinates $x_1, x_2$ are free, so if $\pi(0) \ne 0$ that expression takes
all $n$ values on $P$ and the image cannot be contained in the plane $x_0 = 0$;
hence $\pi(0) = 0$,
and then the image's first coordinate is the constant $\tau(\varepsilon_0(0))$, which
must be $0$.  Conversely every such $g$ maps $P$ into $P$, and an injection of a
finite set into itself is onto it.  $H$ is the stabiliser of a subset, hence a
subgroup, and $[G:H]$ is the size of the $G$-orbit of $P$.  That orbit is
$\lbrace \lbrace x_a = c\rbrace \rbrace$ with $a$ any of the three axes and $c$ any value of $\tau(0)$ or
$\tau(n-1)$: since $\tau$ permutes the pairs $\lbrace t, n-1-t\rbrace$ as blocks, $c$ ranges
over every value except the middle one of an odd $n$, so over $2\lfloor n/2 \rfloor$
values.  Hence $[G:H] = 3 \cdot 2\lfloor n/2 \rfloor = 24$ and $\lvert H\rvert = 384$.  (`symmetry group` recomputes both by construction: it selects the elements of
$G$ that fix $P$, exhibits the 24 planes, and forms all $384^2 = 147\,456$
products of pairs of selected elements, requiring each to be **an element of the
selected set**, found by lookup.  Checking instead that a product preserves $P$
would check nothing, since a composition of two plane-preserving maps preserves
the plane whatever else it does.)

(b)  With $\pi(0) = 0$, $\pi$ restricts to a permutation of $\lbrace1, 2\rbrace$, which gives
the two displayed forms.  The graph $\lbrace(j, p(j))\rbrace$ is carried to
$\lbrace(u(j), v(p(j)))\rbrace$, the graph of $v \circ p \circ u^{-1}$, or to $\lbrace(u(p(j)), v(j))\rbrace$, the
graph of $v \circ p^{-1} \circ u^{-1}$.  For admissibility, note that $u$ and $v$ are each
$\tau$ or $\tau \circ r$, and that $\tau$ commutes with $r$.  Take the first form and
put $s = u^{-1}(t)$.  Then $t$ is a fixed point of $p^h$ iff $v(p(s)) = u(s)$,
and a reflected point iff $v(p(s)) = r(u(s))$.  Cancelling $\tau$ from both
sides — legitimate because $\tau$ is a bijection commuting with $r$ — turns
those two conditions into $p(s) = s$ and $p(s) = r(s)$, in that order when
$\varepsilon_1 = \varepsilon_2$ and in the opposite order when not.  Each has exactly one
solution because $p$ is admissible, so $p^h$ has exactly one fixed and exactly
one reflected point.  The second form is the same computation after the
substitution $s = p^{-1}(u^{-1}(t))$.  So $p^h$ is admissible, and
$h \to (p \to p^h)$ is an action of $H$ on the $48\,912$ admissible permutations.
(`symmetry orbits` checks all $48\,912 \times 384 = 18\,782\,208$ images one by one —
every row against every subgroup element, not only the 157 representatives —
requiring each to be admissible and to be present in the brute-force shard
list.  It reports the number of images it checked, and stops if that is not
$\text{shards} \times \lvert H \rvert$.)

(c)  $h$ is a bijection of the cells with $h(P) = P$, so for any set $T$,
$h(T) \cap P = h(T \cap P)$.  Let $T$ be a support whose row 0 is $p$.  Then
$hT$ is a support, since $G$ maps supports to supports (Lemma 3), and its row 0 is $h(T \cap P) = p^h$ by (b).  Hence
$h(S_p) \subseteq S_{p^h}$.  Applying the same to $h^{-1}$, which is in $H$ because $H$
is a group, gives $h^{-1}(S_{p^h}) \subseteq S_p$, i.e. $S_{p^h} \subseteq h(S_p)$.  The two
are therefore equal. ∎

At $n = 9$ the $48\,912$ shards fall into **157** $H$-orbits, of sizes 384 (107
orbits), 192 (36), 96 (8), 48 (2) and 12 (4); at $n = 8$ the $5\,568$ shards fall
into 25.  By (c) the number of supports in a shard is constant on an orbit, and
that is visible in the finished catalogue: across all 157 orbits, no orbit
contains two shards with different support counts.

`verify.sh symmetric` does exactly what (c) licenses: it enumerates the 157
representatives with the same exact-cover search `full` uses, writes every other
shard as the image of its representative's payload under a recorded element of
$H$, and then checks three things — that at $n = 8$ the same construction from
25 enumerated shards reproduces `data/n8_supports.bin` set-equal *and*
byte-identical; that at $n = 9$ a uniform sample of the *mapped* shards (320 by
default, drawn with a fixed seed from the 48 755 never enumerated) agrees
byte-identically with a direct re-enumeration; and that the assembled catalogue
has the canonical SHA-256, which is a statement about the whole set of
$14\,616\,576$ supports and is where an error anywhere in the mapping would
surface.  A referee who distrusts the shard symmetry (Lemma 5) need not argue with it: `full` and
`symmetric` produce the same $1\,183\,942\,656$ bytes, and `full` never mentions
it.

---

## Part II. The computation

### The catalogue is a canonical object

The expensive artefact is `n9_supports.bin`: all $14\,616\,576$ supports of
$[9]^3$, $81$ bytes each (byte $i \cdot 9 + j$ is $k$ for the cell $`(i,j,k)`$),
$1\,183\,942\,656$ bytes in all.  It is written in a **canonical order**: shards in
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

All three modes run the order-8 controls of step 1.  From there, `full`
enumerates every shard at step 3; `symmetric` uses steps 3a-3d instead, to
enumerate the representatives and reconstruct the other shards; `fast` takes a
supplied catalogue and redoes everything from step 4 on.  Every artefact's
SHA-256 is in `checksums.txt` and is checked by `verify.sh`.  Times below were
measured by the runs recorded in this README (see *How long it takes*), not
estimated.

| # | step | program | input | expected output |
|---|------|---------|-------|-----------------|
| 1 | order-8 controls | `tests/test_n8.py` | — | 39 checks pass (below) |
| 2 | shard universe | `shards 9 list` | — | 48 912 admissible rows = A007016(9) |
| 3 | enumerate `T(9)` | `enum 9 shards` | shard list | 48 912 shards EXHAUSTED, 14 616 576 supports |
| 3a | *(`symmetric` only)* the plane-fixing subgroup | `symmetry 9 group` | — | `\|H\| = 384`, index 24, all 147 456 products elements of `H` |
| 3b | *(`symmetric` only)* shard orbits | `symmetry 9 orbits` | shard list | 157 orbits, all 48 912 x 384 = 18 782 208 images admissible and present |
| 3c | *(`symmetric` only)* enumerate, then map | `enum 9 shards`, `symmetry 9 expand` | 157 representatives | the same 48 912 shards, 14 616 576 supports |
| 3d | *(`symmetric` only)* controls | `catalogue.py sample`, `enum 9 shards`, `catalogue.py setcmpshards` | 320 mapped shards | set-equal and byte-identical to direct enumeration |
| 4 | audit, assemble | `catalogue.py audit`, `pack` | manifest, payloads, shard universe | the universe covered exactly, records in canonical order, then `n9_supports.bin` |
| 5 | definition check | `check.py verify 9` | catalogue | 14 616 576 records, 0 failures against the 301 lines |
| 6 | group closure | `orbits 9 closure` | catalogue | 87 699 456 images checked, 0 missing |
| 7 | roots | `orbits 9 center` | catalogue | 11 821 056 supports through `(4,4,4)` |
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

`n9_roots.bin` (957 505 536 B, the 11 821 056 supports through the center) is a
one-second filter over the catalogue and is not checksummed separately.
`n9_root_results.jsonl` carries no timings, so that it is reproducible byte for
byte; the raw output with wall times is left beside it as `n9_results_raw.jsonl`.

### The order-8 controls

Order 8 is where a cube exists and the census is known, so the whole pipeline is
exercised at an order where the answers were not produced by this run.
`./verify.sh test` runs 39 checks; the substantive ones are

* $T(8)$ enumerated from nothing is **13 056** supports, set-equal *and*
  byte-identical to the shipped `data/n8_supports.bin`;
* all 13 056 pass the definition-level check against the 244 lines, and the
  C and numpy constructions of the lines agree as *sets of lines*, in
  canonical form — at $n = 8$ and at $n = 9$;
* $T(8)$ is closed under the order-9216 cell group and falls into **6** orbits,
  of sizes 768, 768, 2 304, 2 304, 2 304, 4 608;
* the exact cover of $[8]^3$ by supports has exactly **198 624** solutions, with
  node counts by depth `1, 1 632, 26 016, 60 672, 69 513, 147 292, 154 660,
  198 624, 198 624`.  Relabelling the eight classes of a partition by two
  different permutations gives two different arrays, so this also counts the
  labelled cubes: $8! \cdot 198\,624 = 8\,008\,519\,680$ FDLHs of order 8;
* the **root search agrees with the census on every one of the 13 056 supports**
  — for each support, the number of covers containing it, computed once by the
  census and once by the same root machinery that is run at order 9;
* the packing ceiling and the search agree at order 8 too: a support lies in a
  cover **iff** its pool graph has a 7-clique;
* the whole of `verify.sh symmetric` is run at order 8: the plane-fixing
  subgroup has order 384 and index 24 and all $384^2$ of its products are
  elements of it, the 5 568 shards fall into **25** orbits with every one of the
  $5\,568 \times 384$ images admissible and present in the shard list, and mapping the
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
is not a coordinate in $[n]$ must be refused by all five programs that load
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
the checks this repository now runs, which add about half a minute of order-8
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

**Enumeration is the whole expense; the refutation is free.**  Building the
catalogue costs 35-43 core-hours.  Exhausting all 2 049 root cases costs
1.70 CPU-seconds, and the separate maximum-clique computation another 4.5.

**The speed-up.**  The exact part of it is a count, not a time: the 157
representatives cost **2 512 814 852** search nodes against A's measured
7.83 x 10^11 for all 48 912, a factor of **312** — the reduction factor
48 912 / 157 almost exactly.  That agreement is a measurement rather than a
consequence: the shard symmetry (Lemma 5) equates the number of supports in symmetry-related shards,
not their search effort.  In wall-clock terms, comparing within a single run —
the representatives against the per-shard cost of the 320-shard control drawn
from the same sweep, at the same moment — the three sessions give **326**,
**313** and **214**.  Counting both controls as part of the price, as they
should be since they are what makes the mode believable, the factor is **89**
to **103**.  End to end, `symmetric` took 1 426 s against `full`'s projection of
46 900 s on the same laptop in the same session: **33x**.

Two thirds of `fast` is the numpy re-derivation of every companion pool, which
is bound by memory bandwidth and page cache, and it is where most of the
run-to-run spread lives (258 s in the fastest session against 511 s here).

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
run wants roughly $3.6\,\text{GB} \times \text{pool workers} + 1.2\,\text{GB}$: about 23 GB at 6, which is
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

`--workers N` defaults to the machine's core count; `--pool-workers N` sizes the
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
`CLOCK_MONOTONIC` will not build.  Python 3.8 is the floor because of
`math.comb`, and numpy is the only Python dependency; if it lives in a
virtualenv, set `PYTHON=/path/to/python`.  Tested on macOS 15 arm64
with Apple Clang 21 and Python 3.14.5 / numpy 2.5.1, and on Linux arm64 with
GCC.

`fast` checks the supplied catalogue's SHA-256 **before** anything reads it and
refuses to continue on a mismatch, so a wrong or truncated catalogue cannot
propagate downstream; `full` and `symmetric` check the same hash on the
catalogue they have just assembled, likewise before anything reads it.

Every binary takes `--help`.  Nothing in this repository downloads anything or
sends anything anywhere.  All output goes under `--work` (default `./work`),
with two exceptions: `make` writes the six binaries into `bin/`, and the order-8
suite makes its scratch directory under `--work` when `verify.sh` runs it and
under the system temporary directory when run on its own.

---

## What a referee still has to take on trust

Part I reduces the theorem to enumerating every support, covering every root
orbit, and constructing and bounding the resulting 2 049 companion graphs —
"here are 2 049 finite graphs, count their triangles", which a reader can check
by hand.  Those computational premises remain.  Checking every stored record
establishes validity, not completeness; closure under the group does not
exclude a whole missing orbit; and a matching checksum establishes agreement
with reference bytes rather than exhaustion.

`full` establishes completeness by exhaustive search over every admissible row.
`symmetric` uses exhaustive representative searches together with the shard
symmetry (Lemma 5).
`fast` **assumes** the supplied catalogue is complete, checking its reference
digest before performing anything downstream.

In every mode the sweep is audited against the shard universe: the manifest
must cover it exactly, every record must carry its own shard's row 0, and the
records must be in strict lexicographic order within each shard and the shards
in strict row order between them.  That is an exhaustion audit rather than a
digest comparison, and it explains a discrepancy instead of merely detecting
one.

1. **That the programs implement the definitions.**  The lines are built
   twice — in C from the pattern description, and independently in numpy — and
   the two are compared as *sets of lines*, in canonical form, not by their
   counts; every record ever produced is checked against the numpy version of
   the definition.  But the definition itself is written down twice by the same
   author, and a reader who disagrees with the definition of a line will
   disagree with everything downstream.  Read `src/lines.h` and `main_lines()`
   in `src/check.py`; they are about twenty lines each.

2. **That the exhaustive enumeration is exhaustive.**  At order 8 the identical
   pipeline reproduces a census whose values — 13 056 supports, 198 624 covers,
   the node profile, the six-orbit decomposition — are supplied to the suite as
   regression targets rather than derived by it; they are the author's, quoted
   from an earlier run, so they check that the pipeline still computes what it
   computed and are not an independent authority.  And the enumerated $T(9)$ is
   *closed under the order-9216 group* (87 699 456 images checked, none
   missing), so a missing support would have had to be missed together with its
   whole orbit.

3. **Common-mode error in the model.**  Both the exact-cover enumerator and the
   definition-level checker are built on the same reading of "line" and of
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

5. **Which mode was run.**  `full` and `fast` rest on the clique bound and what leads to it (Lemmata 1-4) only.
   `symmetric` additionally rests on the shard symmetry (Lemma 5) and on `src/symmetry.c`
   implementing it: a wrong subgroup, element index or record image would mean
   48 755 of the 48 912 shards came from an untrusted map rather than a search.
   Three things stand against that, all described under Lemma 5: the n = 8
   census reproduced the same way and 320 mapped shards re-enumerated
   byte-identically, which are samples, and the canonical SHA-256 of the whole
   catalogue, which is not.  A reader who wants no lemma beyond 1-5 should run `full`.

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
src/lines.h          the lines of [n]^3, from the definition
src/enum.c           the enumerator: exact cover by dancing links
src/shards.c         the admissible row-0 permutations (A007016)
src/orbits.c         the order-9216 cell group; closure and orbit classification
src/pools.c          companion pools as a filter over a complete catalogue
src/pack.c           exhaustive exact cover, and the packing ceiling
src/symmetry.c       the plane-fixing subgroup, and the shard orbits (Lemma 5)
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
Recompute them with `./verify.sh full`, or download the catalogue and use
`./verify.sh fast`, which checks its SHA-256 against `checksums.txt` before
doing anything else:

```
curl -O https://files.apriiori.com/fdlh/n9/n9_supports.bin      # 1 183 942 656 bytes
./verify.sh fast --catalogue n9_supports.bin
```

## Provenance

The programs, the verification script and this README were written with
Claude Code (Anthropic), working from the mathematics above; the lemmas, the
checks and the checksums were reviewed by a human before being published, and
the point of the three verification modes is that nothing here need be taken
on trust from either.

## References

* W. Taylor, *On the coloration of cubes*, Discrete Mathematics **2** (1972)
  187–190.  The objects (as the property $`P(m, n)`$), the $n \ge 2^d$ corner
  bound (Proposition 4), and Problems 2 ($d = 3$, $`n = 9`$) and 3.
* E. Gergely, *A simple method for constructing doubly diagonalized Latin
  squares*, J. Combinatorial Theory Ser. A **16** (1974) 266–272.  Doubly
  diagonalized Latin squares of every order except 2 and 3; in particular
  order 12, Taylor's Problem 1.
* A. J. W. Hilton, *On double diagonal and cross Latin squares*, J. London
  Math. Soc. (2) **6** (1973) 679–689.  The earlier resolution of every order
  $m \ge 4$, credited in Taylor's added-in-proof note.
* J. Arkin, V. E. Hoggatt Jr. and E. G. Straus, *Systems of magic Latin
  k-cubes*, Canadian J. Math. **28** (1976) 1153–1161.  The name *completely
  Latin* for these objects, citing Taylor for the concept.
* D. E. Knuth, *Dancing links*, in J. Davies, B. Roscoe, J. Woodcock (eds.),
  *Millennial Perspectives in Computer Science*, Palgrave (2000), 187–214.
  Algorithm X and the doubly linked cover matrix used by `src/enum.c`.
* [OEIS A007016](https://oeis.org/A007016) — permutations with exactly one
  fixed point and exactly one reflected point; the shard count.
