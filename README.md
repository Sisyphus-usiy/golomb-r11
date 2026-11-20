# The shortest real 11-mark Golomb ruler has length 72

**Theorem.** Suppose the distances between different pairs of marks are at least one unit apart, and each distance is at least one. Then every real 11-mark ruler has length at least 72. Equality is attained by

$$
(0,1,4,13,28,33,47,54,64,70,72).
$$

The lower bound is a computer-assisted proof. The independent checker verifies the upper witness, 7,507,832 lower-bound certificates, their coverage, and all smaller bounds used in the proof.

## Verify or regenerate

Two files are too large for a Git repository on GitHub, so they ship as assets of the [`v1` release](https://github.com/Sisyphus-usiy/golomb-r11/releases/tag/v1):

- `certificates/n11.rgc.gz` (148 MiB, above the 100 MiB per-file limit). Download it into `certificates/` before running the full check.
- `real-golomb-11.bundle` (148 MiB), a self-contained Git bundle of this repository including that certificate. It clones offline with `git clone real-golomb-11.bundle golomb-r11`, which needs no further download.

From the directory containing this file, run the following with Python 3. The checker uses only the standard library and arbitrary-precision integer arithmetic.

```sh
python -I check.py --through 10
python -I check.py
```

The second command also verifies the smaller proofs in order. It exits successfully only after checking `n11.rgc.gz`. To build fresh certificates from source, supply a C++17 compiler and run

```sh
python -I reproduce.py --through 11 --out regenerated
```

The driver defaults to `g++` and accepts `--cxx clang++`. It builds `search.cpp` in a temporary directory, generates fresh proofs through ten marks, runs the independent checker, and only then starts the 11-mark search. All proof data needed for verification are included under `certificates/` once `n11.rgc.gz` has been fetched from the release.

| Marks | Target | Tree nodes | Certified leaves | Maximum depth | Search | Check including smaller proofs |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 55 | 20,635 | 10,318 | 28 | 0.033 s | 0.264 s |
| 11 | 72 | 15,015,663 | 7,507,832 | 59 | 52.843 s | 108.381 s |

These wall times were measured on an Intel Core i9-13900HX with GCC 16.1 and CPython 3.14.7, one process per stage. Search times exclude compilation and gzip compression.

Re-checked on 2026-11-20 on a second machine (Windows, CPython 3.14.7): `python -I check.py` printed `PASS through n=11 total_ms=125149`, i.e. 125.1 s, with `nodes=15015663 leaves=7507832 empty=0 max_depth=59` for the eleven-mark tree and exit status 0.

## Proof of coverage and of the lower bound

Use marks numbered $0,\ldots,n-1$. Subtracting $x_0$ preserves every distance and puts the first mark at zero. Write $g_k=x_{k+1}-x_k$, $L=\sum_k g_k$, and $d_{ij}=x_j-x_i$. Each $d_{ij}$ is a sum of consecutive gaps. Reflection reverses the gap vector and preserves $L$. The two end gaps differ by at least one, so exactly one orientation has $g_{n-2}-g_0\ge1$.

A complete ordering of the distances gives a polyhedral chamber. The tree describes these chambers by partial orderings. At a split comparing distinct distances $u,v$, its two children impose $u-v\ge1$ and $v-u\ge1$. Every admissible ruler follows exactly one child. In a reflected chamber the same proof applies with the gap coordinates reversed.

The following inequalities hold wherever the checker admits them.

1. Each proper interval spanning $s$ marks contains a translated $s$-mark ruler. Its distance is at least the lower bound already proved for $s$. For $s=2,3,4$ the counting bound $\binom{s}{2}$ suffices. In fact the $j$-th smallest distance is at least $j$, which also gives $L\ge\binom n2$.
2. For any $r$ distinct distances, their sum is at least $r(r+1)/2$ by the same ranking argument. For proper distances $d$, the numbers $L-d$ are also positive, at least one, and mutually unit-separated. Any $r$ of these satisfy the same sum bound.
3. The reflection inequality holds at the root. Each branch inequality holds throughout its own subtree.

At every length leaf the checker reconstructs rows $a_i\cdot g\ge b_i$ and checks positive integer weights $w_i$ and a positive denominator $q$ satisfying

$$
\sum_i w_i a_i=q\mathbf 1,\qquad \sum_i w_i b_i\ge qT.
$$

Adding the inequalities gives $qL\ge qT$ for **every** ruler in every chamber extending that path. The checker reads both children at each split and requires end-of-file immediately after the root tree. Induction over that finite tree proves the lower bound for the canonical orientation, and reflection proves it for both. Induction in $n$ supplies the sub-ruler bounds $11,17,25,34,44,55$ before the checker accepts the 11-mark bound. Finally it computes all 55 differences of the displayed witness and checks unit separation. Hence $G_{\mathbb R}(11)=72$.

## A tight leaf by hand

Consider these two sets of indexed distances.

$$
\begin{aligned}
A&=\{d_{0,1},d_{0,2},d_{1,2},d_{9,10}\},\\
B&=A\cup\{d_{2,3},d_{4,5},d_{6,7},d_{7,8},d_{8,9},d_{8,10}\}.
\end{aligned}
$$

At the leaf containing the displayed ruler, the five parenthesized differences below are all at least one. Expanding $d_{ij}=x_j-x_i$ gives the identity

$$
\begin{aligned}
L={}&\sum_{u\in A}u+\sum_{v\in B}v
 +2(d_{2,3}-d_{8,10})+(d_{8,9}-d_{4,5})\\
 &+2(d_{5,6}-d_{0,3})+(d_{3,4}-d_{5,6})+(d_{4,5}-d_{0,2}).
\end{aligned}
$$

Thus this leaf proves $L\ge10+55+(2+1+2+1+1)=72$. On the witness, $A$ takes the values $1,\ldots,4$, $B$ takes the values $1,\ldots,10$, and each comparison is tight.

## Search and certificate format

The generator minimizes $L$ over a relaxation at each node. It starts with a dual-feasible unit-gap basis. For a basis matrix $H$ it keeps $|\det H|H^{-1}$ and scales primal coordinates and dual multipliers by $|\det H|$. Exact simplex pivots preserve dual feasibility. The implementation checks integer ranges and exact divisions. An available dual bound of at least $T$ closes the node. Otherwise the generator separates packing inequalities at an optimal LP point by sorting the distances and checking prefix sums in both directions. It branches on a violated comparison. Overlapping intervals can be cancelled to disjoint intervals, and proper containment already gives separation. The checker independently reconstructs the entire proof from the records.

Each `.rgc.gz` file is a gzip-compressed ASCII stream of integer records. Intervals $(i,j)$, $0\le i<j<n$, are numbered from zero in lexicographic order. Write $m=\binom n2$.

| Record | Meaning |
|---|---|
| `RGC1 n T`, then `W x0 ...` | Expected claim and integer upper-bound witness. |
| `P mask` or `C mask` | Declare a packing row using the intervals whose bits are set. `P` uses their distances. `C` uses their complements in $L$ and excludes the full interval. The checker derives every coefficient and the right side. |
| `B u v` | Binary split in preorder. The next subtree assumes $d_u-d_v\ge1$. The following subtree assumes the reverse. |
| `L q ref w ...` | Weighted-row identity above. For `q=0` the checker instead requires zero coefficients and strictly positive right side, certifying an empty node. All leaves in the delivered files have $q>0$. |

Row references $0,\ldots,m-1$ select interval bounds, $m$ selects reflection, and $m+1,m+2,\ldots$ select packing declarations in stream order. Reference $-1-h$ selects the current path's branch at depth $h$ starting from zero. Packing declarations consume no tree nodes. Every leaf carries its own multipliers.

The 11-mark proof has 562,539,124 uncompressed bytes. SHA-256 of the uncompressed record stream is `5e6acc46ba7205190f7079be1528c7502bbf86502ac8cb15df1d11f828c1798d`.

## Prior work and date check

Esther Notik proves the result through ten marks. Sections 2 and 3 give the chamber argument, packing cuts, reflection and induction on smaller rulers. The programs here implement that argument independently, with integer LP pivots and complementary packing cuts. [Notik, arXiv:2611.22656](https://arxiv.org/pdf/2611.22656)

For November 16–20, 2026, the arXiv math.CO submission-date query returned 176 records at the time of the search. A title-and-abstract scan for Golomb, ruler, Sidon, separated or distinct differences, and difference triangles returned only Notik's paper. A separate search of the math.CO listings also found no paper settling eleven marks. [Date-bounded arXiv query](https://export.arxiv.org/api/query?search_query=cat%3Amath.CO%20AND%20submittedDate%3A%5B202611160000%20TO%20202611202359%5D&max_results=500&sortBy=submittedDate&sortOrder=ascending) · [November math.CO index](https://arxiv.org/list/math.CO/2026-11?show=2000&skip=0)
