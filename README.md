# genll

`genll` generates all **locally linear graphs** on a given number of vertices, one
from each isomorphism class. A graph is locally linear when every edge lies in exactly
one triangle. Equivalently, the neighbourhood of every vertex induces a perfect
matching.

Its interface follows `geng` from [nauty](https://pallini.di.uniroma1.it/): the order
is given on the command line, graphs are written to standard output in graph6 format,
and a run can be split into independent parts with `res/mod`. Options restrict the
output to C4-free, C5-free or maximal graphs, or to a minimum degree. These
restrictions are built into the search rather than applied to its output afterwards.

## Building

You need:

* a C++20 compiler: GCC or Clang (tested with GCC 13 and Clang 18);
* nauty, built from source (tested with nauty 2.9.1). genll links against nauty's
  archives `nautyW1.a` and `nautyL1.a`, which nauty's normal build produces:

  ```sh
  cd /path/to/nauty
  ./configure && make
  ```

Then, in this directory:

```sh
make NAUTY=/path/to/nauty           # builds genll and genllL
make check NAUTY=/path/to/nauty     # optional: compares against known counts (~15 s)
```

`NAUTY` defaults to `../nauty`. The build produces two executables from the same
source:

| executable | orders     | nauty archive |
|------------|------------|---------------|
| `genll`    | N ≤ 32     | `nautyW1.a`   |
| `genllL`   | N ≤ 64     | `nautyL1.a`   |

The default flags include `-march=native`, so the binaries are tuned to the build
machine. Set `OPT` to override them, e.g. `make OPT=-O3`. Three compile-time switches
are meant for development:

* `STATS=1` prints counts of where candidates are accepted and rejected;
* `ROUNDS=n` sets the number of colour-refinement rounds in the vertex invariant
  (default 3);
* `VERIFY=1` checks the incrementally maintained invariants against a direct
  computation (slow).

## Usage

```
genll [-u] [-q] [-f] [-5] [-M] [-d#] [res/mod] N
```

| option    | meaning |
|-----------|---------|
| `-f`      | only C4-free graphs |
| `-5`      | only C5-free graphs; these are automatically C4-free (see below); cannot be combined with `-M` |
| `-M`      | only maximal graphs: no triangle can be added while keeping the graph locally linear (with `-f`: locally linear and C4-free) |
| `-d#`     | minimum degree at least `#` (rounded up to an even number, since every degree is even) |
| `-u`      | count only, write no graphs |
| `-q`      | no summary on standard error |
| `res/mod` | generate only part `res` of `mod` (`0 ≤ res < mod`) |

The graphs go to standard output, one graph6 string per line. A summary goes to
standard error, as with `geng`:

```
$ genll 7
>A genll n=7
F????
Fw???
FwCW?
F{c??
F{eCG
F{dAG
>Z 6 graphs generated in 6.6783e-05 sec
```

All locally linear graphs are generated, including disconnected ones and ones with
isolated vertices; the edgeless graph is locally linear.

### Splitting a run

The parts `0/mod`, `1/mod`, …, `(mod-1)/mod` are disjoint and together give exactly
the output of the unsplit run. All parts must be given the same other options. For
example, to use eight cores:

```sh
for r in 0 1 2 3 4 5 6 7; do genll $r/8 20 > ll20.$r.g6 & done; wait
```

### C5-free graphs

A locally linear graph that contains a 4-cycle also contains a 5-cycle. The triangle on
any edge of the 4-cycle has its third vertex outside the 4-cycle (otherwise an edge
would lie in two triangles). So `-5` and `-f -5` generate the same graphs.

## Method

genll uses the canonical construction path method: each graph is built from a smaller
one by one augmentation and kept only if that augmentation is the canonical one.
Details are in the comments in `genll.cpp`.

* **Augmentation.** A new vertex `u`, the *hub*, is added together with `k ≥ 0`
  disjoint pairs `{v_i, w_i}` of existing vertices. The new edges are `u v_i`,
  `u w_i` and `v_i w_i`. The result is locally linear exactly when every pair is at
  distance at least 3 and no edge joins two chosen vertices. For C4-free graphs the
  conditions become distance at least 4 within a pair and at least 3 between pairs.
  Local linearity therefore never has to be tested.
* **Reduction.** The inverse operation deletes a vertex together with the perfect
  matching on its neighbourhood. Every locally linear graph can be reduced this way at
  any vertex, so every graph is reached.
* **Canonicity.** The canonical vertex to delete has minimum degree. Ties are broken
  by a ladder of cheap vertex invariants, and remaining ties by a canonical labelling
  from nauty. nauty is called for a labelling only when the new vertex ties with
  another one on all of the invariants.
* **Automorphisms.** Of each orbit of pair sets under the automorphism group of the
  parent, only the lexicographically least one is used.
* **Restrictions.** `-f` tightens the conditions on the pairs. `-5` rejects a new
  graph as soon as it has a 5-cycle. `-d#` cuts off every branch that can no longer
  reach the minimum degree. `-M` restricts the last augmentation to those that leave
  no room for another triangle.


## Known counts

These are the numbers of graphs on `n` vertices, up to isomorphism, that genll
reproduces. The file `counts.txt` holds the counts that `make check` verifies,
including those of the C5-free graphs.

|  n | locally linear | maximal | C4-free | maximal C4-free |
|---:|---:|---:|---:|---:|
|  5 | 3 | 1 | 3 | 1 |
|  6 | 4 | 2 | 4 | 2 |
|  7 | 6 | 2 | 6 | 2 |
|  8 | 8 | 3 | 7 | 3 |
|  9 | 15 | 4 | 12 | 4 |
| 10 | 22 | 7 | 16 | 7 |
| 11 | 41 | 9 | 25 | 9 |
| 12 | 80 | 16 | 35 | 12 |
| 13 | 201 | 33 | 60 | 20 |
| 14 | 584 | 83 | 91 | 30 |
| 15 | 2,498 | 281 | 169 | 53 |
| 16 | 14,202 | 1,311 | 288 | 80 |
| 17 | 115,152 | 8,527 | 598 | 160 |
| 18 | 1,235,365 | 74,070 | 1,248 | 284 |
| 19 | 17,150,323 | 799,999 | 3,153 | 657 |
| 20 | 299,882,532 | 10,661,060 | 8,576 | 1,469 |

## Files

| file | contents |
|------|----------|
| `genll.cpp`      | the generator (a single translation unit) |
| `nautyBitset.h`  | one-word bitset with nauty's bit order, so the adjacency rows are passed to nauty without conversion |
| `graph.h`        | graph as an array of bitset rows, and the graph6 writer |
| `Makefile`       | build and test |
| `check.sh`, `counts.txt` | the test run by `make check` |
| `LICENSE`        | the GNU General Public License, version 3 |

## Author

Steven Van Overberghe, Ghent University.

The implementation was written with the assistance of Claude (Anthropic), following
the author's design. Its output has been checked against independently generated
reference data.

## License

Copyright (C) 2023–2026 Steven Van Overberghe.

genll is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. See `LICENSE`.

genll links against nauty, which is distributed separately under the Apache License 2.0.
