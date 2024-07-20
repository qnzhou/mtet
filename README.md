# MTet: A Mini Tetrahedral Mesh Data Structure

MTet is mini tetrahedral mesh data structure in C++. It is written by Qingnan Zhou for the Siggraph
2024 paper "[Adaptive Grid Generation for Discretizing Implicit
Complexes](https://jurwen.github.io/Adaptive-grid-for-implicit-complexes/)".

## Features

* Tet key with mirror indices for efficient storage of local vertex-vertex correspondence between adjacent tets.
* O(1) worst-case vertex/tet creation and access based on [slot map](https://github.com/SergeyMakeev/slot_map).
* Efficient edge split support.
* Strong index [type safety](https://github.com/rollbear/strong_type).

## Quick start

To create tet mesh with a single tet:

```cpp
#include <mtet/mtet.h>

mtet::MTetMesh mesh;
auto v0 = mesh.add_vertex(0, 0, 0);
auto v1 = mesh.add_vertex(1, 0, 0);
auto v2 = mesh.add_vertex(0, 1, 0);
auto v3 = mesh.add_vertex(0, 0, 1);
auto tid = mesh.add_tet(v0, v1, v2, v3);

assert(mesh.get_num_vertices() == 4);
assert(mesh.get_num_tets() == 1);
```

To access vertices and tets:
```cpp
mesh.has_vertex(v0);
auto p0 = mesh.get_vertex(v0);

mesh.has_tet(tid);
auto t = mesh.get_tet(tid);

assert(t[0] == v0);
assert(t[1] == v1);
assert(t[2] == v2);
assert(t[3] == v3);
```

To initialize and query tet connectivity:
```cpp
mesh.initialize_connectivity();

auto adj_tid = mesh.get_mirror(tid, 0); // Adjacent tet across facet 0.
```

To split an edge:
```cpp
auto e = mesh.get_edge(tid, 0); // Edge 0 of tet tid.

auto [vid, e0, e1] = mesh.split_edge(e);
```

One can also batch process vertices and tets using callback functions:

```cpp
auto vfn = [](VertexId vid, std::span<const Scalar, 3> position) {
    // Do something with the vertex.
};

mesh.par_foreach_vertex(vfn); // pararllel vertex callback.
mesh.seq_foreach_vertex(vfn); // sequential vertex callback.

auto tfn = [](TetId tid, std::span<const VertexId, 4> vts) {
    // Do something with the tet.
};

mesh.par_foreach_tet(tfn); // pararllel tet callback.
mesh.seq_foreach_tet(tfn); // sequential tet callback.
```

Similarly, once can use callback functions sequentially for each edges within a tet and for each tet
around an edge:

```cpp
auto efn = [](EdgeId eid, VertexId v0, VertexId v1) {
    // Do something with the edge.
};
mesh.foreach_edge_in_tet(tid, efn);

auto tfn = [](TetId tid) {
    // Do something with the tet.
};
mesh.foreach_tet_around_edge(eid, tfn);
```

## Performance

We use 1-to-1,000,000 tet mesh generation via sequential edge split as a benchmark.
```sh
$ ./mtet_tests benchmark
```

The following are the benchmark results measured on my MacBook Pro (2.4 GHz 8-Core Intel Core i9).
```
benchmark name                       samples       iterations    estimated
                                     mean          low mean      high mean
                                     std dev       low std dev   high std dev
-------------------------------------------------------------------------------
1-1M tets                                      100             1     32.5669 s
                                        303.581 ms     300.74 ms     307.49 ms
                                        16.8615 ms     12.811 ms    25.9629 ms
```


## Technical details

The following contains some technical details about the data structure design. It serves more like a
note for developers, and it is not necessary to understand these details to use the library.

### Mirror index

Mirror index is a small data structure that stores both tet-tet adjacency as well as the local
vertex-vertex "mirror correspondence" between the two facet-adjacent tets, all in a 64-bit integer. It is
motivated by the need of local mesh operations (e.g. edge collapse, split) which requires adjacent
elements to be updated jointly.


<img src=https://github.com/user-attachments/assets/573112c4-c9e2-4619-9cbc-97b5e98a375e width=50% />

In the above example, two adjacent tets are shown. The left tet has vertices `v0`, `v1`, `v2` and
`v3`. The right tet has vertices `u0`, `u1`, `u2` and `u3`. The mirror correspondence of local vertices between
these two tets are the following:

| v | u |
|---|---|
| 0 | 0 |
| 1 | 3 |
| 2 | 1 |
| 3 | 2 |

With mirror correspondence, we can quickly identify the local index of the shared vertex/edge/facet
between adjacent tets. For example, edge (`v2`, `v3`) in the first tet is the same as edge (`u1`,
`u2`) in the second tet. Since the valid range of a local vertex index is from 0 to 3, it can be
stored in 2 bits, and we can pack all four local vertex indices into the tags of a 64-bit tet slot
map key.

Tet key:
```
+--------------------------------+--------------------+----+--+--+--+--+
|        index (32 bits)         |  version (20 bits) | ei |m0|m1|m2|m3|
+--------------------------------+--------------------+----+--+--+--+--+
```

* `index` (32 bits): The index of the adjacent tet.
* `version` (20 bits): The version of the adjacent tet. It is used to detect whether the adjacent tet has been
  deleted/replaced.
* `ei` (4 bits): The local edge index. (valid range: 0-5)
* `m0`, `m1`, `m2`, `m3` (2 bits): The local vertex indices of the adjacent tet that mirrors the local
  vertices of the current tet. (valid range: 0-3)


We call such tet key with local mirror correspondence a "mirror index". Since each tet can have at
most 4 facet-adjacent tets, we only need 4 mirror indices to store the complete correspondence with
all facet-adjacent tets. If a tet is on the boundary, one or more of its facets are boundary facets.
We use a special invalid tet key as mirror index to indicate the corresponding tet facet is a
boundary facet. A total of 4x64=256 bits are needed to store all 4 mirror indices, which is much
smaller than half-edge/half-face data structure.

### Tet Id and Edge Id

Note that in addition to using tet key as mirror index. We can also use tet key to just specify a
particular tet or an edge of that tet depending on the context. We use [strong
type](https://github.com/rollbear/strong_type) to distinguish the different uses case. The type
`TetId` is used to specify a tet or a mirror index, `EdgeId` is used to specify an edge.

When specifying an edge using `EdgeId`, the index part of the tet key yields a tet containing the
target edge, and the local edge index part of the tet key specifies the local edge index within that
tet. The local edge index convention is illustrated bellow:

```
                   v2
                 ,/|`\
               ,/  |  `\
             ,2    |.   `1
           ,/       5     `\
         ,/         |       `\
       v0-------0-- |. -------v1
         `\.         |      ,/
            `\.      |    ,4
               `3.   |. ,/
                  `\. |/
                     `v3
```
