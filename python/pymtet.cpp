#include <mtet/grid.h>
#include <mtet/io.h>
#include <mtet/mtet.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <ankerl/unordered_dense.h>

#include <fmt/core.h>
#include <cstdint>
#include <span>

namespace nb = nanobind;

NB_MODULE(pymtet, m)
{
    using namespace nb::literals;

    nb::class_<mtet::VertexId>(
        m, "VertexId", "Unique identifier for a vertex in an :class:`MTetMesh`.")
        .def(nb::init<>(), "Construct an invalid vertex id.")
        .def_prop_ro(
            "value",
            [](mtet::VertexId& self) { return value_of(self); },
            "The underlying integer value of this id.")
        .def(nb::self == nb::self, "Check equality with another :class:`VertexId`.")
        .def("__repr__", [](mtet::VertexId& self) {
            return fmt::format("Vertex id: {}", value_of(self));
        });
    nb::class_<mtet::TetId>(
        m, "TetId", "Unique identifier for a tetrahedron in an :class:`MTetMesh`.")
        .def(nb::init<>(), "Construct an invalid tet id.")
        .def_prop_ro(
            "value",
            [](mtet::TetId& self) { return value_of(self); },
            "The underlying integer value of this id.")
        .def(nb::self == nb::self, "Check equality with another :class:`TetId`.")
        .def("__repr__", [](mtet::TetId& self) {
            return fmt::format("Tet id: {}", value_of(self));
        });
    nb::class_<mtet::EdgeId>(
        m, "EdgeId", "Unique identifier for an edge in an :class:`MTetMesh`.")
        .def(nb::init<>(), "Construct an invalid edge id.")
        .def_prop_ro(
            "value",
            [](mtet::EdgeId& self) { return value_of(self); },
            "The underlying integer value of this id.")
        .def("__repr__", [](mtet::EdgeId& self) {
            return fmt::format("Edge id: {}", value_of(self));
        });
    nb::class_<mtet::MTetMesh>(
        m,
        "MTetMesh",
        "A dynamic tetrahedral mesh that supports incremental construction and edge splits.")
        .def(nb::init<>(), "Construct an empty tetrahedral mesh.")
        .def(
            "add_vertex",
            &mtet::MTetMesh::add_vertex,
            "x"_a,
            "y"_a,
            "z"_a,
            R"(Add a vertex at the given coordinates.

:param x: X coordinate.
:param y: Y coordinate.
:param z: Z coordinate.
:return: The id of the newly added vertex.
:rtype: VertexId)")
        .def(
            "add_tet",
            &mtet::MTetMesh::add_tet,
            "v0"_a,
            "v1"_a,
            "v2"_a,
            "v3"_a,
            R"(Add a tetrahedron connecting the 4 given vertices.

:param v0: First vertex id.
:param v1: Second vertex id.
:param v2: Third vertex id.
:param v3: Fourth vertex id.
:return: The id of the newly added tetrahedron.
:rtype: TetId)")
        .def(
            "initialize_connectivity",
            &mtet::MTetMesh::initialize_connectivity,
            R"(Build tet-tet adjacency information.

Must be called once after the mesh has been constructed and before using
connectivity-dependent queries such as :meth:`get_mirror` or
:meth:`foreach_tet_around_edge`.)")
        .def(
            "has_vertex",
            &mtet::MTetMesh::has_vertex,
            "vertex_id"_a,
            R"(Check whether a vertex id refers to a vertex that exists in the mesh.

:param vertex_id: Vertex id to check.
:return: True if the vertex exists.
:rtype: bool)")
        .def(
            "has_tet",
            &mtet::MTetMesh::has_tet,
            "tet_id"_a,
            R"(Check whether a tet id refers to a tetrahedron that exists in the mesh.

:param tet_id: Tet id to check.
:return: True if the tetrahedron exists.
:rtype: bool)")
        .def(
            "has_edge",
            &mtet::MTetMesh::has_edge,
            "edge_id"_a,
            R"(Check whether an edge id refers to an edge that exists in the mesh.

:param edge_id: Edge id to check.
:return: True if the edge exists.
:rtype: bool)")
        .def(
            "get_vertex",
            [](mtet::MTetMesh& self, mtet::VertexId vertex_id) -> std::array<mtet::Scalar, 3> {
                auto v = self.get_vertex(vertex_id);
                return {v[0], v[1], v[2]};
            },
            "vertex_id"_a,
            R"(Get the coordinates of a vertex.

:param vertex_id: Id of the vertex to query.
:return: The ``[x, y, z]`` coordinates of the vertex.
:rtype: list[float]
:raises RuntimeError: If ``vertex_id`` does not exist.)")
        .def(
            "get_tet",
            [](mtet::MTetMesh& self, mtet::TetId tet_id) -> std::array<mtet::VertexId, 4> {
                auto t = self.get_tet(tet_id);
                return {t[0], t[1], t[2], t[3]};
            },
            "tet_id"_a,
            R"(Get the 4 vertices of a tetrahedron.

:param tet_id: Id of the tetrahedron to query.
:return: The 4 vertex ids ``[v0, v1, v2, v3]`` of the tetrahedron.
:rtype: list[VertexId]
:raises RuntimeError: If ``tet_id`` does not exist.)")
        .def(
            "get_edge",
            &mtet::MTetMesh::get_edge,
            "tet_id"_a,
            "local_index"_a,
            R"(Get the id of a local edge of a tetrahedron.

:param tet_id: Id of the tetrahedron.
:param local_index: Local edge index in ``[0, 6)``. See :meth:`split_edge` for the local
    edge numbering convention.
:return: The id of the requested edge.
:rtype: EdgeId
:raises RuntimeError: If ``tet_id`` does not exist or ``local_index`` is not in ``[0, 6)``.)")
        .def(
            "print",
            [](mtet::MTetMesh& self, mtet::TetId tet_id) {
                auto vts = self.get_tet(tet_id);
                fmt::print(
                    "Tet: {} {} {} {}\n",
                    value_of(vts[0]),
                    value_of(vts[1]),
                    value_of(vts[2]),
                    value_of(vts[3]));
            },
            "tet_id"_a,
            R"(Print the 4 vertex ids of a tetrahedron to stdout.

:param tet_id: Id of the tetrahedron to print.
:raises RuntimeError: If ``tet_id`` does not exist.)")
        .def(
            "get_edge_vertices",
            &mtet::MTetMesh::get_edge_vertices,
            "edge_id"_a,
            R"(Get the 2 endpoint vertices of an edge.

:param edge_id: Id of the edge to query.
:return: The 2 vertex ids of the edge's endpoints.
:rtype: list[VertexId]
:raises RuntimeError: If ``edge_id`` does not exist.)")
        .def(
            "get_edge_tet",
            &mtet::MTetMesh::get_edge_tet,
            "edge_id"_a,
            R"(Get a tetrahedron incident to the given edge.

:param edge_id: Id of the edge to query.
:return: The id of a tetrahedron containing the edge.
:rtype: TetId)")
        .def(
            "get_mirror",
            &mtet::MTetMesh::get_mirror,
            "tet_id"_a,
            "local_index"_a,
            R"(Get the tetrahedron on the other side of a local face.

Requires :meth:`initialize_connectivity` to have been called.

:param tet_id: Id of the tetrahedron.
:param local_index: Local face index of ``tet_id`` in ``[0, 4)`` whose mirror to look up.
:return: The id of the tetrahedron adjacent to ``tet_id`` across the given face.
:rtype: TetId
:raises RuntimeError: If ``tet_id`` does not exist or ``local_index`` is not in ``[0, 4)``.)")
        .def(
            "get_num_vertices",
            &mtet::MTetMesh::get_num_vertices,
            R"(Get the number of vertices currently in the mesh.

:rtype: int)")
        .def(
            "get_num_tets",
            &mtet::MTetMesh::get_num_tets,
            R"(Get the number of tetrahedra currently in the mesh.

:rtype: int)")
        .def(
            "split_edge",
            nb::overload_cast<mtet::EdgeId>(&mtet::MTetMesh::split_edge),
            "edge_id"_a,
            R"(Split an edge at its midpoint.

Inserts a new vertex at the midpoint of ``edge_id`` and splits it into two edges.

:param edge_id: Id of the edge to split.
:return: A tuple ``(vertex_id, edge_id_0, edge_id_1)``: the id of the new vertex and the
    ids of the two halves of the split edge.
:rtype: tuple[VertexId, EdgeId, EdgeId]
:raises RuntimeError: If ``edge_id`` does not exist.)")
        .def(
            "split_edge",
            nb::overload_cast<mtet::TetId, uint8_t>(&mtet::MTetMesh::split_edge),
            "tet_id"_a,
            "local_edge_id"_a,
            R"(Split a local edge of a tetrahedron at its midpoint.

Inserts a new vertex at the midpoint of the specified edge and splits it into two edges.

With the oriented tet given by ``[v0, v1, v2, v3]``, the local edges are numbered::

    0: [v0, v1]
    1: [v1, v2]
    2: [v2, v0]
    3: [v0, v3]
    4: [v1, v3]
    5: [v2, v3]

:param tet_id: Id of the tetrahedron whose edge to split.
:param local_edge_id: Local edge index in ``[0, 6)`` of the edge to split, as numbered above.
:return: A tuple ``(vertex_id, edge_id_0, edge_id_1)``: the id of the new vertex and the
    ids of the two halves of the split edge.
:rtype: tuple[VertexId, EdgeId, EdgeId]
:raises RuntimeError: If ``tet_id`` does not exist or ``local_edge_id`` is not in ``[0, 6)``.)")
        .def(
            "par_foreach_vertex",
            [](mtet::MTetMesh& self, std::function<void(mtet::VertexId)> f) {
                self.par_foreach_vertex(
                    [&](mtet::VertexId vid, std::span<const mtet::Scalar, 3>) { f(vid); });
            },
            "callback"_a,
            R"(Visit every vertex in parallel, in unspecified order.

:param callback: A callable invoked as ``callback(vertex_id)`` for each vertex. It may be
    called concurrently from multiple threads and must be thread-safe.)")
        .def(
            "seq_foreach_vertex",
            [](mtet::MTetMesh& self, std::function<void(mtet::VertexId)> f) {
                self.seq_foreach_vertex(
                    [&](mtet::VertexId vid, std::span<const mtet::Scalar, 3>) { f(vid); });
            },
            "callback"_a,
            R"(Visit every vertex sequentially.

:param callback: A callable invoked as ``callback(vertex_id)`` for each vertex, in an
    unspecified but deterministic order.)")
        .def(
            "par_foreach_tet",
            [](mtet::MTetMesh& self, std::function<void(mtet::TetId)> f) {
                self.par_foreach_tet(
                    [&](mtet::TetId tid, std::span<const mtet::VertexId, 4>) { f(tid); });
            },
            "callback"_a,
            R"(Visit every tetrahedron in parallel, in unspecified order.

:param callback: A callable invoked as ``callback(tet_id)`` for each tetrahedron. It may
    be called concurrently from multiple threads and must be thread-safe.)")
        .def(
            "seq_foreach_tet",
            [](mtet::MTetMesh& self, std::function<void(mtet::TetId)> f) {
                self.seq_foreach_tet(
                    [&](mtet::TetId tid, std::span<const mtet::VertexId, 4>) { f(tid); });
            },
            "callback"_a,
            R"(Visit every tetrahedron sequentially.

:param callback: A callable invoked as ``callback(tet_id)`` for each tetrahedron, in an
    unspecified but deterministic order.)")
        .def(
            "foreach_tet_around_edge",
            &mtet::MTetMesh::foreach_tet_around_edge,
            "edge_id"_a,
            "callback"_a,
            R"(Visit every tetrahedron incident to an edge.

Requires :meth:`initialize_connectivity` to have been called.

:param edge_id: Id of the edge whose incident tetrahedra to visit.
:param callback: A callable invoked as ``callback(tet_id)`` for each incident tetrahedron.)")
        .def(
            "export",
            [](mtet::MTetMesh& self) {
                using Vertices =
                    nb::ndarray<mtet::Scalar, nb::numpy, nb::shape<-1, 3>, nb::c_contig>;
                using Tets = nb::ndarray<int64_t, nb::numpy, nb::shape<-1, 4>, nb::c_contig>;

                size_t num_vertices = self.get_num_vertices();
                size_t num_tets = self.get_num_tets();

                using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
                IndexMap vertex_tag_map;
                vertex_tag_map.reserve(num_vertices);

                struct ExportData {
                    std::vector<mtet::Scalar> v_data;
                    std::vector<int64_t> t_data;
                };
                auto* export_data = new ExportData();
                nb::capsule owner(export_data, [](void* p) noexcept {
                    delete (ExportData*)p;
                });
                std::vector<mtet::Scalar>& v_data = export_data->v_data;
                std::vector<int64_t>& t_data = export_data->t_data;
                v_data.reserve(num_vertices * 3);
                t_data.reserve(num_tets * 4);
                self.seq_foreach_vertex(
                    [&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> data) {
                        size_t vertex_tag = vertex_tag_map.size();
                        vertex_tag_map[value_of(vid)] = vertex_tag;

                        v_data.push_back(data[0]);
                        v_data.push_back(data[1]);
                        v_data.push_back(data[2]);
                    });
                self.seq_foreach_tet([&](mtet::TetId, std::span<const mtet::VertexId, 4> data) {
                    t_data.push_back(vertex_tag_map[value_of(data[0])]);
                    t_data.push_back(vertex_tag_map[value_of(data[1])]);
                    t_data.push_back(vertex_tag_map[value_of(data[2])]);
                    t_data.push_back(vertex_tag_map[value_of(data[3])]);
                });

                Vertices vertices(v_data.data(), {num_vertices, 3}, owner);
                Tets tets(t_data.data(), {num_tets, 4}, owner);
                return std::make_tuple(std::move(vertices), std::move(tets));
            },
            R"(Export the mesh as vertex and tetrahedron arrays.

:return: A tuple ``(vertices, tets)``. ``vertices`` is an array of shape
    ``(n_vertices, 3)`` containing vertex coordinates. ``tets`` is an array
    of shape ``(n_tets, 4)`` where row ``i`` holds the 4 vertex indices
    (into ``vertices``) of tetrahedron ``i``.
:rtype: tuple[numpy.ndarray, numpy.ndarray])");

    m.def(
        "load_mesh",
        &mtet::load_mesh,
        "filename"_a,
        R"(Load a tetrahedral mesh from file.

:param filename: Path of the mesh file to load.
:return: The loaded mesh.
:rtype: MTetMesh)");
    m.def(
        "save_mesh",
        nb::overload_cast<std::string, const mtet::MTetMesh&>(&mtet::save_mesh),
        "filename"_a,
        "mesh"_a,
        R"(Save a tetrahedral mesh to file.

:param filename: Output file path.
:param mesh: Mesh to save.)");
    m.def(
        "save_mesh",
        nb::overload_cast<std::string, const mtet::MTetMesh&, const std::vector<mtet::TetId>&>(
            &mtet::save_mesh),
        "filename"_a,
        "mesh"_a,
        "active_tets"_a,
        R"(Save a subset of a tetrahedral mesh to file.

:param filename: Output file path.
:param mesh: Mesh to save.
:param active_tets: Only the tetrahedra listed here are written out; all vertices of
    ``mesh`` are included regardless.)");
    m.def(
        "save_mesh",
        [](std::string filename,
           const mtet::MTetMesh& mesh,
           std::string name,
           std::vector<mtet::Scalar>& values) {
            mtet::save_mesh(filename, mesh, name, {values.data(), values.size()});
        },
        "filename"_a,
        "mesh"_a,
        "name"_a,
        "values"_a,
        R"(Save a tetrahedral mesh together with a per-vertex scalar field.

:param filename: Output file path.
:param mesh: Mesh to save.
:param name: Name of the scalar field.
:param values: One scalar value per vertex, in the same order as
    :meth:`MTetMesh.seq_foreach_vertex` traversal.)");

    m.def(
        "generate_tet_grid",
        [](const std::array<size_t, 3>& resolution,
           const std::array<float, 3>& bbox_min,
           const std::array<float, 3>& bbox_max,
           int style) {
            mtet::GridStyle grid_style = mtet::GridStyle::TET5;
            switch (style) {
            case 5: grid_style = mtet::GridStyle::TET5; break;
            case 6: grid_style = mtet::GridStyle::TET6; break;
            default: throw std::invalid_argument("Invalid style. Use 5 or 6.");
            }
            return mtet::generate_tet_grid(resolution, bbox_min, bbox_max, grid_style);
        },
        "resolution"_a,
        "bbox_min"_a = std::array<float, 3>{0.0f, 0.0f, 0.0f},
        "bbox_max"_a = std::array<float, 3>{1.0f, 1.0f, 1.0f},
        "style"_a = 5,
        R"(Generate a tetrahedral grid with the specified resolution and bounding box.

:param resolution: Number of divisions along the x, y and z axes.
:type resolution: tuple[int, int, int]
:param bbox_min: Minimum corner of the bounding box. Defaults to ``(0, 0, 0)``.
:type bbox_min: tuple[float, float, float]
:param bbox_max: Maximum corner of the bounding box. Defaults to ``(1, 1, 1)``.
:type bbox_max: tuple[float, float, float]
:param style: Tetrahedralization style: ``5`` splits each grid cell into 5
    tets (``TET5``), ``6`` splits each cell into 6 tets (``TET6``). Defaults
    to ``5``.
:type style: int
:return: The generated tetrahedral grid mesh.
:rtype: MTetMesh
:raises ValueError: If ``style`` is not ``5`` or ``6``.)");
}
