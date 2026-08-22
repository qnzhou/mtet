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

    nb::class_<mtet::VertexId>(m, "VertexId")
        .def(nb::init<>())
        .def_prop_ro("value", [](mtet::VertexId& self) { return value_of(self); })
        .def(nb::self == nb::self)
        .def("__repr__", [](mtet::VertexId& self) {
            return fmt::format("Vertex id: {}", value_of(self));
        });
    nb::class_<mtet::TetId>(m, "TetId")
        .def(nb::init<>())
        .def_prop_ro("value", [](mtet::TetId& self) { return value_of(self); })
        .def(nb::self == nb::self)
        .def("__repr__", [](mtet::TetId& self) {
            return fmt::format("Tet id: {}", value_of(self));
        });
    nb::class_<mtet::EdgeId>(m, "EdgeId")
        .def(nb::init<>())
        .def_prop_ro("value", [](mtet::EdgeId& self) { return value_of(self); })
        .def("__repr__", [](mtet::EdgeId& self) {
            return fmt::format("Edge id: {}", value_of(self));
        });
    nb::class_<mtet::MTetMesh>(m, "MTetMesh")
        .def(nb::init<>())
        .def("add_vertex", &mtet::MTetMesh::add_vertex)
        .def("add_tet", &mtet::MTetMesh::add_tet)
        .def("initialize_connectivity", &mtet::MTetMesh::initialize_connectivity)
        .def("has_vertex", &mtet::MTetMesh::has_vertex)
        .def("has_tet", &mtet::MTetMesh::has_tet)
        .def("has_edge", &mtet::MTetMesh::has_edge)
        .def(
            "get_vertex",
            [](mtet::MTetMesh& self, mtet::VertexId vertex_id) -> std::array<mtet::Scalar, 3> {
                auto v = self.get_vertex(vertex_id);
                return {v[0], v[1], v[2]};
            })
        .def(
            "get_tet",
            [](mtet::MTetMesh& self, mtet::TetId tet_id) -> std::array<mtet::VertexId, 4> {
                auto t = self.get_tet(tet_id);
                return {t[0], t[1], t[2], t[3]};
            })
        .def("get_edge", &mtet::MTetMesh::get_edge)
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
            })
        .def("get_edge_vertices", &mtet::MTetMesh::get_edge_vertices)
        .def("get_edge_tet", &mtet::MTetMesh::get_edge_tet)
        .def("get_mirror", &mtet::MTetMesh::get_mirror)
        .def("get_num_vertices", &mtet::MTetMesh::get_num_vertices)
        .def("get_num_tets", &mtet::MTetMesh::get_num_tets)
        .def("split_edge", nb::overload_cast<mtet::EdgeId>(&mtet::MTetMesh::split_edge))
        .def("split_edge", nb::overload_cast<mtet::TetId, uint8_t>(&mtet::MTetMesh::split_edge))
        .def(
            "par_foreach_vertex",
            [](mtet::MTetMesh& self, std::function<void(mtet::VertexId)> f) {
                self.par_foreach_vertex(
                    [&](mtet::VertexId vid, std::span<const mtet::Scalar, 3>) { f(vid); });
            })
        .def(
            "seq_foreach_vertex",
            [](mtet::MTetMesh& self, std::function<void(mtet::VertexId)> f) {
                self.seq_foreach_vertex(
                    [&](mtet::VertexId vid, std::span<const mtet::Scalar, 3>) { f(vid); });
            })
        .def(
            "par_foreach_tet",
            [](mtet::MTetMesh& self, std::function<void(mtet::TetId)> f) {
                self.par_foreach_tet(
                    [&](mtet::TetId tid, std::span<const mtet::VertexId, 4>) { f(tid); });
            })
        .def(
            "seq_foreach_tet",
            [](mtet::MTetMesh& self, std::function<void(mtet::TetId)> f) {
                self.seq_foreach_tet(
                    [&](mtet::TetId tid, std::span<const mtet::VertexId, 4>) { f(tid); });
            })
        .def("foreach_tet_around_edge", &mtet::MTetMesh::foreach_tet_around_edge)
        .def("export", [](mtet::MTetMesh& self) {
            using Vertices = nb::ndarray<mtet::Scalar, nb::numpy, nb::shape<-1, 3>, nb::c_contig>;
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
            nb::capsule owner(export_data, [](void* p) noexcept { delete (ExportData*)p; });
            std::vector<mtet::Scalar>& v_data = export_data->v_data;
            std::vector<int64_t>& t_data = export_data->t_data;
            v_data.reserve(num_vertices * 3);
            t_data.reserve(num_tets * 4);
            self.seq_foreach_vertex([&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> data) {
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
        });

    m.def("load_mesh", &mtet::load_mesh);
    m.def("save_mesh", nb::overload_cast<std::string, const mtet::MTetMesh&>(&mtet::save_mesh));
    m.def(
        "save_mesh",
        nb::overload_cast<std::string, const mtet::MTetMesh&, std::span<mtet::TetId>>(
            &mtet::save_mesh));
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
        "values"_a);

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

@param resolution The number of divisions along each axis (x, y, z).
@param bbox_min The minimum coordinates of the bounding box.
@param bbox_max The maximum coordinates of the bounding box.
@param style The style of the tetrahedral mesh (5 or 6).

@return A tetrahedral mesh object.)");
}
