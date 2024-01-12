#include <mtet/io.h>
#include <mtet/mtet.h>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <fmt/core.h>
#include <span>

namespace nb = nanobind;

NB_MODULE(pymtet, m)
{
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
        .def("foreach_tet_around_edge", &mtet::MTetMesh::foreach_tet_around_edge);

    m.def("load_mesh", &mtet::load_mesh);
    m.def("save_mesh", nb::overload_cast<std::string, const mtet::MTetMesh&>(&mtet::save_mesh));
    m.def(
        "save_mesh",
        nb::overload_cast<std::string, const mtet::MTetMesh&, std::span<mtet::TetId>>(
            &mtet::save_mesh));
}
