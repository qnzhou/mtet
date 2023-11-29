#include <catch2/catch_test_macros.hpp>

#include <mtet/mtet.h>

TEST_CASE("basics", "[mtet]")
{
    mtet::MTetMesh mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    auto t0 = mesh.add_tet(v0, v1, v2, v3);
    REQUIRE(mesh.has_vertex(v0));
    REQUIRE(mesh.has_vertex(v1));
    REQUIRE(mesh.has_vertex(v2));
    REQUIRE(mesh.has_vertex(v3));
    REQUIRE(mesh.has_tet(t0));
    REQUIRE(mesh.get_num_vertices() == 4);
    REQUIRE(mesh.get_num_tets() == 1);

    mesh.initialize_connectivity();

    mesh.split_edge(t0, 0);
    REQUIRE(mesh.get_num_vertices() == 5);
    REQUIRE(mesh.get_num_tets() == 2);
}
