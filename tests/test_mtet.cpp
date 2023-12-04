#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <mtet/io.h>
#include <mtet/mtet.h>

#include <algorithm>
#include <iostream>
#include <vector>

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

TEST_CASE("loops", "[mtet]")
{
    mtet::MTetMesh mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    auto v4 = mesh.add_vertex(1, 1, 1);
    auto t0 = mesh.add_tet(v0, v1, v2, v3);
    auto t2 = mesh.add_tet(v4, v3, v2, v1);
    mesh.initialize_connectivity();

    mesh.seq_foreach_vertex([&](mtet::VertexId vid, std::span<const mtet::Scalar, 3> coord) {
        REQUIRE(mesh.has_vertex(vid));
        auto v = mesh.get_vertex(vid);
        REQUIRE(v[0] == coord[0]);
        REQUIRE(v[1] == coord[1]);
        REQUIRE(v[2] == coord[2]);
    });

    mesh.seq_foreach_tet([&](mtet::TetId tid, std::span<const mtet::VertexId, 4> vertices) {
        REQUIRE(mesh.has_tet(tid));
        REQUIRE(mesh.has_vertex(vertices[0]));
        REQUIRE(mesh.has_vertex(vertices[1]));
        REQUIRE(mesh.has_vertex(vertices[2]));
        REQUIRE(mesh.has_vertex(vertices[3]));

        auto tet_vertices = mesh.get_tet(tid);
        REQUIRE(tet_vertices[0] == vertices[0]);
        REQUIRE(tet_vertices[1] == vertices[1]);
        REQUIRE(tet_vertices[2] == vertices[2]);
        REQUIRE(tet_vertices[3] == vertices[3]);

        mesh.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid, mtet::VertexId v0, mtet::VertexId v1) {
            REQUIRE(mesh.has_edge(eid));
            REQUIRE(std::find(vertices.begin(), vertices.end(), v0) != vertices.end());
            REQUIRE(std::find(vertices.begin(), vertices.end(), v1) != vertices.end());

            std::vector<mtet::TetId> edge_1_ring;
            mesh.foreach_tet_around_edge(eid, [&](mtet::TetId tid2) {
                REQUIRE(mesh.has_tet(tid2));
                edge_1_ring.push_back(tid2);
            });
            REQUIRE(std::find(edge_1_ring.begin(), edge_1_ring.end(), tid) != edge_1_ring.end());
        });
    });
}

TEST_CASE("benchmark", "[mtet][.benchmark]")
{
    auto generate_base_mesh = []() {
        mtet::MTetMesh mesh;
        auto v0 = mesh.add_vertex(0, 0, 0);
        auto v1 = mesh.add_vertex(1, 0, 0);
        auto v2 = mesh.add_vertex(0, 1, 0);
        auto v3 = mesh.add_vertex(0, 0, 1);
        auto t0 = mesh.add_tet(v0, v1, v2, v3);
        mesh.initialize_connectivity();
        REQUIRE(mesh.get_num_tets() == 1);
        return mesh;
    };

    auto split_until = [](mtet::MTetMesh& mesh, size_t N) {
        size_t count = 0;
        while (mesh.get_num_tets() < N) {
            mesh.seq_foreach_tet([&](mtet::TetId tet_id, std::span<const mtet::VertexId, 4>) {
                if (mesh.has_tet(tet_id)) {
                    mesh.split_edge(tet_id, value_of(tet_id) % 4);
                    count++;
                }
            });
        }
        REQUIRE(mesh.get_num_tets() >= N);
        return count;
    };

    constexpr size_t N = 1000000;
    auto mesh = generate_base_mesh();
    auto num_splits = split_until(mesh, N);
    save_mesh("debug.msh", mesh);
    std::cout << "num splits: " << num_splits << std::endl;
    std::cout << "  num tets: " << mesh.get_num_tets() << std::endl;

    BENCHMARK_ADVANCED("1-1M tets")(Catch::Benchmark::Chronometer meter)
    {
        auto mesh = generate_base_mesh();
        meter.measure([&] { return split_until(mesh, N); });
    };
}
