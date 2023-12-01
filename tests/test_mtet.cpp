#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <mtet/mtet.h>

#include <iostream>

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
            mesh.seq_foreach_tet([&](uint64_t tet_id, std::span<const uint64_t, 4>) {
                if (mesh.has_tet(tet_id)) {
                    mesh.split_edge(tet_id, tet_id % 4);
                    count++;
                }
            });
        }
        REQUIRE(mesh.get_num_tets() >= N);
        return count;
    };

    BENCHMARK_ADVANCED("1M split")(Catch::Benchmark::Chronometer meter)
    {
        constexpr size_t N = 1000000;
        auto mesh = generate_base_mesh();
        meter.measure([&] { return split_until(mesh, N); });
    };
}
