#include <catch2/catch_test_macros.hpp>

#include <MTetMeshImpl.h>

#include <array>

TEST_CASE("TriangleEqual", "[unordered_dense]")
{
    std::array<uint64_t, 3> t0 = {{0, 1, 2}};
    std::array<uint64_t, 3> t1 = {{2, 0, 1}};
    std::array<uint64_t, 3> t2 = {{2, 1, 0}};
    std::array<uint64_t, 3> t3 = {{1, 1, 1}};

    mtet::TriangleEqual eq;

    REQUIRE(eq(t0, t0));
    REQUIRE(eq(t0, t1));
    REQUIRE(eq(t1, t1));
    REQUIRE(eq(t2, t0));
    REQUIRE(eq(t3, t3));
    REQUIRE(!eq(t3, t0));
    REQUIRE(!eq(t3, t2));
}

TEST_CASE("TriangleHash", "[unordered_dense]")
{
    std::array<uint64_t, 3> t0 = {{0, 1, 2}};
    std::array<uint64_t, 3> t1 = {{2, 0, 1}};
    std::array<uint64_t, 3> t2 = {{2, 1, 0}};
    std::array<uint64_t, 3> t3 = {{1, 1, 1}};

    mtet::TriangleHash hash;

    REQUIRE(hash(t0) == hash(t0));
    REQUIRE(hash(t0) == hash(t1));
    REQUIRE(hash(t1) == hash(t1));
    REQUIRE(hash(t2) == hash(t0));
    REQUIRE(hash(t3) == hash(t3));
    REQUIRE(hash(t3) != hash(t0));
    REQUIRE(hash(t3) != hash(t2));
}

TEST_CASE("invalid_key", "[slotmap]")
{
    REQUIRE(mtet::invalid_key == mtet::MTetMeshImpl::VertexKey::invalid());
    REQUIRE(mtet::invalid_key == mtet::MTetMeshImpl::TetKey::invalid());
}

TEST_CASE("tag", "[slotmap]")
{
    mtet::MTetMeshImpl::TetMap tet_map;
    auto key = tet_map.emplace(mtet::MTet());
    auto key_copy = key;
    REQUIRE(key == key_copy);

    for (uint8_t i = 0; i < 4; i++) {
        mtet::set_tag(key, i, i);
        REQUIRE(mtet::get_tag(key, i) == i);
    }
    REQUIRE(key != key_copy);
}

TEST_CASE("unordered_dense", "[unordered_dense]")
{
    using HashMap = ankerl::unordered_dense::map<
        std::array<uint64_t, 3>,
        std::array<uint64_t, 2>,
        mtet::TriangleHash,
        mtet::TriangleEqual>;

    std::array<uint64_t, 3> t0 = {{0, 1, 2}};
    std::array<uint64_t, 3> t1 = {{2, 0, 1}};
    std::array<uint64_t, 3> t2 = {{2, 1, 0}};
    std::array<uint64_t, 3> t3 = {{1, 1, 1}};

    HashMap triangle_map;
    triangle_map[t0] = {{0, 1}};
    REQUIRE(triangle_map.contains(t0));
    REQUIRE(triangle_map.contains(t1));
    REQUIRE(triangle_map.contains(t2));
    REQUIRE(!triangle_map.contains(t3));
    REQUIRE(triangle_map[t0] == std::array<uint64_t, 2>({0, 1}));
}

TEST_CASE("adjacency", "[mtet]")
{
    mtet::MTetMeshImpl mesh;
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

    REQUIRE(mesh.get_adjacent_tet(t0, 0) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t0, 1) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t0, 2) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t0, 3) == mtet::invalid_key);

    // Add second tet.
    auto v4 = mesh.add_vertex(1, 1, 1);
    auto t1 = mesh.add_tet(v4, v3, v2, v1);
    REQUIRE(mesh.has_vertex(v4));
    REQUIRE(mesh.has_tet(t1));
    mesh.initialize_connectivity();

    using TetKey = mtet::MTetMeshImpl::TetKey;
    REQUIRE(mesh.get_adjacent_tet(t0, 1) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t0, 2) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t0, 3) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t1, 1) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t1, 2) == mtet::invalid_key);
    REQUIRE(mesh.get_adjacent_tet(t1, 3) == mtet::invalid_key);

    // The opposite tet of t0 across face 0 is t1.
    TetKey t0_0(mesh.get_adjacent_tet(t0, 0));
    TetKey t1_0(mesh.get_adjacent_tet(t1, 0));

    REQUIRE(TetKey::toIndex(t0_0) == TetKey::toIndex(TetKey(t1)));
    REQUIRE(TetKey::toIndex(t1_0) == TetKey::toIndex(TetKey(t0)));

    REQUIRE(mtet::get_tag(t0_0, 0) == 0);
    REQUIRE(mtet::get_tag(t0_0, 1) == 3);
    REQUIRE(mtet::get_tag(t0_0, 2) == 2);
    REQUIRE(mtet::get_tag(t0_0, 3) == 1);

    REQUIRE(mtet::get_tag(t1_0, 0) == 0);
    REQUIRE(mtet::get_tag(t1_0, 1) == 3);
    REQUIRE(mtet::get_tag(t1_0, 2) == 2);
    REQUIRE(mtet::get_tag(t1_0, 3) == 1);
}
