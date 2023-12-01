#include <catch2/catch_test_macros.hpp>

#include <MTetMeshImpl.h>

#include <array>

void validate_mesh(const mtet::MTetMeshImpl& mesh)
{
    mesh.par_foreach_tet([&](uint64_t tet_id, std::span<const uint64_t, 4> tet_vertices) {
        using VertexKey = mtet::MTetMeshImpl::VertexKey;
        using TetKey = mtet::MTetMeshImpl::TetKey;
        TetKey key(tet_id);

        const auto& vertices = mesh.get_vertices();
        const auto& tets = mesh.get_tets();
        REQUIRE(tets.has_key(key));

        auto ptr = tets.get(key);
        REQUIRE(ptr != nullptr);
        const auto& tet = *ptr;

        REQUIRE(tet.vertices[0] == tet_vertices[0]);
        REQUIRE(tet.vertices[1] == tet_vertices[1]);
        REQUIRE(tet.vertices[2] == tet_vertices[2]);
        REQUIRE(tet.vertices[3] == tet_vertices[3]);

        REQUIRE(vertices.has_key(VertexKey(tet.vertices[0])));
        REQUIRE(vertices.has_key(VertexKey(tet.vertices[1])));
        REQUIRE(vertices.has_key(VertexKey(tet.vertices[2])));
        REQUIRE(vertices.has_key(VertexKey(tet.vertices[3])));

        auto validate_adjacency = [&](uint8_t local_index) {
            if (tet.mirrors[local_index] == mtet::invalid_key) return;

            auto key2 = TetKey(tet.mirrors[local_index]);
            REQUIRE(tets.has_key(key2));
            const auto& tet2 = *tets.get(key2);

            uint8_t sum = 0;
            for (uint8_t i = 0; i < 4; i++) {
                uint8_t j = mtet::get_tag(key2, i);
                sum += j;
                if (i != local_index) {
                    REQUIRE(tet.vertices[i] == tet2.vertices[j]);
                } else {
                    REQUIRE(tet.vertices[i] != tet2.vertices[j]);
                }
            }
            REQUIRE(sum == 6);
        };

        validate_adjacency(0);
        validate_adjacency(1);
        validate_adjacency(2);
        validate_adjacency(3);
    });
}

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

TEST_CASE("slot_map", "[slot_map]")
{
    using TetMap = mtet::MTetMeshImpl::TetMap;
    TetMap tet_map;

    SECTION("Empty")
    {
        REQUIRE(tet_map.empty());
        REQUIRE(tet_map.size() == 0);
        auto begin = tet_map.begin();
        auto end = tet_map.end();
        REQUIRE(begin == end);
    }

    SECTION("Empty2")
    {
        auto key = tet_map.emplace();
        tet_map.erase(key);

        REQUIRE(tet_map.empty());
        REQUIRE(tet_map.size() == 0);
        auto begin = tet_map.begin();
        auto end = tet_map.end();
        REQUIRE(begin == end);
    }

    SECTION("Repeated add/remove")
    {
        auto key = tet_map.emplace(mtet::MTet());
        REQUIRE(tet_map.has_key(key));

        TetMap::key key2 = key;
        for (size_t i = 0; i < 3000; i++) {
            tet_map.emplace(mtet::MTet());
            REQUIRE(tet_map.has_key(key2));
            tet_map.erase(key2);
            key2 = tet_map.emplace(mtet::MTet());
            tet_map.emplace(mtet::MTet());
        }
    }
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

    validate_mesh(mesh);
}

TEST_CASE("split_edge", "[mtet]")
{
    mtet::MTetMeshImpl mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    auto t0 = mesh.add_tet(v0, v1, v2, v3);
    mesh.initialize_connectivity();

    SECTION("case 1")
    {
        auto vid = mesh.split_edge(t0, 0);
        validate_mesh(mesh);
        REQUIRE(mesh.has_vertex(vid));
        REQUIRE(!mesh.has_tet(t0));
    }

    SECTION("case 2")
    {
        auto v4 = mesh.add_vertex(1, 1, 1);
        auto t1 = mesh.add_tet(v4, v3, v2, v1);
        mesh.initialize_connectivity();
        validate_mesh(mesh);
        auto vid = mesh.split_edge(t0, 1);
        REQUIRE(mesh.has_vertex(vid));
        REQUIRE(mesh.get_num_tets() == 4);
        validate_mesh(mesh);
    }

    SECTION("case 3")
    {
        // split in the loop.
        mesh.seq_foreach_tet([&](uint64_t tet_id, std::span<const uint64_t, 4>) {
            if (mesh.has_tet(tet_id)) {
                mesh.split_edge(mtet::MTetMeshImpl::TetKey(tet_id), 1);
            }
        });
        validate_mesh(mesh);
    }

    SECTION("case 4")
    {
        // Repeated split
        std::vector<uint64_t> tet_ids;
        for (size_t i = 0; i < 10; i++) {
            tet_ids.clear();
            size_t num_tets = mesh.get_num_tets();
            REQUIRE(num_tets > 0);
            tet_ids.reserve(num_tets);
            mesh.seq_foreach_tet(
                [&](uint64_t tet_id, std::span<const uint64_t, 4>) { tet_ids.push_back(tet_id); });
            for (auto tet_id : tet_ids) {
                if (mesh.has_tet(tet_id)) {
                    mesh.split_edge(mtet::MTetMeshImpl::TetKey(tet_id), i % 6);
                    validate_mesh(mesh);
                }
            }
        }
    }
}

