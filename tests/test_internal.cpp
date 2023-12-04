#include <catch2/catch_test_macros.hpp>

#include <MTetMeshImpl.h>

#include <array>

void validate_mesh(const mtet::MTetMeshImpl& mesh)
{
    using VertexId = mtet::VertexId;
    using TetId = mtet::TetId;

    mesh.par_foreach_tet([&](TetId tet_id, std::span<const VertexId, 4> tet_vertices) {
        using VertexKey = mtet::MTetMeshImpl::VertexKey;
        using TetKey = mtet::MTetMeshImpl::TetKey;
        TetKey key(value_of(tet_id));

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

        REQUIRE(vertices.has_key(VertexKey(value_of(tet.vertices[0]))));
        REQUIRE(vertices.has_key(VertexKey(value_of(tet.vertices[1]))));
        REQUIRE(vertices.has_key(VertexKey(value_of(tet.vertices[2]))));
        REQUIRE(vertices.has_key(VertexKey(value_of(tet.vertices[3]))));

        auto validate_adjacency = [&](uint8_t local_index) {
            if (value_of(tet.mirrors[local_index]) == mtet::invalid_key) return;

            auto key2 = TetKey(value_of(tet.mirrors[local_index]));
            REQUIRE(tets.has_key(key2));
            const auto& tet2 = *tets.get(key2);

            uint8_t sum = 0;
            for (uint8_t i = 0; i < 4; i++) {
                uint8_t j = mtet::get_mirror_index(key2, i);
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
    mtet::VertexId v0(0);
    mtet::VertexId v1(1);
    mtet::VertexId v2(2);
    mtet::VertexId v3(3);
    std::array<mtet::VertexId, 3> t0 = {{v0, v1, v2}};
    std::array<mtet::VertexId, 3> t1 = {{v2, v0, v1}};
    std::array<mtet::VertexId, 3> t2 = {{v2, v1, v0}};
    std::array<mtet::VertexId, 3> t3 = {{v1, v1, v1}};

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
    mtet::VertexId v0(0);
    mtet::VertexId v1(1);
    mtet::VertexId v2(2);
    mtet::VertexId v3(3);
    std::array<mtet::VertexId, 3> t0 = {{v0, v1, v2}};
    std::array<mtet::VertexId, 3> t1 = {{v2, v0, v1}};
    std::array<mtet::VertexId, 3> t2 = {{v2, v1, v0}};
    std::array<mtet::VertexId, 3> t3 = {{v1, v1, v1}};

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
        mtet::set_mirror_index(key, i, i);
    }
    mtet::set_edge_index(key, 3);
    REQUIRE(mtet::get_edge_index(key) == 3);
    mtet::set_edge_index(key, 3, 2);

    for (uint8_t i = 0; i < 4; i++) {
        REQUIRE(mtet::get_mirror_index(key, i) == i);
    }
    REQUIRE(mtet::get_edge_index(key) == 5);

    REQUIRE(key != key_copy);
}

TEST_CASE("unordered_dense", "[unordered_dense]")
{
    using HashMap = ankerl::unordered_dense::map<
        std::array<mtet::VertexId, 3>,
        std::array<mtet::TetId, 2>,
        mtet::TriangleHash,
        mtet::TriangleEqual>;

    mtet::VertexId v0(0);
    mtet::VertexId v1(1);
    mtet::VertexId v2(2);
    mtet::VertexId v3(3);
    std::array<mtet::VertexId, 3> t0 = {{v0, v1, v2}};
    std::array<mtet::VertexId, 3> t1 = {{v2, v0, v1}};
    std::array<mtet::VertexId, 3> t2 = {{v2, v1, v0}};
    std::array<mtet::VertexId, 3> t3 = {{v1, v1, v1}};

    mtet::TetId tet_0(0);
    mtet::TetId tet_1(1);

    HashMap triangle_map;
    triangle_map[t0] = {{tet_0, tet_1}};
    REQUIRE(triangle_map.contains(t0));
    REQUIRE(triangle_map.contains(t1));
    REQUIRE(triangle_map.contains(t2));
    REQUIRE(!triangle_map.contains(t3));
    REQUIRE(is_same_tet(triangle_map[t0][0], tet_0));
    REQUIRE(is_same_tet(triangle_map[t0][1], tet_1));
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
    using TetKey = mtet::MTetMeshImpl::TetKey;

    mtet::MTetMeshImpl mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    auto t0 = mesh.add_tet(v0, v1, v2, v3);
    TetKey t0_key(value_of(t0));
    REQUIRE(mesh.has_vertex(v0));
    REQUIRE(mesh.has_vertex(v1));
    REQUIRE(mesh.has_vertex(v2));
    REQUIRE(mesh.has_vertex(v3));
    REQUIRE(mesh.has_tet(t0));
    REQUIRE(mesh.get_num_vertices() == 4);
    REQUIRE(mesh.get_num_tets() == 1);

    mesh.initialize_connectivity();

    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 0)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 1)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 2)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 3)));

    // Add second tet.
    auto v4 = mesh.add_vertex(1, 1, 1);
    auto t1 = mesh.add_tet(v4, v3, v2, v1);
    TetKey t1_key(value_of(t1));
    REQUIRE(mesh.has_vertex(v4));
    REQUIRE(mesh.has_tet(t1));
    mesh.initialize_connectivity();

    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 1)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 2)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t0, 3)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t1, 1)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t1, 2)));
    REQUIRE(mtet::is_invalid_tet(mesh.get_adjacent_tet(t1, 3)));

    // The opposite tet of t0 across face 0 is t1.
    TetKey t0_0_key(value_of(mesh.get_adjacent_tet(t0, 0)));
    TetKey t1_0_key(value_of(mesh.get_adjacent_tet(t1, 0)));

    REQUIRE(TetKey::toIndex(t0_0_key) == TetKey::toIndex(t1_key));
    REQUIRE(TetKey::toIndex(t1_0_key) == TetKey::toIndex(t0_key));

    REQUIRE(mtet::get_mirror_index(t0_0_key, 0) == 0);
    REQUIRE(mtet::get_mirror_index(t0_0_key, 1) == 3);
    REQUIRE(mtet::get_mirror_index(t0_0_key, 2) == 2);
    REQUIRE(mtet::get_mirror_index(t0_0_key, 3) == 1);

    REQUIRE(mtet::get_mirror_index(t1_0_key, 0) == 0);
    REQUIRE(mtet::get_mirror_index(t1_0_key, 1) == 3);
    REQUIRE(mtet::get_mirror_index(t1_0_key, 2) == 2);
    REQUIRE(mtet::get_mirror_index(t1_0_key, 3) == 1);

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
        auto [vid, e0, e1] = mesh.split_edge(t0, 0);
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
        auto r = mesh.split_edge(t0, 1);
        auto vid = std::get<0>(r);
        REQUIRE(mesh.has_vertex(vid));
        REQUIRE(mesh.get_num_tets() == 4);
        validate_mesh(mesh);
    }

    SECTION("case 3")
    {
        // split in the loop.
        mesh.seq_foreach_tet([&](mtet::TetId tet_id, std::span<const mtet::VertexId, 4>) {
            if (mesh.has_tet(tet_id)) {
                mesh.split_edge(tet_id, 1);
            }
        });
        validate_mesh(mesh);
    }

    SECTION("case 4")
    {
        // Repeated split
        std::vector<mtet::TetId> tet_ids;
        for (size_t i = 0; i < 10; i++) {
            tet_ids.clear();
            size_t num_tets = mesh.get_num_tets();
            REQUIRE(num_tets > 0);
            tet_ids.reserve(num_tets);
            mesh.seq_foreach_tet([&](mtet::TetId tet_id, std::span<const mtet::VertexId, 4>) {
                tet_ids.push_back(tet_id);
            });
            for (auto tet_id : tet_ids) {
                if (mesh.has_tet(tet_id)) {
                    mesh.split_edge(tet_id, i % 6);
                    validate_mesh(mesh);
                }
            }
        }
    }
}

TEST_CASE("foreach_edge_in_tet", "[mtet]")
{
    mtet::MTetMeshImpl mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    auto t0 = mesh.add_tet(v0, v1, v2, v3);
    mesh.initialize_connectivity();

    mesh.foreach_edge_in_tet(t0, [&](mtet::EdgeId edge_id, mtet::VertexId v0, mtet::VertexId v1) {
        REQUIRE(mesh.has_edge(edge_id));
        REQUIRE(mesh.has_vertex(v0));
        REQUIRE(mesh.has_vertex(v1));
        auto ev = mesh.get_edge_vertices(edge_id);
        REQUIRE(ev[0] == v0);
        REQUIRE(ev[1] == v1);
        auto tet_id = mesh.get_edge_tet(edge_id);
        REQUIRE(mesh.has_tet(tet_id));
    });
}

