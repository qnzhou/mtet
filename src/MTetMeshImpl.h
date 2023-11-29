#pragma once

#include <array>
#include <cassert>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include <ankerl/unordered_dense.h>
#include <nanothread/nanothread.h>
#include <slot_map.h>

#include <mtet/mtet.h>

namespace dr = drjit; // For nanothread

namespace mtet {

const uint64_t invalid_key = dod::slot_map<int>::key::invalid();

using MVertex = std::array<Scalar, 3>;

struct MTet
{
    uint64_t vertices[4]{invalid_key, invalid_key, invalid_key, invalid_key};
    uint64_t mirrors[4]{invalid_key, invalid_key, invalid_key, invalid_key};
};


struct TriangleHash
{
    using is_avalanching = void;
    using is_transparent = void;

    [[nodiscard]] auto operator()(std::array<uint64_t, 3> const& x) const noexcept -> uint64_t
    {
        ankerl::unordered_dense::hash<uint64_t> hash_fn;
        return hash_fn(x[0]) + hash_fn(x[1]) + hash_fn(x[2]);
    }
};

struct TriangleEqual
{
    using is_transparent = void;

    bool operator()(std::array<uint64_t, 3> const& lhs, std::array<uint64_t, 3> const& rhs)
        const noexcept
    {
        return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]) ||
               (lhs[0] == rhs[1] && lhs[1] == rhs[2] && lhs[2] == rhs[0]) ||
               (lhs[0] == rhs[2] && lhs[1] == rhs[0] && lhs[2] == rhs[1]) ||
               (lhs[0] == rhs[0] && lhs[1] == rhs[2] && lhs[2] == rhs[1]) ||
               (lhs[0] == rhs[1] && lhs[1] == rhs[0] && lhs[2] == rhs[2]) ||
               (lhs[0] == rhs[2] && lhs[1] == rhs[1] && lhs[2] == rhs[0]);
    }
};

inline void set_tag(dod::slot_map<MTet>::key& key, uint8_t local_index, uint8_t value)
{
    auto tag = key.get_tag();
    uint16_t mask = 3 << (10 - local_index * 2);
    tag &= ~mask;
    tag |= (value << (10 - local_index * 2));
    key.set_tag(tag);
}

inline uint8_t get_tag(dod::slot_map<MTet>::key key, uint8_t local_index)
{
    auto tag = key.get_tag();
    uint16_t mask = 3 << (10 - local_index * 2);
    return static_cast<uint8_t>((tag & mask) >> (10 - local_index * 2));
}

inline uint32_t get_index(dod::slot_map<MTet>::key key)
{
    return dod::slot_map<MTet>::key::toIndex(key);
}

inline bool is_same_key(dod::slot_map<MTet>::key key1, dod::slot_map<MTet>::key key2)
{
    // Zero out the tag bits before comparison.
    key1.set_tag(0);
    key2.set_tag(0);
    return key1 == key2;
}

class MTetMeshImpl
{
public:
    using VertexMap = dod::slot_map<MVertex>;
    using TetMap = dod::slot_map<MTet>;
    using VertexKey = VertexMap::key;

    /// TetKey allocates 12 bits to store user-defined tags. We use the following convention:
    /// * bits [10, 12): local index of the adjacent tet that mirrors to corner 0.
    /// * bits [8, 10) : local index of the adjacent tet that mirrors to corner 1.
    /// * bits [6, 8)  : local index of the adjacent tet that mirrors to corner 2.
    /// * bits [4, 6)  : local index of the adjacent tet that mirrors to corner 3.
    /// * bits [0, 4)  : local index of an edge.
    using TetKey = TetMap::key;

public:
    uint64_t add_vertex(Scalar x, Scalar y, Scalar z)
    {
        return m_vertices.emplace(MVertex({{x, y, z}}));
    }

    uint64_t add_tet(uint64_t v0, uint64_t v1, uint64_t v2, uint64_t v3)
    {
        assert(v0 != v1);
        assert(v0 != v2);
        assert(v0 != v3);
        assert(v1 != v2);
        assert(v1 != v3);
        assert(v2 != v3);
        assert(m_vertices.has_key(VertexKey{v0}));
        assert(m_vertices.has_key(VertexKey{v1}));
        assert(m_vertices.has_key(VertexKey{v2}));
        assert(m_vertices.has_key(VertexKey{v3}));
        MTet tet;
        tet.vertices[0] = v0;
        tet.vertices[1] = v1;
        tet.vertices[2] = v2;
        tet.vertices[3] = v3;
        return m_tets.emplace(std::move(tet));
    }

    void initialize_connectivity()
    {
        const size_t n = m_tets.size();

        // Step 1: build a triangle->tet_id map.
        using HashMap = ankerl::unordered_dense::
            map<std::array<uint64_t, 3>, std::array<uint64_t, 2>, TriangleHash, TriangleEqual>;
        HashMap triangle_map;
        triangle_map.reserve(3 * n);
        std::vector<TetKey> tet_ids;
        tet_ids.reserve(n);

        std::array<uint64_t, 3> f123, f032, f021, f013;
        for (const auto& [tet_id, tet_ref] : m_tets.items()) {
            tet_ids.push_back(tet_id);
            const MTet& tet = tet_ref.get();
            f123 = {tet.vertices[1], tet.vertices[2], tet.vertices[3]};
            f032 = {tet.vertices[0], tet.vertices[3], tet.vertices[2]};
            f021 = {tet.vertices[0], tet.vertices[2], tet.vertices[1]};
            f013 = {tet.vertices[0], tet.vertices[1], tet.vertices[3]};

            auto itr_123 = triangle_map.find(f123);
            if (itr_123 == triangle_map.end()) {
                triangle_map.insert({f123, {tet_id, invalid_key}});
            } else {
                assert(itr_123->second[1] == invalid_key);
                itr_123->second[1] = tet_id;
            }

            auto itr_032 = triangle_map.find(f032);
            if (itr_032 == triangle_map.end()) {
                triangle_map.insert({f032, {tet_id, invalid_key}});
            } else {
                assert(itr_032->second[1] == invalid_key);
                itr_032->second[1] = tet_id;
            }

            auto itr_021 = triangle_map.find(f021);
            if (itr_021 == triangle_map.end()) {
                triangle_map.insert({f021, {tet_id, invalid_key}});
            } else {
                assert(itr_021->second[1] == invalid_key);
                itr_021->second[1] = tet_id;
            }

            auto itr_013 = triangle_map.find(f013);
            if (itr_013 == triangle_map.end()) {
                triangle_map.insert({f013, {tet_id, invalid_key}});
            } else {
                assert(itr_013->second[1] == invalid_key);
                itr_013->second[1] = tet_id;
            }
        }

        // Step 2: build mirror mapping between adjacent tets.

        /**
         * Compute the mirror mapping of an adjacent tet with the current tet.
         *
         * @param tet_id     The id of the current tet.
         * @param local_fid  The local face id of the current tet, which is shared with the adjcent
         *                   tet.
         *
         * The current tet will be updated with the correct mirror mapping for this adjacency.
         */
        auto process_tet_triangle = [&](TetKey tet_id, uint8_t local_fid) {
            MTet& tet = *m_tets.get(tet_id);

            // Get the common triangle face.
            const uint8_t i0 = local_fid;
            const uint8_t i1 = (i0 + 1) % 4;
            const uint8_t i2 = (i0 + 2) % 4;
            const uint8_t i3 = (i0 + 3) % 4;
            std::array<uint64_t, 3> tri = {tet.vertices[i1], tet.vertices[i2], tet.vertices[i3]};
            assert(triangle_map.contains(tri));

            // Retrieve the adjacent tet id.
            const auto& adj_tets = triangle_map[tri];
            assert(adj_tets[0] == tet_id || adj_tets[1] == tet_id);
            TetKey adj_tet_id(adj_tets[0] == tet_id ? adj_tets[1] : adj_tets[0]);

            // Update adjacency.
            if (adj_tet_id != invalid_key) {
                const MTet& adj_tet = *m_tets.get(adj_tet_id);

                // Update tag to store mirror corner mapping.
                for (uint8_t i = 0; i < 4; i++) {
                    for (uint8_t j = 0; j < 4; j++) {
                        if (tet.vertices[i] == adj_tet.vertices[j]) {
                            set_tag(adj_tet_id, i, j);
                        }
                    }
                }
                tet.mirrors[i0] = adj_tet_id;
            } else {
                tet.mirrors[i0] = invalid_key;
            }
        };

        /**
         * Compute `tet.mirrors` for the current tet.
         *
         * @param tet_id  The id of the current tet.
         */
        auto process_tet = [&](TetKey tet_id) {
            assert(m_tets.has_key(tet_id));

            process_tet_triangle(tet_id, 0);
            process_tet_triangle(tet_id, 1);
            process_tet_triangle(tet_id, 2);
            process_tet_triangle(tet_id, 3);
        };

        dr::parallel_for(dr::blocked_range<size_t>(0, n, 1), [&](dr::blocked_range<size_t> range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                process_tet(tet_ids[i]);
            }
        });
    }

public:
    bool has_vertex(uint64_t vertex_id) const
    {
        VertexKey key{vertex_id};
        return m_vertices.has_key(key);
    }

    bool has_tet(uint64_t tet_id) const
    {
        TetKey key{tet_id};
        return m_tets.has_key(key);
    }

    std::span<Scalar, 3> get_vertex(uint64_t vertex_id)
    {
        VertexKey key{vertex_id};
        auto ptr = m_vertices.get(key);
        if (ptr != nullptr) {
            return *ptr;
        } else {
            throw std::runtime_error("Vertex not found");
        }
    }

    std::span<const Scalar, 3> get_vertex(uint64_t vertex_id) const
    {
        VertexKey key{vertex_id};
        auto ptr = m_vertices.get(key);
        if (ptr != nullptr) {
            return *ptr;
        } else {
            throw std::runtime_error("Vertex not found");
        }
    }

    std::span<uint64_t, 4> get_tet(uint64_t tet_id)
    {
        TetKey key{tet_id};
        auto ptr = m_tets.get(key);
        if (ptr != nullptr) {
            return std::span<uint64_t, 4>(ptr->vertices, 4);
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

    std::span<const uint64_t, 4> get_tet(uint64_t tet_id) const
    {
        TetKey key{tet_id};
        auto ptr = m_tets.get(key);
        if (ptr != nullptr) {
            return std::span<const uint64_t, 4>(ptr->vertices, 4);
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

    size_t get_num_vertices() const { return m_vertices.size(); }
    size_t get_num_tets() const { return m_tets.size(); }

    uint64_t split_edge(uint64_t tet_id, uint8_t local_index)
    {
        constexpr std::array<std::array<uint8_t, 2>, 6> edge_map = {
            {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}}};

        auto get_next_tet_id =
            [&](TetKey curr_tet_id, uint8_t llv0, uint8_t llv1, uint8_t llv2, uint8_t llv3) {
                assert(m_tets.has_key(curr_tet_id));
                assert(llv0 < 4);
                assert(llv1 < 4);
                assert(llv2 < 4);
                assert(llv3 < 4);

                const auto& curr_tet = *m_tets.get(curr_tet_id);
                TetKey next_tet_id = TetKey(curr_tet.mirrors[llv3]);
                llv0 = get_tag(next_tet_id, llv0);
                llv1 = get_tag(next_tet_id, llv1);
                llv2 = get_tag(next_tet_id, llv2);
                llv3 = get_tag(next_tet_id, llv3);

                // Note that we swap llv2 and llv3 here intentionally.
                return std::make_tuple(next_tet_id, llv0, llv1, llv3, llv2);
            };


        /**
         * Split the current tet along the edge [llv0, llv1] at vm.
         */
        auto split_tet = [&](TetKey curr_tet_id,
                             uint64_t vm,
                             uint8_t llv0,
                             uint8_t llv1,
                             uint8_t llv2,
                             uint8_t llv3) -> std::array<TetKey, 2> {
            assert(m_tets.has_key(curr_tet_id));
            const auto& curr_tet = *m_tets.get(curr_tet_id);

            const auto o0 = TetKey(curr_tet.mirrors[llv0]);
            const auto o1 = TetKey(curr_tet.mirrors[llv1]);

            // Add two new tets.
            std::array<uint64_t, 4> tet_0_vertices{
                curr_tet.vertices[0],
                curr_tet.vertices[1],
                curr_tet.vertices[2],
                curr_tet.vertices[3],
            };
            tet_0_vertices[llv1] = vm;
            std::array<uint64_t, 4> tet_1_vertices{
                curr_tet.vertices[0],
                curr_tet.vertices[1],
                curr_tet.vertices[2],
                curr_tet.vertices[3],
            };
            tet_1_vertices[llv0] = vm;

            auto t0 = TetKey(add_tet(
                tet_0_vertices[0],
                tet_0_vertices[1],
                tet_0_vertices[2],
                tet_0_vertices[3]));
            auto t1 = TetKey(add_tet(
                tet_1_vertices[0],
                tet_1_vertices[1],
                tet_1_vertices[2],
                tet_1_vertices[3]));

            auto& tet_0 = *m_tets.get(t0);
            auto& tet_1 = *m_tets.get(t1);

            // Update mirror information for new tets.
            set_tag(t1, llv0, llv1);
            set_tag(t1, llv1, llv0);
            set_tag(t1, llv2, llv2);
            set_tag(t1, llv3, llv3);

            set_tag(t0, llv0, llv1);
            set_tag(t0, llv1, llv0);
            set_tag(t0, llv2, llv2);
            set_tag(t0, llv3, llv3);

            tet_0.mirrors[llv0] = t1;
            tet_0.mirrors[llv1] = o1;
            tet_0.mirrors[llv2] = invalid_key;
            tet_0.mirrors[llv3] = invalid_key;

            tet_1.mirrors[llv0] = t0;
            tet_1.mirrors[llv1] = o0;
            tet_1.mirrors[llv2] = invalid_key;
            tet_1.mirrors[llv3] = invalid_key;

            // Update mirror information for tet adjacent to the new tets.
            if (o0 != invalid_key) {
                auto& tet_o0 = *m_tets.get(o0);
                auto o0_llv0 = get_tag(o0, llv0);
                auto old_t0_mirror = TetKey(tet_o0.mirrors[o0_llv0]);
                t1.set_tag(old_t0_mirror.get_tag());
                tet_o0.mirrors[o0_llv0] = t1;
            }
            if (o1 != invalid_key) {
                auto& tet_o1 = *m_tets.get(o1);
                auto o1_llv1 = get_tag(o1, llv1);
                auto old_t1_mirror = TetKey(tet_o1.mirrors[o1_llv1]);
                t0.set_tag(old_t1_mirror.get_tag());
                tet_o1.mirrors[o1_llv1] = t0;
            }

            return {t0, t1};
        };

        // Compute local indices.
        // lv0 and lv1 are vertices of the edge with local_index.
        // lv2 and lv3 are the indices of vertices in the opposite edge.
        uint8_t lv0 = edge_map[local_index][0];
        uint8_t lv1 = edge_map[local_index][1];
        uint8_t lv2 = ~(lv0 & lv1) & 3;
        uint8_t lv3 = (lv0 ^ lv1 ^ lv2) & 3;

        TetKey key{tet_id};
        auto ptr = m_tets.get(key);
        if (ptr != nullptr) {
            // Create new vertex
            const auto& curr_tet = *ptr;
            const auto v0 = curr_tet.vertices[lv0];
            const auto v1 = curr_tet.vertices[lv1];
            const auto p0 = get_vertex(v0);
            const auto p1 = get_vertex(v1);
            const auto vm =
                add_vertex((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2, (p0[2] + p1[2]) / 2);

            // Compute the staring tet if the edge is on the boundary.
            TetKey curr_key = key;
            bool on_boundary = false;
            do {
                assert(curr_key != invalid_key);

                TetKey next_key;
                uint8_t llv0, llv1, llv2, llv3;
                // Note that we are traversing in the lv2 direction. I.e. the next tet is the tet
                // opposite to the lv2 vertex.
                std::tie(next_key, llv0, llv1, llv3, llv2) =
                    get_next_tet_id(curr_key, lv0, lv1, lv3, lv2);
                if (next_key == invalid_key) {
                    on_boundary = true;
                    break;
                }
                curr_key = next_key;
                lv0 = llv0;
                lv1 = llv1;
                lv2 = llv2;
                lv3 = llv3;
            } while (is_same_key(curr_key, key));
            // Current key should hold a valid key that serve as the first tet to traverse in the
            // lv3 direction. lv0 to lv3 should be valid local indices for curr_key.
            assert(m_tets.has_key(curr_key));

            // Gather 1-ring tets around the edge
            std::vector<TetKey> old_one_ring, new_one_ring_0, new_one_ring_1;
            std::vector<uint8_t> local_indices;
            old_one_ring.reserve(16);
            new_one_ring_0.reserve(16);
            new_one_ring_1.reserve(16);
            local_indices.reserve(16 * 4);

            do {
                old_one_ring.push_back(curr_key);
                auto [t0, t1] = split_tet(curr_key, vm, lv0, lv1, lv2, lv3);
                new_one_ring_0.push_back(t0);
                new_one_ring_1.push_back(t1);
                local_indices.push_back(lv0);
                local_indices.push_back(lv1);
                local_indices.push_back(lv2);
                local_indices.push_back(lv3);
                std::tie(curr_key, lv0, lv1, lv2, lv3) =
                    get_next_tet_id(curr_key, lv0, lv1, lv2, lv3);
            } while (curr_key != invalid_key && is_same_key(curr_key, key));
            size_t one_ring_size = old_one_ring.size();

            for (size_t i = 0; i < one_ring_size; i++) {
                if (on_boundary && i == one_ring_size - 1) break;
                const size_t j = (i + 1) % one_ring_size;

                const auto curr_lv0 = local_indices[i * 4 + 0];
                const auto curr_lv1 = local_indices[i * 4 + 1];
                const auto curr_lv2 = local_indices[i * 4 + 2];
                const auto curr_lv3 = local_indices[i * 4 + 3];

                const auto next_lv0 = local_indices[j * 4 + 0];
                const auto next_lv1 = local_indices[j * 4 + 1];
                const auto next_lv2 = local_indices[j * 4 + 2];
                const auto next_lv3 = local_indices[j * 4 + 3];

                auto& curr_tet_id_0 = new_one_ring_0[i];
                auto& next_tet_id_0 = new_one_ring_0[j];
                auto& curr_tet_0 = *m_tets.get(curr_tet_id_0);
                auto& next_tet_0 = *m_tets.get(next_tet_id_0);

                assert(curr_tet_0.vertices[curr_lv0] == next_tet_0.vertices[next_lv0]);
                assert(curr_tet_0.vertices[curr_lv1] == next_tet_0.vertices[next_lv1]);
                assert(curr_tet_0.vertices[curr_lv2] == next_tet_0.vertices[next_lv3]);

                set_tag(next_tet_id_0, curr_lv0, next_lv0);
                set_tag(next_tet_id_0, curr_lv1, next_lv1);
                set_tag(next_tet_id_0, curr_lv2, next_lv3);
                set_tag(next_tet_id_0, curr_lv3, next_lv2);
                curr_tet_0.mirrors[curr_lv3] = next_tet_id_0;

                set_tag(curr_tet_id_0, next_lv0, curr_lv0);
                set_tag(curr_tet_id_0, next_lv1, curr_lv1);
                set_tag(curr_tet_id_0, next_lv2, curr_lv3);
                set_tag(curr_tet_id_0, next_lv3, curr_lv2);
                next_tet_0.mirrors[next_lv2] = curr_tet_id_0;

                auto& curr_tet_id_1 = new_one_ring_1[i];
                auto& next_tet_id_1 = new_one_ring_1[j];
                auto& curr_tet_1 = *m_tets.get(curr_tet_id_1);
                auto& next_tet_1 = *m_tets.get(next_tet_id_1);

                assert(curr_tet_1.vertices[curr_lv0] == next_tet_1.vertices[next_lv0]);
                assert(curr_tet_1.vertices[curr_lv1] == next_tet_1.vertices[next_lv1]);
                assert(curr_tet_1.vertices[curr_lv2] == next_tet_1.vertices[next_lv3]);

                set_tag(next_tet_id_1, curr_lv0, next_lv0);
                set_tag(next_tet_id_1, curr_lv1, next_lv1);
                set_tag(next_tet_id_1, curr_lv2, next_lv3);
                set_tag(next_tet_id_1, curr_lv3, next_lv2);
                curr_tet_1.mirrors[curr_lv3] = next_tet_id_1;

                set_tag(curr_tet_id_1, next_lv0, curr_lv0);
                set_tag(curr_tet_id_1, next_lv1, curr_lv1);
                set_tag(curr_tet_id_1, next_lv2, curr_lv3);
                set_tag(curr_tet_id_1, next_lv3, curr_lv2);
                next_tet_1.mirrors[next_lv2] = curr_tet_id_1;
            }

            for (auto tid : old_one_ring) {
                m_tets.erase(tid);
            }

            return vm;
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

public:
    /**
     * Internal: get the adjacent tet of a tet across one of its triangles.
     *
     * @param key          The key of the current tet.
     * @param local_index  The local index of the triangle face of the current tet.
     *
     * @return The key of the adjacent tet.
     */
    uint64_t get_adjacent_tet(uint64_t key, uint8_t local_index) const
    {
        auto ptr = m_tets.get(static_cast<TetKey>(key));
        if (ptr != nullptr) {
            return ptr->mirrors[local_index];
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

private:
    VertexMap m_vertices;
    TetMap m_tets;
};

} // namespace mtet
