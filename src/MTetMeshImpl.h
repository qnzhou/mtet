#pragma once

#include <array>
#include <cassert>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include <SmallVector.h>
#include <ankerl/unordered_dense.h>
#include <nanothread/nanothread.h>
#include <slot_map.h>

#include <mtet/mtet.h>

namespace dr = drjit; // For nanothread

namespace mtet {

/// Invalid key for slot_map.
/// This is equivalent to dod::slot_map_key64<T>::invalid() for any `T`.
constexpr uint64_t invalid_key = 0;

using MVertex = std::array<Scalar, 3>;

constexpr VertexId invalid_vertex_id = VertexId(invalid_key);
constexpr TetId invalid_tet_id = TetId(invalid_key);

constexpr std::array<std::array<uint8_t, 4>, 6> edge_map = {
    {{0, 1, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {0, 3, 1, 2}, {1, 3, 2, 0}, {2, 3, 0, 1}}};

/**
 * Tet data structrue.
 */
struct MTet
{
    /// Vertex keys of a tet. The vertex order determines the tet orientation.
    /// vertex[0] should be on the negative side of the face (vertex[1], vertex[2], vertex[3]).
    VertexId vertices[4]{
        invalid_vertex_id,
        invalid_vertex_id,
        invalid_vertex_id,
        invalid_vertex_id};

    /// Tet keys of the mirror tets.
    ///
    /// The tag of each tet key is used to store mirror indices of the local indices of
    /// the current tet. Use `set_mirror_index` and `get_mirror_index` to access the
    /// mirror indices.
    TetId mirrors[4]{invalid_tet_id, invalid_tet_id, invalid_tet_id, invalid_tet_id};
};


/**
 * Triangle hash function.
 *
 * This hash function is agnostic to the vertex order and triangle orientation. I.e. all the of the
 * following triangles will have the same hash: [0, 1, 2], [1, 2, 0], [2, 1, 0].
 */
struct TriangleHash
{
    using is_avalanching = void;
    using is_transparent = void;

    [[nodiscard]] auto operator()(std::array<VertexId, 3> const& x) const noexcept -> uint64_t
    {
        ankerl::unordered_dense::hash<uint64_t> hash_fn;
        return hash_fn(value_of(x[0])) + hash_fn(value_of(x[1])) + hash_fn(value_of(x[2]));
    }
};

/**
 * Triangle equality function.
 *
 * This hash function is agnostic to the vertex order and triangle orientation. I.e. all the of the
 * following triangles are considered as equal: [0, 1, 2], [1, 2, 0], [2, 1, 0].
 */
struct TriangleEqual
{
    using is_transparent = void;

    bool operator()(std::array<VertexId, 3> const& lhs, std::array<VertexId, 3> const& rhs)
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


/**
 * Set the mirror local index in a key.
 *
 * @param key[in/out] The target tet id.
 * @param local_index The local index of the current tet, which is adjacent to the target tet.
 * @param value       The local index of the target tet that mirrors to the `local_index` in the
 *                    current tet.
 *
 * The target `key` will be updated with the new mirror index.
 */
inline void set_mirror_index(dod::slot_map<MTet>::key& key, uint8_t local_index, uint8_t value)
{
    assert(local_index < 4);
    assert(value < 4);
    auto tag = key.get_tag();
    uint16_t mask = 3 << (10 - local_index * 2);
    tag &= ~mask;
    tag |= (value << (10 - local_index * 2));
    key.set_tag(tag);
}

/**
 * Get the mirror local index stored in a key.
 *
 * @param key[in]       The target tet id.
 * @param local_index   The local index of the current tet, which is adjacent to the target tet.
 *
 * @return The local index of the target tet that mirrors to the `local_index` in the current tet.
 */
inline uint8_t get_mirror_index(dod::slot_map<MTet>::key key, uint8_t local_index)
{
    assert(local_index < 4);
    auto tag = key.get_tag();
    uint16_t mask = 3 << (10 - local_index * 2);
    return static_cast<uint8_t>((tag & mask) >> (10 - local_index * 2));
}

/**
 * Set the local edge index in a key.
 *
 * @param key[in/out] The target tet id.
 * @param edge_index   The local edge index of the target tet.
 */
inline void set_edge_index(dod::slot_map<MTet>::key& key, uint8_t edge_index)
{
    assert(edge_index < 6);
    auto tag = key.get_tag();
    constexpr uint16_t mask = 15;
    tag &= ~mask;
    tag |= edge_index;
    key.set_tag(tag);
}

/**
 * Se the local edge index in a key.
 *
 * @param key[in/out] The target tet id.
 * @param lv0          The first local vertex index of the edge.
 * @param lv1          The second local vertex index of the edge.
 *
 * Here is the edge index to local vertices index mapping:
 *   0: [v0, v1],
 *   1: [v1, v2],
 *   2: [v2, v0],
 *   3: [v0, v3],
 *   4: [v1, v3],
 *   5: [v2, v3]
 */
inline void set_edge_index(dod::slot_map<MTet>::key& key, uint8_t lv0, uint8_t lv1)
{
    constexpr uint8_t invalid_edge_index = 0xff;
    constexpr uint8_t edge_index_table[13]{
        invalid_edge_index, // 0: 0b0000
        invalid_edge_index, // 1: 0b0001
        invalid_edge_index, // 2: 0b0010
        0, // 3: 0b0011
        invalid_edge_index, // 4: 0b0100
        2, // 5: 0b0101
        1, // 6: 0b0110
        invalid_edge_index, // 7: 0b0111
        invalid_edge_index, // 8: 0b1000
        3, // 9: 0b1001
        4, // 10: 0b1010
        invalid_edge_index, // 11: 0b1011
        5 // 12: 0b1100
    };

    assert(lv0 < 4);
    assert(lv1 < 4);
    uint8_t bitmap = 0 | (1 << lv0) | (1 << lv1);
    assert(bitmap < 13);
    assert(edge_index_table[bitmap] != invalid_edge_index);
    set_edge_index(key, edge_index_table[bitmap]);
}

/**
 * Get the local edge index stored in a key.
 */
inline uint8_t get_edge_index(dod::slot_map<MTet>::key key)
{
    auto tag = key.get_tag();
    constexpr uint16_t mask = 15;
    return static_cast<uint8_t>(tag & mask);
}

/**
 * Check if two keys represent the same tet.
 *
 * @param key1 The first key.
 * @param key2 The second key.
 *
 * @return True if the two keys are the same.
 *
 * @note This is different from `key1 == key2`. This function ignores the tag bits.
 */
inline bool is_same_key(dod::slot_map<MTet>::key key1, dod::slot_map<MTet>::key key2)
{
    // Zero out the tag bits before comparison.
    key1.set_tag(0);
    key2.set_tag(0);
    return key1 == key2;
}

inline bool is_same_tet(TetId id1, TetId id2)
{
    using TetKey = dod::slot_map<MTet>::key;
    return is_same_key(TetKey(value_of(id1)), TetKey(value_of(id2)));
}

inline bool is_invalid_tet(TetId id)
{
    using TetKey = dod::slot_map<MTet>::key;
    return value_of(id) == invalid_key;
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
    VertexId add_vertex(Scalar x, Scalar y, Scalar z)
    {
        return VertexId(m_vertices.emplace(MVertex({{x, y, z}})));
    }

    TetId add_tet(VertexId v0, VertexId v1, VertexId v2, VertexId v3)
    {
        assert(v0 != v1);
        assert(v0 != v2);
        assert(v0 != v3);
        assert(v1 != v2);
        assert(v1 != v3);
        assert(v2 != v3);
        assert(m_vertices.has_key(VertexKey{value_of(v0)}));
        assert(m_vertices.has_key(VertexKey{value_of(v1)}));
        assert(m_vertices.has_key(VertexKey{value_of(v2)}));
        assert(m_vertices.has_key(VertexKey{value_of(v3)}));
        MTet tet;
        tet.vertices[0] = v0;
        tet.vertices[1] = v1;
        tet.vertices[2] = v2;
        tet.vertices[3] = v3;
        return TetId(m_tets.emplace(std::move(tet)));
    }

    void initialize_connectivity()
    {
        const size_t n = m_tets.size();

        // Step 1: build a triangle->tet_id map.
        using HashMap = ankerl::unordered_dense::
            map<std::array<VertexId, 3>, std::array<TetId, 2>, TriangleHash, TriangleEqual>;
        HashMap triangle_map;
        triangle_map.reserve(3 * n);
        std::vector<TetId> tet_ids;
        tet_ids.reserve(n);

        std::array<VertexId, 3> f123, f032, f021, f013;
        for (const auto& [tet_key, tet_ref] : m_tets.items()) {
            TetId tet_id(tet_key);
            tet_ids.push_back(tet_id);
            const MTet& tet = tet_ref.get();
            f123 = {tet.vertices[1], tet.vertices[2], tet.vertices[3]};
            f032 = {tet.vertices[0], tet.vertices[3], tet.vertices[2]};
            f021 = {tet.vertices[0], tet.vertices[2], tet.vertices[1]};
            f013 = {tet.vertices[0], tet.vertices[1], tet.vertices[3]};

            auto itr_123 = triangle_map.find(f123);
            if (itr_123 == triangle_map.end()) {
                triangle_map.insert({f123, {tet_id, invalid_tet_id}});
            } else {
                assert(itr_123->second[1] == invalid_tet_id);
                itr_123->second[1] = tet_id;
            }

            auto itr_032 = triangle_map.find(f032);
            if (itr_032 == triangle_map.end()) {
                triangle_map.insert({f032, {tet_id, invalid_tet_id}});
            } else {
                assert(itr_032->second[1] == invalid_tet_id);
                itr_032->second[1] = tet_id;
            }

            auto itr_021 = triangle_map.find(f021);
            if (itr_021 == triangle_map.end()) {
                triangle_map.insert({f021, {tet_id, invalid_tet_id}});
            } else {
                assert(itr_021->second[1] == invalid_tet_id);
                itr_021->second[1] = tet_id;
            }

            auto itr_013 = triangle_map.find(f013);
            if (itr_013 == triangle_map.end()) {
                triangle_map.insert({f013, {tet_id, invalid_tet_id}});
            } else {
                assert(itr_013->second[1] == invalid_tet_id);
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
        auto process_tet_triangle = [&](TetId tet_id, uint8_t local_fid) {
            TetKey tet_key(value_of(tet_id));
            MTet& tet = *m_tets.get(tet_key);

            // Get the common triangle face.
            const uint8_t i0 = local_fid;
            const uint8_t i1 = (i0 + 1) % 4;
            const uint8_t i2 = (i0 + 2) % 4;
            const uint8_t i3 = (i0 + 3) % 4;
            std::array<VertexId, 3> tri = {tet.vertices[i1], tet.vertices[i2], tet.vertices[i3]};
            assert(triangle_map.contains(tri));

            // Retrieve the adjacent tet id.
            const auto& adj_tets = triangle_map[tri];
            assert(adj_tets[0] == tet_id || adj_tets[1] == tet_id);
            TetId adj_tet_id(is_same_tet(adj_tets[0], tet_id) ? adj_tets[1] : adj_tets[0]);

            // Update adjacency.
            if (!is_same_tet(adj_tet_id, invalid_tet_id)) {
                TetKey adj_tet_key(value_of(adj_tet_id));
                const MTet& adj_tet = *m_tets.get(adj_tet_key);
                uint8_t sum = 0; // Local indices in mirror should sum to 6.

                // Update tag to store mirror corner mapping.
                for (uint8_t i = 0; i < 4; i++) {
                    for (uint8_t j = 0; j < 4; j++) {
                        if (tet.vertices[i] == adj_tet.vertices[j]) {
                            set_mirror_index(adj_tet_key, i, j);
                            sum += j;
                        }
                    }
                }
                set_mirror_index(adj_tet_key, local_fid, 6 - sum);
                assert(
                    get_mirror_index(adj_tet_key, 0) + get_mirror_index(adj_tet_key, 1) +
                        get_mirror_index(adj_tet_key, 2) + get_mirror_index(adj_tet_key, 3) ==
                    6);
                tet.mirrors[i0] = TetId(adj_tet_key);
            } else {
                tet.mirrors[i0] = invalid_tet_id;
            }
        };

        /**
         * Compute `tet.mirrors` for the current tet.
         *
         * @param tet_id  The id of the current tet.
         */
        auto process_tet = [&](TetId tet_id) {
            assert(has_tet(tet_id));

            process_tet_triangle(tet_id, 0);
            process_tet_triangle(tet_id, 1);
            process_tet_triangle(tet_id, 2);
            process_tet_triangle(tet_id, 3);
        };

        dr::parallel_for(dr::blocked_range<size_t>(0, n, 1), [&](dr::blocked_range<size_t> range) {
            for (auto i = range.begin(); i != range.end(); ++i) {
                process_tet(tet_ids[i]);
            }
        });
    }

public:
    bool has_vertex(VertexId vertex_id) const
    {
        VertexKey key{value_of(vertex_id)};
        return m_vertices.has_key(key);
    }

    bool has_tet(TetId tet_id) const
    {
        TetKey key{value_of(tet_id)};
        return m_tets.has_key(key);
    }

    bool has_edge(EdgeId edge_id) const
    {
        TetKey key{value_of(edge_id)};
        return m_tets.has_key(key) && get_edge_index(key) < 6;
    }

    std::span<Scalar, 3> get_vertex(VertexId vertex_id)
    {
        VertexKey key{value_of(vertex_id)};
        auto ptr = m_vertices.get(key);
        if (ptr != nullptr) {
            return *ptr;
        } else {
            throw std::runtime_error("Vertex not found");
        }
    }

    std::span<const Scalar, 3> get_vertex(VertexId vertex_id) const
    {
        VertexKey key{value_of(vertex_id)};
        auto ptr = m_vertices.get(key);
        if (ptr != nullptr) {
            return *ptr;
        } else {
            throw std::runtime_error("Vertex not found");
        }
    }

    std::span<VertexId, 4> get_tet(TetId tet_id)
    {
        TetKey key{value_of(tet_id)};
        auto ptr = m_tets.get(key);
        if (ptr != nullptr) {
            return std::span<VertexId, 4>(ptr->vertices, 4);
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

    std::span<const VertexId, 4> get_tet(TetId tet_id) const
    {
        TetKey key{value_of(tet_id)};
        auto ptr = m_tets.get(key);
        if (ptr != nullptr) {
            return std::span<const VertexId, 4>(ptr->vertices, 4);
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

    std::array<VertexId, 2> get_edge_vertices(EdgeId edge_id) const
    {
        TetKey key{value_of(edge_id)};
        if (!m_tets.has_key(key)) {
            throw std::runtime_error("Edge not found");
        }
        const auto& tet = *m_tets.get(key);
        const auto& local_indices = edge_map[get_edge_index(key)];
        return {tet.vertices[local_indices[0]], tet.vertices[local_indices[1]]};
    }

    TetId get_edge_tet(EdgeId edge_id) const
    {
        TetKey key{value_of(edge_id)};
        key.set_tag(0);
        return TetId(key);
    }

    size_t get_num_vertices() const { return m_vertices.size(); }
    size_t get_num_tets() const { return m_tets.size(); }

    std::tuple<VertexId, EdgeId, EdgeId> split_edge(EdgeId edge_id)
    {
        TetKey key{value_of(edge_id)};
        if (!m_tets.has_key(key)) {
            throw std::runtime_error("Edge not found");
        }
        return split_edge(TetId(key), get_edge_index(key));
    }

    std::tuple<VertexId, EdgeId, EdgeId> split_edge(TetId tet_id, uint8_t local_index)
    {
        if (local_index >= 6) {
            throw std::runtime_error("Invalid local index");
        }


        /**
         * Split the current tet along the edge [llv0, llv1] at vm_id.
         *
         * @param curr_tet_id  The key of the current tet.
         * @param vm_id        The id of the edge mid-point vertex.
         * @param llv0         The local index of the tet.
         * @param llv1         The local index of the tet.
         * @param llv2         The local index of the tet.
         * @param llv3         The local index of the tet.
         *
         * The local indices llv0 and llv1 specifies the edge that we are splitting.
         *
         * @return The ids of the two new tets resulting from the split.
         */
        auto split_tet = [&](TetId curr_tet_id,
                             VertexId vm_id,
                             uint8_t llv0,
                             uint8_t llv1,
                             uint8_t llv2,
                             uint8_t llv3) -> std::array<TetId, 2> {
            TetKey curr_tet_key(value_of(curr_tet_id));
            assert(m_tets.has_key(curr_tet_key));
            const auto& curr_tet = *m_tets.get(curr_tet_key);

            const auto o0_id = curr_tet.mirrors[llv0];
            const auto o1_id = curr_tet.mirrors[llv1];

            // Add two new tets.
            std::array<VertexId, 4> tet_0_vertices{
                curr_tet.vertices[0],
                curr_tet.vertices[1],
                curr_tet.vertices[2],
                curr_tet.vertices[3],
            };
            tet_0_vertices[llv1] = vm_id;
            std::array<VertexId, 4> tet_1_vertices{
                curr_tet.vertices[0],
                curr_tet.vertices[1],
                curr_tet.vertices[2],
                curr_tet.vertices[3],
            };
            tet_1_vertices[llv0] = vm_id;

            auto t0_id =
                add_tet(tet_0_vertices[0], tet_0_vertices[1], tet_0_vertices[2], tet_0_vertices[3]);
            auto t1_id =
                add_tet(tet_1_vertices[0], tet_1_vertices[1], tet_1_vertices[2], tet_1_vertices[3]);

            auto t0_key = TetKey(value_of(t0_id));
            auto t1_key = TetKey(value_of(t1_id));

            auto& tet_0 = *m_tets.get(t0_key);
            auto& tet_1 = *m_tets.get(t1_key);

            // Update mirror information for new tets.
            set_mirror_index(t1_key, llv0, llv1);
            set_mirror_index(t1_key, llv1, llv0);
            set_mirror_index(t1_key, llv2, llv2);
            set_mirror_index(t1_key, llv3, llv3);

            set_mirror_index(t0_key, llv0, llv1);
            set_mirror_index(t0_key, llv1, llv0);
            set_mirror_index(t0_key, llv2, llv2);
            set_mirror_index(t0_key, llv3, llv3);

            tet_0.mirrors[llv0] = TetId(t1_key);
            tet_0.mirrors[llv1] = o1_id;
            tet_0.mirrors[llv2] = TetId(invalid_key);
            tet_0.mirrors[llv3] = TetId(invalid_key);

            tet_1.mirrors[llv0] = o0_id;
            tet_1.mirrors[llv1] = TetId(t0_key);
            tet_1.mirrors[llv2] = TetId(invalid_key);
            tet_1.mirrors[llv3] = TetId(invalid_key);

            // Update mirror information for tet adjacent to the new tets.
            if (value_of(o0_id) != invalid_key) {
                TetKey o0_key = TetKey(value_of(o0_id));
                assert(m_tets.has_key(o0_key));
                auto& tet_o0 = *m_tets.get(o0_key);
                auto o0_llv0 = get_mirror_index(o0_key, llv0);
                auto old_t0_mirror_key = TetKey(value_of(tet_o0.mirrors[o0_llv0]));
                t1_key.set_tag(old_t0_mirror_key.get_tag());
                tet_o0.mirrors[o0_llv0] = TetId(t1_key);
            }
            if (value_of(o1_id) != invalid_key) {
                TetKey o1_key = TetKey(value_of(o1_id));
                assert(m_tets.has_key(o1_key));
                auto& tet_o1 = *m_tets.get(o1_key);
                auto o1_llv1 = get_mirror_index(o1_key, llv1);
                auto old_t1_mirror_key = TetKey(value_of(tet_o1.mirrors[o1_llv1]));
                t0_key.set_tag(old_t1_mirror_key.get_tag());
                tet_o1.mirrors[o1_llv1] = TetId(t0_key);
            }

            return {t0_id, t1_id};
        };

        // Compute local indices.
        // lv0 and lv1 are vertices of the edge with local_index.
        // lv2 and lv3 are the indices of vertices in the opposite edge.
        uint8_t lv0 = edge_map[local_index][0];
        uint8_t lv1 = edge_map[local_index][1];
        uint8_t lv2 = edge_map[local_index][2];
        uint8_t lv3 = edge_map[local_index][3];

        TetKey key{value_of(tet_id)};
        auto ptr = m_tets.get(key);
        if (ptr != nullptr) {
            // Create new vertex
            const auto& curr_tet = *ptr;
            const auto v0_id = curr_tet.vertices[lv0];
            const auto v1_id = curr_tet.vertices[lv1];
            const auto p0 = get_vertex(v0_id);
            const auto p1 = get_vertex(v1_id);
            const auto vm_id =
                add_vertex((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2, (p0[2] + p1[2]) / 2);

            // Compute the staring tet if the edge is on the boundary.
            TetId curr_id = tet_id;
            bool on_boundary = false;
            do {
                TetKey curr_key = TetKey(value_of(curr_id));
                assert(curr_key != invalid_key);
                assert(m_tets.get(curr_key)->vertices[lv0] == v0_id);
                assert(m_tets.get(curr_key)->vertices[lv1] == v1_id);

                TetId next_id;
                uint8_t llv0, llv1, llv2, llv3;
                // Note that we are traversing in the lv2 direction. I.e. the next tet is the tet
                // opposite to the lv2 vertex.
                std::tie(next_id, llv0, llv1, llv3, llv2) =
                    get_next_tet_id(curr_id, lv0, lv1, lv3, lv2);
                if (TetKey(value_of(next_id)) == invalid_key) {
                    on_boundary = true;
                    break;
                }
                curr_id = next_id;
                lv0 = llv0;
                lv1 = llv1;
                lv2 = llv2;
                lv3 = llv3;
            } while (!is_same_tet(curr_id, tet_id));
            // Current key should hold a valid key that serve as the first tet to traverse in the
            // lv3 direction. lv0 to lv3 should be valid local indices for curr_key.
            assert(has_tet(curr_id));
            assert(get_tet(curr_id)[lv0] == v0_id);
            assert(get_tet(curr_id)[lv1] == v1_id);

            // Gather 1-ring tets around the edge
            llvm_vecsmall::SmallVector<TetId, 16> old_one_ring, new_one_ring_0, new_one_ring_1;
            llvm_vecsmall::SmallVector<uint8_t, 16 * 4> local_indices;

            TetId init_id = curr_id;
            do {
                TetKey curr_key = TetKey(value_of(curr_id));
                const auto& curr_tet = *m_tets.get(curr_key);
                assert(curr_tet.vertices[lv0] == v0_id);
                assert(curr_tet.vertices[lv1] == v1_id);

                old_one_ring.push_back(curr_id);
                local_indices.push_back(lv0);
                local_indices.push_back(lv1);
                local_indices.push_back(lv2);
                local_indices.push_back(lv3);
                std::tie(curr_id, lv0, lv1, lv2, lv3) =
                    get_next_tet_id(curr_id, lv0, lv1, lv2, lv3);
            } while (TetKey(value_of(curr_id) != invalid_key && !is_same_tet(curr_id, init_id)));
            size_t one_ring_size = old_one_ring.size();

            // Split 1-ring tets
            for (size_t i = 0; i < one_ring_size; i++) {
                auto curr_id = old_one_ring[i];
                lv0 = local_indices[i * 4 + 0];
                lv1 = local_indices[i * 4 + 1];
                lv2 = local_indices[i * 4 + 2];
                lv3 = local_indices[i * 4 + 3];
                auto [t0_id, t1_id] = split_tet(curr_id, vm_id, lv0, lv1, lv2, lv3);
                new_one_ring_0.push_back(t0_id);
                new_one_ring_1.push_back(t1_id);
            }

            // Update connectivity of the split 1-ring tets
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

                auto curr_tet_id_0 = new_one_ring_0[i];
                auto next_tet_id_0 = new_one_ring_0[j];
                auto curr_tet_key_0 = TetKey(value_of(curr_tet_id_0));
                auto next_tet_key_0 = TetKey(value_of(next_tet_id_0));
                auto& curr_tet_0 = *m_tets.get(curr_tet_key_0);
                auto& next_tet_0 = *m_tets.get(next_tet_key_0);

                assert(curr_tet_0.vertices[curr_lv0] == next_tet_0.vertices[next_lv0]);
                assert(curr_tet_0.vertices[curr_lv1] == next_tet_0.vertices[next_lv1]);
                assert(curr_tet_0.vertices[curr_lv2] == next_tet_0.vertices[next_lv3]);

                set_mirror_index(next_tet_key_0, curr_lv0, next_lv0);
                set_mirror_index(next_tet_key_0, curr_lv1, next_lv1);
                set_mirror_index(next_tet_key_0, curr_lv2, next_lv3);
                set_mirror_index(next_tet_key_0, curr_lv3, next_lv2);
                curr_tet_0.mirrors[curr_lv3] = TetId(next_tet_key_0);

                set_mirror_index(curr_tet_key_0, next_lv0, curr_lv0);
                set_mirror_index(curr_tet_key_0, next_lv1, curr_lv1);
                set_mirror_index(curr_tet_key_0, next_lv2, curr_lv3);
                set_mirror_index(curr_tet_key_0, next_lv3, curr_lv2);
                next_tet_0.mirrors[next_lv2] = TetId(curr_tet_key_0);

                auto curr_tet_id_1 = new_one_ring_1[i];
                auto next_tet_id_1 = new_one_ring_1[j];
                auto curr_tet_key_1 = TetKey(value_of(curr_tet_id_1));
                auto next_tet_key_1 = TetKey(value_of(next_tet_id_1));
                auto& curr_tet_1 = *m_tets.get(curr_tet_key_1);
                auto& next_tet_1 = *m_tets.get(next_tet_key_1);

                assert(curr_tet_1.vertices[curr_lv0] == next_tet_1.vertices[next_lv0]);
                assert(curr_tet_1.vertices[curr_lv1] == next_tet_1.vertices[next_lv1]);
                assert(curr_tet_1.vertices[curr_lv2] == next_tet_1.vertices[next_lv3]);

                set_mirror_index(next_tet_key_1, curr_lv0, next_lv0);
                set_mirror_index(next_tet_key_1, curr_lv1, next_lv1);
                set_mirror_index(next_tet_key_1, curr_lv2, next_lv3);
                set_mirror_index(next_tet_key_1, curr_lv3, next_lv2);
                curr_tet_1.mirrors[curr_lv3] = TetId(next_tet_key_1);

                set_mirror_index(curr_tet_key_1, next_lv0, curr_lv0);
                set_mirror_index(curr_tet_key_1, next_lv1, curr_lv1);
                set_mirror_index(curr_tet_key_1, next_lv2, curr_lv3);
                set_mirror_index(curr_tet_key_1, next_lv3, curr_lv2);
                next_tet_1.mirrors[next_lv2] = TetId(curr_tet_key_1);
            }

            for (auto tid : old_one_ring) {
                m_tets.erase(TetKey(value_of(tid)));
            }

            auto t0_key = TetKey(value_of(new_one_ring_0.front()));
            auto t1_key = TetKey(value_of(new_one_ring_1.front()));
            set_edge_index(t0_key, local_indices[0], local_indices[1]);
            set_edge_index(t1_key, local_indices[0], local_indices[1]);
            return {vm_id, EdgeId(t0_key), EdgeId(t1_key)};
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

    void par_foreach_vertex(
        const std::function<void(VertexId, std::span<const Scalar, 3>)>& callback) const
    {
        const size_t max_valid_index = m_vertices.getMaxValidIndex();
        dr::parallel_for(
            dr::blocked_range<size_t>(0, max_valid_index, 1),
            [&](dr::blocked_range<size_t> range) {
                for (auto index : range) {
                    if (m_vertices.isValidIndex(index)) {
                        VertexMap::const_kv_iterator itr(&m_vertices, index);
                        const auto& [key, value] = *itr;
                        callback(VertexId(key), std::span<const Scalar, 3>(value.get().data(), 3));
                    }
                }
            });
    }

    void seq_foreach_vertex(
        const std::function<void(VertexId, std::span<const Scalar, 3>)>& callback) const
    {
        for (const auto& [key, value] : m_vertices.items()) {
            callback(VertexId(key), std::span<const Scalar, 3>(value.get().data(), 3));
        }
    }

    void par_foreach_tet(
        const std::function<void(TetId, std::span<const VertexId, 4>)>& callback) const
    {
        const size_t max_valid_index = m_tets.getMaxValidIndex();
        dr::parallel_for(
            dr::blocked_range<size_t>(0, max_valid_index, 1),
            [&](dr::blocked_range<size_t> range) {
                for (auto index : range) {
                    if (m_tets.isValidIndex(index)) {
                        TetMap::const_kv_iterator itr(&m_tets, index);
                        const auto& [key, value] = *itr;
                        assert(&value.get() == m_tets.get(key));
                        callback(TetId(key), std::span<const VertexId, 4>(value.get().vertices, 4));
                    }
                }
            });
    }

    void seq_foreach_tet(
        const std::function<void(TetId, std::span<const VertexId, 4>)>& callback) const
    {
        for (const auto& [key, value] : m_tets.items()) {
            assert(&value.get() == m_tets.get(key));
            callback(TetId(key), std::span<const VertexId, 4>(value.get().vertices, 4));
        }
    }

    void foreach_edge_in_tet(
        TetId tet_id,
        const std::function<void(EdgeId, VertexId, VertexId)>& callback)
    {
        TetKey tet_key(value_of(tet_id));
        tet_key.set_tag(0);
        if (!m_tets.has_key(tet_key)) {
            throw std::runtime_error("Tet not found");
        }
        const auto& tet = *m_tets.get(tet_key);
        for (uint8_t i = 0; i < 6; i++) {
            set_edge_index(tet_key, i);
            const auto& lv = edge_map[i];
            callback(EdgeId(tet_key), tet.vertices[lv[0]], tet.vertices[lv[1]]);
        }
    }

    void foreach_tet_around_edge(EdgeId edge_id, const std::function<void(TetId)>& callback) const
    {
        TetKey tet_key(value_of(edge_id));
        if (!m_tets.has_key(tet_key)) {
            throw std::runtime_error("Edge not found");
        }

        auto lvs = edge_map[get_edge_index(tet_key)];
        uint8_t lv0 = lvs[0], lv1 = lvs[1], lv2 = lvs[2], lv3 = lvs[3];

        tet_key.set_tag(0);
        TetId tet_id(tet_key);

        // Compute the staring tet if the edge is on the boundary.
        TetId curr_id = tet_id;
        bool on_boundary = false;
        do {
            TetKey curr_key = TetKey(value_of(curr_id));
            assert(curr_key != invalid_key);

            TetId next_id;
            uint8_t llv0, llv1, llv2, llv3;
            // Note that we are traversing in the lv2 direction. I.e. the next tet is the tet
            // opposite to the lv2 vertex.
            std::tie(next_id, llv0, llv1, llv3, llv2) =
                get_next_tet_id(curr_id, lv0, lv1, lv3, lv2);
            if (TetKey(value_of(next_id)) == invalid_key) {
                on_boundary = true;
                break;
            }
            curr_id = next_id;
            lv0 = llv0;
            lv1 = llv1;
            lv2 = llv2;
            lv3 = llv3;
        } while (!is_same_tet(curr_id, tet_id));
        assert(!is_invalid_tet(curr_id));

        // Iterator over all tets adjacent to the edge and call the callback.
        TetId init_id = curr_id;
        do {
            callback(curr_id);
            std::tie(curr_id, lv0, lv1, lv2, lv3) = get_next_tet_id(curr_id, lv0, lv1, lv2, lv3);
        } while (!is_invalid_tet(curr_id) && !is_same_tet(init_id, curr_id));
    }


public:
    /**
     * Internal: get the adjacent tet of a tet across one of its triangles.
     *
     * @param tet_id       The id of the current tet.
     * @param local_index  The local index of the triangle face of the current tet.
     *
     * @return The id of the adjacent tet.
     */
    TetId get_adjacent_tet(TetId tet_id, uint8_t local_index) const
    {
        TetKey key = TetKey(value_of(tet_id));
        auto ptr = m_tets.get(static_cast<TetKey>(key));
        if (ptr != nullptr) {
            return ptr->mirrors[local_index];
        } else {
            throw std::runtime_error("Tet not found");
        }
    }

    const auto& get_vertices() const { return m_vertices; }
    const auto& get_tets() const { return m_tets; }

private:
    /**
     * Get the next tet adjacent to the current tet around a given edge.
     *
     * @param curr_tet_id  The id of the current tet.
     * @param llv0         The local index of the tet.
     * @param llv1         The local index of the tet.
     * @param llv2         The local index of the tet.
     * @param llv3         The local index of the tet.
     *
     * The local indices llv0 and llv1 specifies the edge that we are traversing around.
     * The local index llv3 specifies the vertex that is opposite to the face shared by the
     * current tet and the next tet.
     *
     * @return The id of the next tet as well as the consistent local indices of the next tet.
     */
    std::tuple<TetId, uint8_t, uint8_t, uint8_t, uint8_t>
    get_next_tet_id(TetId curr_tet_id, uint8_t llv0, uint8_t llv1, uint8_t llv2, uint8_t llv3) const
    {
        TetKey curr_tet_key(value_of(curr_tet_id));
        assert(m_tets.has_key(curr_tet_key));
        assert(llv0 < 4);
        assert(llv1 < 4);
        assert(llv2 < 4);
        assert(llv3 < 4);
        assert(llv0 != llv1);
        assert(llv0 != llv2);
        assert(llv0 != llv3);
        assert(llv1 != llv2);
        assert(llv1 != llv3);
        assert(llv2 != llv3);

        const auto& curr_tet = *m_tets.get(curr_tet_key);
        TetId next_tet_id = curr_tet.mirrors[llv3];
        TetKey next_tet_key = TetKey(value_of(next_tet_id));
        llv0 = get_mirror_index(next_tet_key, llv0);
        llv1 = get_mirror_index(next_tet_key, llv1);
        llv2 = get_mirror_index(next_tet_key, llv2);
        llv3 = get_mirror_index(next_tet_key, llv3);

        // Note that we swap llv2 and llv3 here intentionally.
        return std::make_tuple(next_tet_id, llv0, llv1, llv3, llv2);
    }


private:
    VertexMap m_vertices;
    TetMap m_tets;
};

} // namespace mtet
