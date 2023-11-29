#pragma once

#include <functional>
#include <span>

#include <nonstd/indirect_value.hpp>

namespace mtet {

using Scalar = double;
class MTetMeshImpl;

class MTetMesh
{
public:
    MTetMesh();
    ~MTetMesh();

    MTetMesh(const MTetMesh&) noexcept;
    MTetMesh& operator=(const MTetMesh&) noexcept;

    MTetMesh(MTetMesh&&) noexcept;
    MTetMesh& operator=(MTetMesh&&) noexcept;

public:
    uint64_t add_vertex(Scalar x, Scalar y, Scalar z);
    uint64_t add_tet(uint64_t v0, uint64_t v1, uint64_t v2, uint64_t v3);
    void initialize_connectivity();

public:
    bool has_vertex(uint64_t vertex_id) const;
    bool has_tet(uint64_t tet_id) const;

    std::span<Scalar, 3> get_vertex(uint64_t vertex_id);
    std::span<uint64_t, 4> get_tet(uint64_t tet_id);

    std::span<const Scalar, 3> get_vertex(uint64_t vertex_id) const;
    std::span<const uint64_t, 4> get_tet(uint64_t tet_id) const;

    size_t get_num_vertices() const;
    size_t get_num_tets() const;

public:
    /**
     * Split the edge of the given tet with the given local edge id.
     *
     * The split will insert a new vertex at the mid point of the specified edge.
     *
     * @param tet_id        The id of the tet to split.
     * @param local_edge_id The local edge id of the edge to split.
     *
     * @return The id of the new vertex.
     *
     * Let the oriented tet be [v0, v1, v2, v3]. The local edges order are the following:
     *   0: [v0, v1],
     *   1: [v1, v2],
     *   2: [v2, v0],
     *   3: [v0, v3],
     *   4: [v1, v3],
     *   5: [v2, v3]
     *
     * The indices are illustrated in the ASCII drawing below.
     *
     *               v2
     *             ,/|`\
     *           ,/  |  `\
     *         ,2    |.   `1
     *       ,/       5     `\
     *     ,/         |       `\
     *   v0-------0-- |. -------v1
     *     `\.         |      ,/
     *        `\.      |    ,4
     *           `3.   |. ,/
     *              `\. |/
     *                 `v3
     */
    uint64_t split_edge(uint64_t tet_id, uint8_t local_edge_id);

    void par_foreach_vertex(const std::function<void(uint64_t)>& callback) const;
    void seq_foreach_vertex(const std::function<void(uint64_t)>& callback) const;
    void par_foreach_tet(const std::function<void(uint64_t)>& callback) const;
    void seq_foreach_tet(const std::function<void(uint64_t)>& callback) const;

private:
    nonstd::indirect_value<MTetMeshImpl> m_impl;
};


} // namespace mtet
