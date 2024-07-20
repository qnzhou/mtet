#include <mtet/mtet.h>

#include "MTetMeshImpl.h"

namespace mtet {

bool operator==(TetId t0, TetId t1)
{
    return is_same_tet(t0, t1);
}

MTetMesh::MTetMesh()
    : m_impl(nonstd::make_indirect_value<MTetMeshImpl>())
{}

MTetMesh::~MTetMesh() = default;

MTetMesh::MTetMesh(const MTetMesh&) noexcept = default;
MTetMesh& MTetMesh::operator=(const MTetMesh&) noexcept = default;

MTetMesh::MTetMesh(MTetMesh&&) noexcept = default;
MTetMesh& MTetMesh::operator=(MTetMesh&&) noexcept = default;

VertexId MTetMesh::add_vertex(Scalar x, Scalar y, Scalar z)
{
    return m_impl->add_vertex(x, y, z);
}

TetId MTetMesh::add_tet(VertexId v0, VertexId v1, VertexId v2, VertexId v3)
{
    return m_impl->add_tet(v0, v1, v2, v3);
}

void MTetMesh::initialize_connectivity()
{
    m_impl->initialize_connectivity();
}

bool MTetMesh::has_vertex(VertexId vertex_id) const
{
    return m_impl->has_vertex(vertex_id);
}

bool MTetMesh::has_tet(TetId tet_id) const
{
    return m_impl->has_tet(tet_id);
}

bool MTetMesh::has_edge(EdgeId edge_id) const
{
    return m_impl->has_edge(edge_id);
}

std::span<Scalar, 3> MTetMesh::get_vertex(VertexId vertex_id)
{
    return m_impl->get_vertex(vertex_id);
}

std::span<const Scalar, 3> MTetMesh::get_vertex(VertexId vertex_id) const
{
    return m_impl->get_vertex(vertex_id);
}

std::span<VertexId, 4> MTetMesh::get_tet(TetId tet_id)
{
    return m_impl->get_tet(tet_id);
}

std::span<const VertexId, 4> MTetMesh::get_tet(TetId tet_id) const
{
    return m_impl->get_tet(tet_id);
}

EdgeId MTetMesh::get_edge(TetId tet_id, uint8_t local_index) const
{
    return m_impl->get_edge(tet_id, local_index);
}

std::array<VertexId, 2> MTetMesh::get_edge_vertices(EdgeId edge_id) const
{
    return m_impl->get_edge_vertices(edge_id);
}

TetId MTetMesh::get_edge_tet(EdgeId edge_id) const
{
    return m_impl->get_edge_tet(edge_id);
}

TetId MTetMesh::get_mirror(TetId tet_id, uint8_t local_index) const
{
    return m_impl->get_mirror(tet_id, local_index);
}

size_t MTetMesh::get_num_vertices() const
{
    return m_impl->get_num_vertices();
}

size_t MTetMesh::get_num_tets() const
{
    return m_impl->get_num_tets();
}

std::tuple<VertexId, EdgeId, EdgeId> MTetMesh::split_edge(EdgeId edge_id)
{
    return m_impl->split_edge(edge_id);
}

std::tuple<VertexId, EdgeId, EdgeId> MTetMesh::split_edge(TetId tet_id, uint8_t local_index)
{
    return m_impl->split_edge(tet_id, local_index);
}

void MTetMesh::par_foreach_vertex(
    const std::function<void(VertexId, std::span<const Scalar, 3>)>& func) const
{
    m_impl->par_foreach_vertex(func);
}

void MTetMesh::seq_foreach_vertex(
    const std::function<void(VertexId, std::span<const Scalar, 3>)>& func) const
{
    m_impl->seq_foreach_vertex(func);
}

void MTetMesh::par_foreach_tet(
    const std::function<void(TetId, std::span<const VertexId, 4>)>& func) const
{
    m_impl->par_foreach_tet(func);
}

void MTetMesh::seq_foreach_tet(
    const std::function<void(TetId, std::span<const VertexId, 4>)>& func) const
{
    m_impl->seq_foreach_tet(func);
}

void MTetMesh::foreach_edge_in_tet(
    TetId tet_id,
    const std::function<void(EdgeId, VertexId, VertexId)>& callback)
{
    m_impl->foreach_edge_in_tet(tet_id, callback);
}

void MTetMesh::foreach_tet_around_edge(EdgeId edge_id, const std::function<void(TetId)>& callback)
    const
{
    m_impl->foreach_tet_around_edge(edge_id, callback);
}

} // namespace mtet
