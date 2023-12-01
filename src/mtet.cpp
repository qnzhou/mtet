#include <mtet/mtet.h>

#include "MTetMeshImpl.h"

namespace mtet {

MTetMesh::MTetMesh()
    : m_impl(nonstd::make_indirect_value<MTetMeshImpl>())
{}

MTetMesh::~MTetMesh() = default;

MTetMesh::MTetMesh(const MTetMesh&) noexcept = default;
MTetMesh& MTetMesh::operator=(const MTetMesh&) noexcept = default;

MTetMesh::MTetMesh(MTetMesh&&) noexcept = default;
MTetMesh& MTetMesh::operator=(MTetMesh&&) noexcept = default;

uint64_t MTetMesh::add_vertex(Scalar x, Scalar y, Scalar z)
{
    return m_impl->add_vertex(x, y, z);
}

uint64_t MTetMesh::add_tet(uint64_t v0, uint64_t v1, uint64_t v2, uint64_t v3)
{
    return m_impl->add_tet(v0, v1, v2, v3);
}

void MTetMesh::initialize_connectivity()
{
    m_impl->initialize_connectivity();
}

bool MTetMesh::has_vertex(uint64_t vertex_id) const
{
    return m_impl->has_vertex(vertex_id);
}

bool MTetMesh::has_tet(uint64_t tet_id) const
{
    return m_impl->has_tet(tet_id);
}

std::span<Scalar, 3> MTetMesh::get_vertex(uint64_t vertex_id)
{
    return m_impl->get_vertex(vertex_id);
}

std::span<const Scalar, 3> MTetMesh::get_vertex(uint64_t vertex_id) const
{
    return m_impl->get_vertex(vertex_id);
}

std::span<uint64_t, 4> MTetMesh::get_tet(uint64_t tet_id)
{
    return m_impl->get_tet(tet_id);
}

std::span<const uint64_t, 4> MTetMesh::get_tet(uint64_t tet_id) const
{
    return m_impl->get_tet(tet_id);
}

size_t MTetMesh::get_num_vertices() const
{
    return m_impl->get_num_vertices();
}

size_t MTetMesh::get_num_tets() const
{
    return m_impl->get_num_tets();
}

uint64_t MTetMesh::split_edge(uint64_t tet_id, uint8_t local_index)
{
    return m_impl->split_edge(tet_id, local_index);
}

void MTetMesh::par_foreach_vertex(
    const std::function<void(uint64_t, std::span<const Scalar, 3>)>& func) const
{
    m_impl->par_foreach_vertex(func);
}

void MTetMesh::seq_foreach_vertex(
    const std::function<void(uint64_t, std::span<const Scalar, 3>)>& func) const
{
    m_impl->seq_foreach_vertex(func);
}

void MTetMesh::par_foreach_tet(
    const std::function<void(uint64_t, std::span<const uint64_t, 4>)>& func) const
{
    m_impl->par_foreach_tet(func);
}

void MTetMesh::seq_foreach_tet(
    const std::function<void(uint64_t, std::span<const uint64_t, 4>)>& func) const
{
    m_impl->seq_foreach_tet(func);
}

} // namespace mtet
