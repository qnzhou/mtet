#include <mtet/io.h>

#include <ankerl/unordered_dense.h>
#include <mshio/mshio.h>

namespace mtet {
void save_mesh(std::string filename, const mtet::MTetMesh& mesh)
{
    mshio::MshSpec spec;
    spec.mesh_format.file_type = 0; // binary

    auto& nodes = spec.nodes;
    nodes.num_entity_blocks = 1;
    nodes.num_nodes = mesh.get_num_vertices();
    nodes.min_node_tag = 1;
    nodes.max_node_tag = nodes.num_nodes;
    nodes.entity_blocks.resize(1);

    auto& node_block = nodes.entity_blocks[0];
    node_block.entity_dim = 3;
    node_block.entity_tag = 1;
    node_block.parametric = 0;
    node_block.num_nodes_in_block = nodes.num_nodes;
    node_block.tags.reserve(nodes.num_nodes);
    node_block.data.reserve(nodes.num_nodes * 3);

    using IndexMap = ankerl::unordered_dense::map<uint64_t, size_t>;
    IndexMap vertex_tag_map;
    vertex_tag_map.reserve(mesh.get_num_vertices());

    mesh.seq_foreach_vertex([&](uint64_t vid, std::span<const mtet::Scalar, 3> data) {
        size_t vertex_tag = vertex_tag_map.size() + 1;
        vertex_tag_map[vid] = vertex_tag;
        node_block.tags.push_back(vertex_tag);
        node_block.data.insert(node_block.data.end(), data.begin(), data.end());
    });

    auto& elements = spec.elements;
    elements.num_entity_blocks = 1;
    elements.num_elements = mesh.get_num_tets();
    elements.min_element_tag = 1;
    elements.max_element_tag = elements.num_elements;
    elements.entity_blocks.resize(1);

    auto& element_block = elements.entity_blocks[0];
    element_block.entity_dim = 3;
    element_block.entity_tag = 1;
    element_block.element_type = 4;
    element_block.num_elements_in_block = elements.num_elements;
    element_block.data.reserve(elements.num_elements * 5);

    size_t tet_tag = 1;
    mesh.seq_foreach_tet([&](uint64_t tet_id, std::span<const uint64_t, 4> data) {
        element_block.data.insert(
            element_block.data.end(),
            {tet_tag,
             vertex_tag_map[data[0]],
             vertex_tag_map[data[1]],
             vertex_tag_map[data[2]],
             vertex_tag_map[data[3]]});
        tet_tag++;
    });

    mshio::save_msh(filename, spec);
}

} // namespace mtet
