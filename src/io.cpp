#include <mtet/io.h>

#include <ankerl/unordered_dense.h>
#include <mshio/mshio.h>

#include <cassert>

namespace mtet {
void save_mesh(std::string filename, const MTetMesh& mesh)
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

    mesh.seq_foreach_vertex([&](uint64_t vid, std::span<const Scalar, 3> data) {
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

MTetMesh load_mesh(std::string filename)
{
    mshio::MshSpec spec = mshio::load_msh(filename);
    std::vector<uint64_t> vertex_ids;
    vertex_ids.reserve(spec.nodes.num_nodes);

    MTetMesh mesh;
    for (const auto& node_block : spec.nodes.entity_blocks) {
        assert(node_block.entity_dim == 3);
        assert(node_block.data.size() % 3 == 0);
        for (size_t i = 0; i < node_block.data.size(); i += 3) {
            auto vid =
                mesh.add_vertex(node_block.data[i], node_block.data[i + 1], node_block.data[i + 2]);
            vertex_ids.push_back(vid);
        }
    }

    for (const auto& element_block : spec.elements.entity_blocks) {
        assert(element_block.entity_dim == 3);
        assert(element_block.element_type == 4);
        assert(element_block.data.size() % 5 == 0);
        for (size_t i = 0; i < element_block.data.size(); i += 5) {
            mesh.add_tet(vertex_ids[element_block.data[i + 1] - 1],
                         vertex_ids[element_block.data[i + 2] - 1],
                         vertex_ids[element_block.data[i + 3] - 1],
                         vertex_ids[element_block.data[i + 4] - 1]);
        }
    }

    mesh.initialize_connectivity();
    return mesh;
}

} // namespace mtet