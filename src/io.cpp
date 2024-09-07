#include <mtet/io.h>

#include <ankerl/unordered_dense.h>
#include <mshio/mshio.h>

#include <cassert>

namespace mtet {

namespace internal {
mshio::MshSpec generate_spec(const MTetMesh& mesh)
{
    mshio::MshSpec spec;
    spec.mesh_format.file_type = 1; // binary

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

    mesh.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data) {
        size_t vertex_tag = vertex_tag_map.size() + 1;
        vertex_tag_map[value_of(vid)] = vertex_tag;
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
    mesh.seq_foreach_tet([&](TetId, std::span<const VertexId, 4> data) {
        element_block.data.insert(
            element_block.data.end(),
            {tet_tag,
             vertex_tag_map[value_of(data[0])],
             vertex_tag_map[value_of(data[1])],
             vertex_tag_map[value_of(data[2])],
             vertex_tag_map[value_of(data[3])]});
        tet_tag++;
    });

    return spec;
}

void add_scalar_field(mshio::MshSpec& spec, std::string name, std::span<Scalar> field)
{
    mshio::Data node_data;
    node_data.header.string_tags.push_back(name);
    node_data.header.int_tags.push_back(0);
    node_data.header.int_tags.push_back(1);
    node_data.header.int_tags.push_back(field.size());

    auto& entries = node_data.entries;
    entries.resize(field.size());
    for (size_t i=0; i<field.size(); i++) {
        auto& entry = entries[i];
        entry.tag = i+1;
        entry.data = {field[i]};
    }

    spec.node_data.push_back(std::move(node_data));
}

} // namespace internal

void save_mesh(std::string filename, const MTetMesh& mesh)
{
    mshio::MshSpec spec = internal::generate_spec(mesh);
    mshio::save_msh(filename, spec);
}

void save_mesh(std::string filename, const MTetMesh& mesh, std::span<TetId> active_tets)
{
    mshio::MshSpec spec;
    spec.mesh_format.file_type = 1; // binary

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

    mesh.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data) {
        size_t vertex_tag = vertex_tag_map.size() + 1;
        vertex_tag_map[value_of(vid)] = vertex_tag;
        node_block.tags.push_back(vertex_tag);
        node_block.data.insert(node_block.data.end(), data.begin(), data.end());
    });

    auto& elements = spec.elements;
    elements.num_entity_blocks = 1;
    elements.num_elements = active_tets.size();
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
    for (auto tet_id : active_tets) {
        assert(mesh.has_tet(tet_id));
        auto data = mesh.get_tet(tet_id);
        element_block.data.insert(
            element_block.data.end(),
            {tet_tag,
             vertex_tag_map[value_of(data[0])],
             vertex_tag_map[value_of(data[1])],
             vertex_tag_map[value_of(data[2])],
             vertex_tag_map[value_of(data[3])]});
        tet_tag++;
    }

    mshio::save_msh(filename, spec);
}

void save_mesh(
    std::string filename,
    const MTetMesh& mesh,
    std::string scalar_field_name,
    std::span<Scalar> scalar_field)
{
    mshio::MshSpec spec = internal::generate_spec(mesh);
    internal::add_scalar_field(spec, scalar_field_name, scalar_field);
    mshio::save_msh(filename, spec);
}

MTetMesh load_mesh(std::string filename)
{
    mshio::MshSpec spec = mshio::load_msh(filename);
    std::vector<VertexId> vertex_ids;
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
            mesh.add_tet(
                vertex_ids[element_block.data[i + 1] - 1],
                vertex_ids[element_block.data[i + 2] - 1],
                vertex_ids[element_block.data[i + 3] - 1],
                vertex_ids[element_block.data[i + 4] - 1]);
        }
    }

    mesh.initialize_connectivity();
    return mesh;
}

} // namespace mtet
