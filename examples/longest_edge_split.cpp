#include <mtet/io.h>
#include <mtet/mtet.h>

#include <algorithm>
#include <vector>

int main(int argc, char** argv)
{
    mtet::MTetMesh mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    mesh.add_tet(v0, v1, v2, v3);

    mtet::save_mesh("init.msh", mesh);

    auto comp = [](std::pair<mtet::Scalar, mtet::EdgeId> e0,
                   std::pair<mtet::Scalar, mtet::EdgeId> e1) { return e0.first < e1.first; };
    std::vector<std::pair<mtet::Scalar, mtet::EdgeId>> Q;
    Q.reserve(1024);

    auto push_longest_edge = [&](mtet::TetId tid) {
        mtet::EdgeId longest_edge;
        mtet::Scalar longest_edge_length = 0;
        mesh.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid, mtet::VertexId v0, mtet::VertexId v1) {
            auto p0 = mesh.get_vertex(v0);
            auto p1 = mesh.get_vertex(v1);
            mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                             (p0[2] - p1[2]) * (p0[2] - p1[2]);
            if (l > longest_edge_length) {
                longest_edge_length = l;
                longest_edge = eid;
            }
        });
        Q.emplace_back(longest_edge_length, longest_edge);
    };

    // Initialize priority queue.
    mesh.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> vs) {
        push_longest_edge(tid);
    });
    std::make_heap(Q.begin(), Q.end(), comp);

    // Keep splitting the longest edge until the longest edge is shorter than 0.01.
    while (!Q.empty()) {
        std::pop_heap(Q.begin(), Q.end(), comp);
        auto [edge_length, eid] = Q.back();
        Q.pop_back();
        if (!mesh.has_edge(eid)) continue;
        if (edge_length < 0.01) break;

        auto [vid, eid0, eid1] = mesh.split_edge(eid);
        mesh.foreach_tet_around_edge(eid0, [&](mtet::TetId tid) {
            push_longest_edge(tid);
            std::push_heap(Q.begin(), Q.end(), comp);
        });
        mesh.foreach_tet_around_edge(eid1, [&](mtet::TetId tid) {
            push_longest_edge(tid);
            std::push_heap(Q.begin(), Q.end(), comp);
        });
    }

    mtet::save_mesh("final.msh", mesh);

    return 0;
}
