#include <mtet/io.h>
#include <mtet/mtet.h>

int main(int argc, char** argv)
{
    mtet::MTetMesh mesh;
    auto v0 = mesh.add_vertex(0, 0, 0);
    auto v1 = mesh.add_vertex(1, 0, 0);
    auto v2 = mesh.add_vertex(0, 1, 0);
    auto v3 = mesh.add_vertex(0, 0, 1);
    mesh.add_tet(v0, v1, v2, v3);

    mtet::save_mesh("output.msh", mesh);

    return 0;
}
