#pragma once

#include <span>
#include <string>
#include "mtet.h"

namespace mtet {

void save_mesh(std::string filename, const MTetMesh& mesh);
void save_mesh(std::string filename, const MTetMesh& mesh, std::span<TetId> active_tets);
void save_mesh(
    std::string filename,
    const MTetMesh& mesh,
    std::string scalar_field_name,
    std::span<Scalar> scalar_field);
MTetMesh load_mesh(std::string filename);

} // namespace mtet

