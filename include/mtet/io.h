#pragma once

#include <span>
#include <string>
#include "mtet.h"

namespace mtet {

void save_mesh(std::string filename, const MTetMesh& mesh);
void save_mesh(std::string filename, const MTetMesh& mesh, std::span<TetId> active_tets);
MTetMesh load_mesh(std::string filename);

} // namespace mtet

