#pragma once

#include <string>
#include "mtet.h"

namespace mtet {

void save_mesh(std::string filename, const MTetMesh& mesh);
MTetMesh load_mesh(std::string filename);

} // namespace mtet

