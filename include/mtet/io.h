#pragma once

#include <string>
#include "mtet.h"

namespace mtet {

void save_mesh(std::string filename, const mtet::MTetMesh& mesh);
mtet::MTetMesh load_mesh(std::string filename);

} // namespace mtet

