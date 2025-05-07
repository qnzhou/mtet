#pragma once

#include <mtet/mtet.h>

#include <array>

namespace mtet {

enum GridStyle {
    TET5, /// 5 tetrahedrons per grid cell
    TET6 /// 6 tetrahedrons per grid cell
};

/**
 * @brief Generate a tetrahedral grid with the specified resolution and bounding box.
 *
 * @param resolution The number of divisions along each axis (x, y, z).
 * @param bbox_min The minimum coordinates of the bounding box.
 * @param bbox_max The maximum coordinates of the bounding box.
 * @param style The style of the tetrahedral mesh (TET5 or TET6).
 *
 * @return A tetrahedral mesh object.
 */
mtet::MTetMesh generate_tet_grid(
    const std::array<size_t, 3>& resolution,
    const std::array<float, 3>& bbox_min,
    const std::array<float, 3>& bbox_max,
    GridStyle style = TET5);

} // namespace mtet
