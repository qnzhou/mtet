#include <catch2/catch_test_macros.hpp>

#include <mtet/grid.h>

TEST_CASE("grid", "[mtet]")
{
    SECTION("tet5")
    {
        auto grid = mtet::generate_tet_grid({1, 1, 1}, {0, 0, 0}, {1, 1, 1}, mtet::TET5);
        REQUIRE(grid.get_num_vertices() == 8);
        REQUIRE(grid.get_num_tets() == 5);
    }
    SECTION("tet6")
    {
        auto grid = mtet::generate_tet_grid({1, 1, 1}, {0, 0, 0}, {1, 1, 1}, mtet::TET6);
        REQUIRE(grid.get_num_vertices() == 8);
        REQUIRE(grid.get_num_tets() == 6);
    }
}
