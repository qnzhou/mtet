import pytest
import mtet
import numpy as np


class TestMTet:
    def test_simple(self):
        mesh = mtet.MTetMesh()
        v0 = mesh.add_vertex(0, 0, 0)
        v1 = mesh.add_vertex(1, 0, 0)
        assert np.all(mesh.get_vertex(v0) == [0, 0, 0])
        assert np.all(mesh.get_vertex(v1) == [1, 0, 0])

    def test_repeated_access(self):
        mesh = mtet.MTetMesh()
        v0 = mesh.add_vertex(0, 0, 0)
        v1 = mesh.add_vertex(1, 0, 0)
        v2 = mesh.add_vertex(0, 1, 0)
        v3 = mesh.add_vertex(0, 0, 1)
        t0 = mesh.add_tet(v0, v1, v2, v3)
        assert mesh.has_tet(t0)
        assert np.all(mesh.get_tet(t0) == [v0, v1, v2, v3])
        assert np.all(mesh.get_tet(t0) == [v0, v1, v2, v3])

    def test_grid(self):
        grid = mtet.generate_tet_grid([1, 1, 1], style=5)
        assert isinstance(grid, mtet.MTetMesh)
        assert grid.get_num_vertices() == 8
        assert grid.get_num_tets() == 5

        grid = mtet.generate_tet_grid([1, 1, 2], style=5)
        assert grid.get_num_vertices() == 12
        assert grid.get_num_tets() == 10 

        grid = mtet.generate_tet_grid([1, 1, 1], style=6)
        assert grid.get_num_vertices() == 8
        assert grid.get_num_tets() == 6

        grid = mtet.generate_tet_grid([2, 2, 2], style=6)
        assert grid.get_num_vertices() == 27
        assert grid.get_num_tets() == 48

    def test_export(self):
        grid = mtet.generate_tet_grid([2, 2, 2], style=6)
        vertices, tets = grid.export()

        assert vertices.shape == (grid.get_num_vertices(), 3)
        assert tets.shape == (grid.get_num_tets(), 4)
        assert not np.all(vertices == 0)
        assert vertices.min() >= 0.0 and vertices.max() <= 1.0
        assert tets.min() == 0
        assert tets.max() == grid.get_num_vertices() - 1
