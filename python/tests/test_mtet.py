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
