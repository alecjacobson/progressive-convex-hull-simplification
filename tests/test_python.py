"""Exercise the installed extension, including its native library dependencies."""

import numpy as np
import pytest

import pchs


def assert_polygon_mesh(vertices, indices, offsets):
    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert np.isfinite(vertices).all()
    assert offsets[0] == 0 and offsets[-1] == len(indices)
    assert (np.diff(offsets) >= 3).all()
    assert ((indices >= 0) & (indices < len(vertices))).all()


def test_halfspace_cube():
    normals = np.vstack((np.eye(3), -np.eye(3)))
    halfspaces = np.column_stack((normals, -np.ones(6)))
    center = pchs.chebyshev_center(halfspaces)
    np.testing.assert_allclose(center, np.zeros(3), atol=1e-10)

    for mesh in (
        pchs.primal_mesh_from_halfspaces(halfspaces),
        pchs.primal_mesh_from_halfspaces(halfspaces, center),
    ):
        assert_polygon_mesh(*mesh)
        vertices, _, offsets = mesh
        assert vertices.shape == (8, 3)
        assert len(offsets) == 7
        np.testing.assert_allclose(np.abs(vertices), np.ones((8, 3)), atol=1e-10)


@pytest.mark.parametrize("cost", [pchs.CostFunction.volume, pchs.CostFunction.area])
def test_simplify_octahedron(cost):
    vertices = np.array([
        [1., 0., 0.], [-1., 0., 0.], [0., 1., 0.],
        [0., -1., 0.], [0., 0., 1.], [0., 0., -1.],
    ])
    faces = np.array([
        [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
        [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5],
    ], dtype=np.int32)
    hull = pchs.ConvexHullSimplification(vertices, faces, cost_function=cost)
    assert hull.num_dual_vertices() == 8
    hull.simplify_to(6)
    assert hull.num_dual_vertices() == 6
    assert len(hull.popped_dual_vertex_ids()) == 2
    assert_polygon_mesh(*hull.get_primal_mesh())
    dual_vertices, dual_faces = hull.get_dual_mesh()
    assert dual_vertices.shape == (6, 3)
    assert dual_faces.shape[1] == 3
    assert np.isfinite(dual_vertices).all()
    assert_polygon_mesh(*pchs.simplify_convex_hull(
        vertices, faces, 6, cost_function=cost
    ))
