import numpy as np

from wundy import first
from wundy.schemas import DIRICHLET, NEUMANN


def test_multidof_axial_equivalence():
    # simple 2-node, 1-element problem with axial DOF only active (local index 0)
    coords = np.array([[0.0], [1.0]])
    blocks = [
        {
            "element": {"properties": {"area": 1.0}},
            "material": "mat",
            "connect": [(0, 1)],
        }
    ]
    materials = {"mat": {"type": "ELASTIC", "parameters": {"E": 10.0}}}
    bcs = [
        {"type": DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},
        {"type": NEUMANN, "nodes": [1], "local_dof": 0, "value": 100.0},
    ]
    block_elem_map = {0: (0, 0)}

    sol1 = first.first_fe_code(coords, blocks, bcs, [], materials, block_elem_map, dof_per_node=1)
    sol2 = first.first_fe_code(coords, blocks, bcs, [], materials, block_elem_map, dof_per_node=2)

    u1 = sol1["dofs"]
    u2 = sol2["dofs"]

    # shapes
    assert u1.shape[0] == 2
    assert u2.shape[0] == 4

    # axial DOFs in the expanded vector should match the single-DOF solution
    assert np.allclose(u2[0], u1[0])
    assert np.allclose(u2[2], u1[1])

    # non-axial DOFs remain zero (uncoupled)
    assert u2[1] == 0.0
    assert u2[3] == 0.0

    # check expanded stiffness has axial entries in the expected locations
    K2 = sol2["stiff"]
    assert np.isclose(K2[0, 0], 10.0)
    assert np.isclose(K2[0, 2], -10.0)
