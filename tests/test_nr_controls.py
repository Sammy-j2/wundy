import numpy as np
from wundy import first


def test_nr_du_tol_triggers_convergence():
    # Simple 2-node bar with Neo-Hookean material; use a very large nr_du_tol
    coords = np.array([[0.0], [1.0]])
    # one block containing one element connecting nodes 0 and 1
    blocks = [
        {
            "element": {"properties": {"area": 1.0}},
            "material": "RUBBER",
            "connect": [(0, 1)],
        }
    ]
    # Neo-Hookean material parameters
    materials = {
        "RUBBER": {
            "type": "neo_hooke",
            "parameters": {"mu": 3.0e5, "lambda": 2.0e5},
            "density": 1200.0,
        }
    }
    # Boundary conditions: fixed at node 0, Neumann at node 1
    bcs = [
        {"type": first.DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},
        {"type": first.NEUMANN, "nodes": [1], "local_dof": 0, "value": 1000.0},
    ]
    dloads = []
    # mapping: global element id 0 -> (block_index=0, local_elem_index=0)
    block_elem_map = {0: (0, 0)}

    # Request nonlinear solve and set nr_du_tol to a very large number so the
    # Newton loop should stop because du_norm < nr_du_tol
    integration = {"nonlinear": "nonlinear", "nr_du_tol": 1e9}

    sol = first.first_fe_code(coords, blocks, bcs, dloads, materials, block_elem_map, integration=integration)

    # Check that we got a convergence dict and that reason is 'du_tol'
    assert "convergence" in sol, "Solution should include convergence info"
    assert sol["convergence"].get("reason") == "du_tol", f"Expected du_tol convergence, got {sol['convergence']}"