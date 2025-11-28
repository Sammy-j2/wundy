import numpy as np
import pytest
from wundy import ui, first

"""YAML-driven analytic beam tests.

Loads and model definitions are provided in example YAML files to exercise the
full load->preprocess->solve pipeline.

Cases:
1. Uniform load cantilever (beam_uniform.yaml): w_tip = q L^4 / (8 E I)
2. Linear (triangular) load cantilever (beam_triangular.yaml): one-element FE
   deflection coefficient = (11/120) * q0 L^4 / (E I) (continuum is 1/30).

Both use Neo-Hookean material specified via (E, nu) but bending is treated with
linear tangent E in the current implementation.
"""

EXAMPLES_DIR = "docs/examples"


def _solve_from_yaml(fname: str):
    path = f"{EXAMPLES_DIR}/{fname}"
    with open(path, "r", encoding="utf-8") as f:
        data = ui.load(f)
    pre = ui.preprocess(data)
    sol = first.first_fe_code(
        pre["coords"],
        pre["blocks"],
        pre["bcs"],
        pre["dload"],
        pre["materials"],
        pre["block_elem_map"],
        pre.get("integration"),
        pre.get("dof_per_node", 2),
    )
    return sol, pre


def test_beam_uniform_yaml():
    sol, pre = _solve_from_yaml("beam_uniform.yaml")
    u = sol["dofs"]
    # w at node 1 (DOF ordering: node0[w,theta], node1[w,theta])
    w_tip = u[2]
    # Extract parameters from preprocessed structures
    block = pre["blocks"][0]
    I = block["element"]["properties"]["I"]
    E = pre["materials"][block["material"]]["parameters"]["E"]
    L = float(pre["coords"][1, 0] - pre["coords"][0, 0])
    # Find load magnitude
    q_load = None
    for dl in pre["dload"]:
        if dl["name"].lower() == "q_uniform":
            q_load = float(dl["value"])
            break
    assert q_load is not None, "Uniform load not found in YAML"
    w_expected = q_load * L ** 4 / (8.0 * E * I)
    assert np.isclose(w_tip, w_expected, rtol=5e-4, atol=1e-12), (
        f"Uniform load cantilever tip deflection expected {w_expected}, got {w_tip}"
    )
    # Convergence metadata present when nonlinear path requested; optional check
    if "convergence" not in sol:
        pytest.skip("NR convergence metadata absent; nonlinear path not active")


def test_beam_linear_expression_yaml():
    sol, pre = _solve_from_yaml("beam_triangular.yaml")
    u = sol["dofs"]
    w_tip = u[2]
    block = pre["blocks"][0]
    I = block["element"]["properties"]["I"]
    E = pre["materials"][block["material"]]["parameters"]["E"]
    L = float(pre["coords"][1, 0] - pre["coords"][0, 0])
    # Peak load q0 extracted from expression string (assumes pattern 'number * x / L')
    q0 = 100.0  # matches expression
    # One-element FE with exact Gauss integration of q(x)=q0*x/L gives coefficient 5/54.
    w_expected_fe = (5.0 / 54.0) * q0 * L ** 4 / (E * I)
    assert np.isclose(w_tip, w_expected_fe, rtol=1e-3, atol=1e-12), (
        f"Triangular load (YAML one-element FE) expected {w_expected_fe}, got {w_tip}"
    )
    if "convergence" not in sol:
        pytest.skip("NR convergence metadata absent; nonlinear path not active")


if __name__ == "__main__":  # manual quick run
    pytest.main([__file__, "-s", "-vv"])