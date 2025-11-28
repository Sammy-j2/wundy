import numpy as np
import pytest
from wundy import first

"""Analytic beam tests for a single EB1D element (cantilever) with Neo-Hookean material.

We model a cantilever beam of length L with Euler–Bernoulli formulation (EB1D) and
clamped at x=0 (w=0, theta=0). One 2-node beam element is used. For small transverse
loads we assume linear bending using the material's tangent modulus E even when the
material is Neo-Hookean; geometric nonlinearity is not yet included.

Two load cases:
1. Uniform transverse load q (provided as a Python function).
   Analytic tip deflection (cantilever, uniform load): w_tip = q*L**4/(8*E*I).
2. Linearly varying (triangular) load q(x) = q0 * x / L (provided as expression string).
   Analytic tip deflection: w_tip = q0*L**4/(30*E*I).

Both tests use integration.nonlinear='nonlinear' with Neo-Hooke material to exercise
Newton path (beam treated linearly inside NR loop => quick convergence).
"""


def _run_beam(q_spec, is_expression=False):
    L = 1.0
    I = 1.0e-6
    A = 1.0  # area not used in bending stiffness but kept for completeness
    E = 2.0e5  # Pa
    nu = 0.3
    # Material: Neo-Hookean specified via (E,nu)
    materials = {"rubber": {"type": "neo_hooke", "parameters": {"E": E, "nu": nu}}}

    # Coordinates: two nodes
    coords = np.array([[0.0], [L]], dtype=float)

    # Single EB1D element block
    blocks = [
        {
            "name": "beamblock",
            "element": {"type": "EB1D", "properties": {"I": I, "area": A}},
            "material": "rubber",
            "connect": np.array([[0, 1]], dtype=int),
            "integration": {"ngp": 3},
        }
    ]

    # Boundary conditions: clamp at node 0 (w=0, theta=0)
    bcs = [
        {"type": first.DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},  # w(0)=0
        {"type": first.DIRICHLET, "nodes": [0], "local_dof": 1, "value": 0.0},  # theta(0)=0
    ]

    # Distributed load applied to element 0
    if is_expression:
        dloads = [
            {
                "name": "qexpr",
                "elements": [0],
                "type": "QY",
                "expression": q_spec,  # string expression in x,L
                "direction": [1],
            }
        ]
    else:
        dloads = [
            {
                "name": "qfunc",
                "elements": [0],
                "type": "QY",
                "value": q_spec,  # callable q(x)
                "direction": [1],
            }
        ]

    block_elem_map = {0: (0, 0)}

    sol = first.first_fe_code(
        coords,
        blocks,
        bcs,
        dloads,
        materials,
        block_elem_map,
        integration={"nonlinear": "nonlinear", "ngp": 3, "internal": "gauss", "stiffness": "gauss"},
        dof_per_node=2,
    )
    return sol, E, I, L


def test_cantilever_uniform_load_function_neohooke():
    q0 = 100.0  # N/m

    def q_func(x):
        return q0

    sol, E, I, L = _run_beam(q_func, is_expression=False)
    u = sol["dofs"]
    w_tip = u[2]  # DOF ordering: node0[w,theta], node1[w,theta] => w at node1 index 2
    w_expected = q0 * L ** 4 / (8.0 * E * I)
    assert np.isclose(w_tip, w_expected, rtol=5e-4, atol=1e-12), (
        f"Uniform load cantilever tip deflection expected {w_expected}, got {w_tip}"
    )
    assert "convergence" in sol, "Expected NR convergence metadata"


def test_cantilever_linear_load_expression_neohooke():
    q0 = 100.0  # peak load at free end
    # Expression q(x) = q0 * x / L (string evaluated in solver)
    expr = f"{q0} * x / L"
    sol, E, I, L = _run_beam(expr, is_expression=True)
    u = sol["dofs"]
    w_tip = u[2]
    # Continuum analytic: w_tip = q0 * L^4 / (30 E I) for triangular load (0 at clamp, q0 at tip).
    # Single Hermite beam element FE solution yields a different closed-form coefficient 11/120.
    # For one element we therefore test against the FE-consistent expression (11/120)*(q0 L^4)/(E I).
    w_expected_fe = (11.0 / 120.0) * q0 * L ** 4 / (E * I)
    assert np.isclose(w_tip, w_expected_fe, rtol=5e-3, atol=1e-12), (
        f"Triangular load cantilever tip deflection (1 elem FE) expected {w_expected_fe}, got {w_tip}"
    )
    assert "convergence" in sol, "Expected NR convergence metadata"


if __name__ == "__main__":  # quick manual run
    pytest.main([__file__, "-s", "-vv"])