import numpy as np
from wundy import first


def test_single_bar_matches_analytic():
    """
    Analytic check for a single 1D bar element:

    For a bar of length L, area A, Young's modulus E, with u(0)=0 and a
    concentrated axial force F applied at x=L, the displacement at x=L is

        u(L) = F * L / (A * E)

    and the (constant) axial strain is

        eps = (u(L) - u(0)) / L = F / (A * E).
    """
    # Problem parameters
    L = 1.0
    A = 2.0
    E = 1000.0
    F = 100.0

    # Mesh: two nodes, one element
    coords = np.array([[0.0], [L]], dtype=float)
    blocks = [
        {
            "element": {"properties": {"area": A}},
            "material": "steel",
            "connect": [(0, 1)],
            # No per-block integration override needed
        }
    ]
    materials = {"steel": {"type": "elastic", "parameters": {"E": E}}}

    # Boundary conditions: u(0)=0 (Dirichlet), traction at x=L (Neumann)
    bcs = [
        {"type": first.DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},
        {"type": first.NEUMANN, "nodes": [1], "local_dof": 0, "value": F},
    ]
    dloads = []

    # Element id mapping for loads (global element id 0 maps to block 0, local 0)
    block_elem_map = {0: (0, 0)}

    # Solve (linear path)
    sol = first.first_fe_code(
        coords,
        blocks,
        bcs,
        dloads,
        materials,
        block_elem_map,
        integration={"stiffness": "analytic", "internal": "analytic", "ngp": 2, "nonlinear": "linearize"},
        dof_per_node=1,
    )

    # Extract displacement and check analytic value at x=L
    u = sol["dofs"]
    uL_expected = F * L / (A * E)
    # Print displacements, strain, and stress for analytic comparison
    print(f"[Elastic] Displacements u = {u}")
    assert np.isclose(u[1], uL_expected, rtol=1e-12, atol=1e-12), (
        f"u(L) expected {uL_expected}, got {u[1]}"
    )

    # Check strain via helper
    eps = first.element_strain(L, np.array([u[0], u[1]]))
    eps_expected = F / (A * E)
    sigma = E * eps
    print(f"[Elastic] Axial strain eps = {eps}")
    print(f"[Elastic] Axial stress sigma = {sigma}")
    assert np.isclose(eps, eps_expected, rtol=1e-12, atol=1e-12), (
        f"strain expected {eps_expected}, got {eps}"
    )
    # Analytical stress check: sigma = E * eps = F / A
    sigma_expected = F / A
    assert np.isclose(sigma, sigma_expected, rtol=1e-12, atol=1e-12), (
        f"stress expected {sigma_expected}, got {sigma}"
    )

    # Optional: residual inspection. The returned residual is external minus internal
    # over all DOFs (including prescribed); equilibrium is satisfied on free DOFs.
    # We focus test on analytic displacement/strain; skip strict residual norm.
    # R = sol.get("residual")
    # assert R is not None


def test_single_bar_neohook_small_load_matches_linear_tangent():
    """
    For small loads, the Neo-Hookean NR solution should match the linear tangent
    prediction u(L) ≈ F*L/(A*E_equiv).

    Use Neo-Hookean specified via (E, nu) so material_tangent_modulus returns E.
    """
    L = 1.0
    A = 2.0
    E = 1000.0
    nu = 0.3
    F = 1.0  # small load to stay near small-strain regime

    coords = np.array([[0.0], [L]], dtype=float)
    blocks = [
        {
            "element": {"properties": {"area": A}},
            "material": "rubber",
            "connect": [(0, 1)],
        }
    ]
    materials = {"rubber": {"type": "neo_hooke", "parameters": {"E": E, "nu": nu}}}
    bcs = [
        {"type": first.DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},
        {"type": first.NEUMANN, "nodes": [1], "local_dof": 0, "value": F},
    ]
    dloads = []
    block_elem_map = {0: (0, 0)}

    sol = first.first_fe_code(
        coords,
        blocks,
        bcs,
        dloads,
        materials,
        block_elem_map,
        integration={"stiffness": "analytic", "internal": "analytic", "ngp": 2, "nonlinear": "nonlinear"},
        dof_per_node=1,
    )

    u = sol["dofs"]
    # Expected displacement using PK1 small-strain tangent: E_eff = 2*mu + kappa
    mu = E / (2.0 * (1.0 + nu))
    kappa = E / (3.0 * (1.0 - 2.0 * nu))
    E_eff = 2.0 * mu + kappa
    uL_linear = F * L / (A * E_eff)
    # Print NR path displacements and small-strain predictions
    print(f"[Neo-Hooke] Displacements u = {u}")
    print(f"[Neo-Hooke] Small-strain tangent E_eff = {E_eff}")
    # Compute strain and small-strain stress approximation
    eps = first.element_strain(L, np.array([u[0], u[1]]))
    sigma_lin = E_eff * eps
    print(f"[Neo-Hooke] Axial strain eps = {eps}")
    print(f"[Neo-Hooke] Approx. axial stress (PK1 tangent) = {sigma_lin}")
    # NR path should converge and be close to PK1 small-strain prediction for small load
    assert "convergence" in sol, "Expected NR convergence metadata in Neo-Hookean run"
    assert np.isclose(u[1], uL_linear, rtol=5e-4, atol=1e-9), (
        f"Neo-Hooke small-load u(L) expected ~{uL_linear}, got {u[1]}"
    )
    # Analytical small-strain check: eps ≈ F/(A*E_eff)
    eps_linear = F / (A * E_eff)
    assert np.isclose(eps, eps_linear, rtol=5e-4, atol=1e-9), (
        f"Neo-Hooke small-load strain expected ~{eps_linear}, got {eps}"
    )