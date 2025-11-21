"""1D finite-element helpers and a simple FE solver used by `wundy`.

This module provides a compact 1D bar/truss solver `first_fe_code` and a
collection of element helpers used for assembly, load application, and
material handling. The implementation intentionally targets simple axial
elements and supports a uniform `dof_per_node` option which expands the
standard 2×2 axial element into a `2 * dof_per_node` element matrix/vector
by placing axial contributions into the local axial DOF (index 0) of each
node. Non-axial DOFs are left uncoupled by design (see docs/material_models.md
for details and rationale).

Nonlinear solver notes:
- The solver includes a Newton–Raphson path for compressible Neo‑Hookean
    materials. Newton controls (`nr_tol`, `nr_max_it`, and the optional
    `nr_du_tol`) are read from the `integration` mapping (global or
    per-block) and can be used to tune convergence behavior.

Key conventions:
- Nodes are indexed 0..N-1.
- `dof_per_node` is a uniform integer across the mesh.
- Local axial DOF at a node is assumed to have index 0 (e.g., `X` when
    `dof_per_node == 2` with local DOFs `[X, Y]`).
- To avoid singular global systems when non-axial DOFs are uncoupled, the
    solver automatically prescribes any DOF whose global stiffness row is
    entirely zero to value `0.0` (automatic zero-prescription).
"""

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .schemas import DIRICHLET
from .schemas import NEUMANN


def first_fe_code(
    coords: NDArray[float],
    blocks: list[dict],
    bcs: list[dict],
    dloads: list[dict],
    materials: dict[str, Any],
    block_elem_map: dict[int, tuple[int, int]],
    integration: dict | None = None,
    dof_per_node: int | None = None,
) -> dict[str, Any]:
    r"""
    Assemble and solve a 1D linear finite element (FE) problem for axial deformation.

    This function performs the full finite element procedure for a 1D bar under
    axial loading. It assembles the global stiffness matrix `K` and force vector `F`
    from element, boundary condition, and material definitions, then applies both
    Neumann (traction/force) and Dirichlet (displacement) boundary conditions.
    The system of equations

        K * u = F

    is solved for the nodal displacements `u` using static condensation of prescribed
    degrees of freedom. The final nodal displacements, stiffness matrix, and force
    vector are returned.

    Parameters
    ----------
    coords : (num_nodes, 1) NDArray[float]
        Array of nodal coordinates along the 1D domain.
        Each row corresponds to the x-position of a node.
    blocks : list of dict
        Element block definitions. Each block contains:
            - "element": dict with "properties" → {"area": float}
            - "material": str name corresponding to `materials` entry
            - "connect": list of tuples, each containing the node indices
              (start, end) defining each element.
    bcs : list of dict
        Boundary condition definitions. Each dictionary includes:
            - "type": either `DIRICHLET` or `NEUMANN`
            - "nodes": list of node indices
            - "local_dof": integer (typically 0 for 1D problems)
            - "value": prescribed displacement (Dirichlet) or force (Neumann)
    dloads : list of dict
        Distributed load definitions. Each dictionary includes:
            - "type": "BX" for uniform body load, or "GRAV" for gravity-based load
            - "direction": array-like of length 1 indicating load direction (+1 or −1)
            - "value": load magnitude
            - "elements": list of element IDs the load applies to
            - "name": descriptive string used for error reporting
    materials : dict[str, Any]
        Material database mapping material names to property dictionaries.
        Each material must define:
            - "parameters": {"E": float}  (Young’s modulus)
            - "density": float (only needed for gravity-type distributed loads)
    block_elem_map : dict[int, tuple[int, int]]
        Maps global element IDs (used in `dloads`) to the tuple
        (block_index, local_element_index_within_block).

    Returns
    -------
    solution : dict[str, Any]
        Dictionary containing:
            - "dofs" : NDArray[float]
                Solved global nodal displacement vector.
            - "stiff" : NDArray[float]
                Assembled global stiffness matrix.
            - "force" : NDArray[float]
                Assembled global force vector (after loads and BCs applied).

    Raises
    ------
    ValueError
        If element length is zero, a distributed load references an unknown
        element, or load direction is invalid.
    NotImplementedError
        If a distributed load type is not supported (only "BX" and "GRAV" are allowed).

    Notes
    -----
    Each 1D bar element contributes the local stiffness matrix

        ke = (A * E / L) * [[ 1, -1 ],
                            [ -1,  1 ]]

    and, for distributed loads `q`, the equivalent nodal force vector

        fe = q * L / 2 * [1, 1]^T.

    Dirichlet boundary conditions are applied via partitioning:
    free and prescribed DOFs are separated, and the reduced system is solved as

        K_ff * u_f = F_f - K_fp * u_p

    before back-substituting to obtain the complete displacement field.

    Examples
    --------
    >>> coords = np.array([[0.0], [1.0]])
    >>> blocks = [{
    ...     "element": {"properties": {"area": 1.0}},
    ...     "material": "steel",
    ...     "connect": [(0, 1)]
    ... }]
    >>> materials = {"steel": {"parameters": {"E": 210e9}, "density": 7850.0}}
    >>> bcs = [
    ...     {"type": DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},
    ...     {"type": NEUMANN, "nodes": [1], "local_dof": 0, "value": 1000.0},
    ... ]
    >>> dloads = []
    >>> block_elem_map = {0: (0, 0)}
        >>> result = first_fe_code(coords, blocks, bcs, dloads, materials, block_elem_map)
        >>> result["dofs"]
        array([0.0, 4.7619e-06])  # approximate

        Additional notes about nonlinear materials
        -----------------------------------------

        If any material in the `materials` database has a type containing both
        the tokens `neo` and `hook` (case-insensitive, underscores or hyphens
        allowed), the solver will switch to a nonlinear Newton–Raphson solution
        path for the whole problem. In that case:

        - Element internal forces are formed using a 1‑D compressible
            Neo‑Hookean model and a consistent tangent (see `_neo_PK1_and_tangent`).
        - The Newton iteration uses a default residual tolerance of 1e-8 and a
            maximum of 25 iterations. If Newton fails to converge the solver
            raises a RuntimeError describing the final residual.

        At present the nonlinear path is automatic and uses the material
        parameters supplied in the YAML (either E/nu or mu/lambda/mu/kappa).
        Future releases may expose NR controls (tolerance, max iterations)
        as explicit API arguments.

        """
    # Degrees of freedom per node: accept caller-provided value or default to 1
    if dof_per_node is None:
        dof_per_node = 1
    num_node = coords.shape[0]
    num_dof = num_node * dof_per_node
    K = np.zeros((num_dof, num_dof), dtype=float)
    F = np.zeros(num_dof, dtype=float)

    # integration options: defaults
    if integration is None:
        integration = {
            "stiffness": "analytic",
            "internal": "analytic",
            "ngp": 2,
            # how to treat nonlinear materials: 'linearize'|'nonlinear'|'auto'
            "nonlinear": "linearize",
        }
    ngp = int(integration.get("ngp", 2))

    # Decide whether a nonlinear (Neo-Hookean) solution is required.
    def is_neo_material(m: dict[str, Any]) -> bool:
        t = str(m.get("type", "")).upper()
        tn = t.replace("-", "").replace("_", "")
        return "NEO" in tn and "HOOK" in tn
    has_neo = any(is_neo_material(mat) for mat in materials.values())

    # Decide actual runtime mode based on integration['nonlinear']:
    # - 'linearize' (default): treat Neo-Hookean materials as linear by using
    #    their small-strain tangent (no NR loop).
    # - 'nonlinear': run full Newton-Raphson when Neo-Hooke materials exist.
    # - 'auto': preserve legacy behaviour (run NR if any Neo material present).
    nonlinear_mode = str(integration.get("nonlinear", "linearize")).lower()
    if nonlinear_mode == "auto":
        do_nonlin = has_neo
    elif nonlinear_mode == "nonlinear":
        do_nonlin = has_neo
    else:
        # 'linearize' or any other value => linear solution using tangent modulus
        do_nonlin = False

    # If requested, solve using Newton-Raphson with consistent tangent stiffness per element
    # (1D compressible Neo-Hookean).
    if do_nonlin:
        # Build constant external force vector (Neumann + dloads) and prescribed lists
        F_ext, prescribed_dofs, prescribed_vals = assemble_external_forces(
            coords, blocks, bcs, dloads, materials, block_elem_map, dof_per_node, integration=integration
        )

        # Initialize displacement vector (prescribed DOFs set to prescribed_vals)
        u = np.zeros(num_dof, dtype=float)
        for pd, pv in zip(prescribed_dofs, prescribed_vals):
            u[pd] = pv

        # Newton-Raphson iterations
        # Newton-Raphson controls (can be overridden via integration dict)
        tol = float(integration.get("nr_tol", 1e-8))
        max_it = int(integration.get("nr_max_it", 25))
        # Optional displacement-increment tolerance: stop when ||du|| < nr_du_tol
        nr_du_tol = integration.get("nr_du_tol", None)
        if nr_du_tol is not None:
            try:
                nr_du_tol = float(nr_du_tol)
            except Exception:
                nr_du_tol = None
        for it in range(max_it):
            # Reset global tangent and internal force
            K[:, :] = 0.0
            F_int = np.zeros_like(F_ext)

            # Assemble internal forces and tangent stiffness
            # Note: when `dof_per_node > 1` each element contribution is
            # expanded from the axial 2-entry/vector and 2x2/matrix into
            # a full `2 * dof_per_node` sized element vector/matrix, placing
            # axial contributions at local axial DOF index 0 of each node.
            for block in blocks:
                # Merge per-block integration options (block may contain its own 'integration')
                block_integration = dict(integration or {})
                block_integration.update(block.get("integration", {}))
                ngp_block = int(block_integration.get("ngp", 2))
                A = block["element"]["properties"]["area"]
                mat = materials[block["material"]]
                for nodes in block["connect"]:
                    eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]
                    xe = coords[list(nodes)]
                    he = xe[1, 0] - xe[0, 0]
                    if np.isclose(he, 0.0):
                        raise ValueError(f"Zero-length element detected between nodes {nodes}")

                    # Element nodal displacements in current configuration
                    ue = u[eft]
                    nd = int(dof_per_node)
                    # axial DOF is assumed to be local index 0 at each node
                    u1 = float(ue[0])
                    u2 = float(ue[nd])
                    # Deformation gradient for 1D element (constant over element)
                    F_e = 1.0 + (u2 - u1) / float(he)

                    # Compute first Piola-Kirchhoff stress P and its derivative dP/dF
                    if is_neo_material(mat):
                        P, dPdF = _neo_PK1_and_tangent(mat, F_e)
                    else:
                        # Material is linear elastic; use P = E*(F-1), dP/dF = E
                        E = material_tangent_modulus(mat)
                        P = E * (F_e - 1.0)
                        dPdF = float(E)

                    # Element internal nodal force (reference configuration)
                    # Element internal nodal force (reference configuration)
                    # Build axial (2-entry) internal force then expand to multi-DOF
                    fe_ax = A * P * np.array([-1.0, 1.0], dtype=float)
                    # Expand to full element DOF vector (2*dof_per_node)
                    nd = dof_per_node
                    fe_int_big = np.zeros(2 * nd, dtype=float)
                    # place axial contributions into local axial DOF (index 0)
                    fe_int_big[0] = fe_ax[0]
                    fe_int_big[nd] = fe_ax[1]
                    F_int[eft] += fe_int_big

                    # Consistent element tangent: 2x2 axial ke, then expand
                    ke_ax = (A * dPdF / float(he)) * np.array([[1.0, -1.0], [-1.0, 1.0]], dtype=float)
                    ke_big = np.zeros((2 * nd, 2 * nd), dtype=float)
                    ke_big[0, 0] = ke_ax[0, 0]
                    ke_big[0, nd] = ke_ax[0, 1]
                    ke_big[nd, 0] = ke_ax[1, 0]
                    ke_big[nd, nd] = ke_ax[1, 1]
                    K[np.ix_(eft, eft)] += ke_big

            # Residual
            R = F_ext - F_int
            # Automatically constrain DOFs with zero stiffness (uncoupled DOFs)
            # Any global DOF that has a completely zero stiffness row is
            # effectively uncoupled by the element formulation (common when
            # using `dof_per_node>1` with axial-only elements). To avoid
            # singular reduced systems we conservatively treat such DOFs as
            # prescribed to zero (same as applying a Dirichlet 0.0).
            all_dofs = np.arange(num_dof)
            zero_rows = np.where(np.all(np.isclose(K, 0.0, atol=1e-12), axis=1))[0]
            for zr in zero_rows:
                if zr not in prescribed_dofs:
                    prescribed_dofs.append(int(zr))
                    prescribed_vals.append(0.0)

            # Restrict to free DOFs and solve for increment
            free_dofs = np.setdiff1d(all_dofs, prescribed_dofs)
            Kff = K[np.ix_(free_dofs, free_dofs)]
            Rf = R[free_dofs]

            # Solve linear system for increment on free DOFs
            try:
                du_f = np.linalg.solve(Kff, Rf)
            except np.linalg.LinAlgError as exc:
                raise np.linalg.LinAlgError(f"Tangent stiffness singular at iteration {it}") from exc

            # Update solution
            u[free_dofs] += du_f

            # Convergence checks
            rnorm = np.linalg.norm(Rf)
            # compute full increment norm (including zeros for prescribed DOFs)
            du_full = np.zeros(num_dof, dtype=float)
            du_full[free_dofs] = du_f
            du_norm = np.linalg.norm(du_full)

            # Check residual tolerance first, then (optionally) displacement increment
            if rnorm < tol:
                # compute full residual (external - internal) for return
                F_int_full = assemble_internal_forces(
                    u, coords, blocks, materials, block_elem_map, dof_per_node, integration=integration
                )
                F_ext_full, _, _ = assemble_external_forces(
                    coords, blocks, bcs, dloads, materials, block_elem_map, dof_per_node, integration=integration
                )
                R_full = F_ext_full - F_int_full
                solution = {"dofs": u, "stiff": K, "force": F_ext, "residual": R_full}
                solution.setdefault("convergence", {})
                solution["convergence"]["reason"] = "residual"
                solution["convergence"]["rnorm"] = float(rnorm)
                return solution
            if nr_du_tol is not None and du_norm < nr_du_tol:
                # Converged by displacement-increment tolerance
                F_int_full = assemble_internal_forces(
                    u, coords, blocks, materials, block_elem_map, dof_per_node, integration=integration
                )
                F_ext_full, _, _ = assemble_external_forces(
                    coords, blocks, bcs, dloads, materials, block_elem_map, dof_per_node, integration=integration
                )
                R_full = F_ext_full - F_int_full
                solution = {"dofs": u, "stiff": K, "force": F_ext, "residual": R_full}
                solution.setdefault("convergence", {})
                solution["convergence"]["reason"] = "du_tol"
                solution["convergence"]["du_norm"] = float(du_norm)
                solution["convergence"]["rnorm"] = float(rnorm)
                return solution

        # If we get here, NR did not converge
        raise RuntimeError(f"Newton-Raphson did not converge after {max_it} iterations; residual={rnorm}")

    # Nonlinear not requested: previous linear assembly path
    # Assemble global stiffness
    for block in blocks:
        # merge per-block integration options
        block_integration = dict(integration or {})
        block_integration.update(block.get("integration", {}))
        ngp_block = int(block_integration.get("ngp", 2))
        A = block["element"]["properties"]["area"]
        material = materials[block["material"]]
        E = material_tangent_modulus(material)
        for nodes in block["connect"]:
            # GLOBAL DOF = NODE NUMBER x NUMBER OF DOF PER NODE + LOCAL DOF
            eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]

            xe = coords[list(nodes)]
            he = xe[1, 0] - xe[0, 0]
            if np.isclose(he, 0.0):
                raise ValueError(f"Zero-length element detected between nodes {nodes}")
            # obtain the 2x2 axial stiffness and expand to multi-DOF element
            # (axial contributions placed at local DOF index 0 of each node)
            ke_ax = element_stiffness(
                A, E, he, integration=block_integration.get("stiffness", "analytic"), ngp=ngp_block
            )
            nd = dof_per_node
            ke_big = np.zeros((2 * nd, 2 * nd), dtype=float)
            ke_big[0, 0] = ke_ax[0, 0]
            ke_big[0, nd] = ke_ax[0, 1]
            ke_big[nd, 0] = ke_ax[1, 0]
            ke_big[nd, nd] = ke_ax[1, 1]
            K[np.ix_(eft, eft)] += ke_big

    # Apply boundary conditions (Neumann forces and Dirichlet prescriptions)
    prescribed_dofs, prescribed_vals = apply_boundary_conditions(K, F, bcs, dof_per_node)

    # Apply distributed loads (moved to helper for clarity and reuse)
    apply_distributed_loads(
        F, dloads, coords, blocks, materials, block_elem_map, dof_per_node, integration=integration
    )

    # Apply Dirchlet boundary conditions using a symmetry preserving elimination
    # Let
    #   Ku = f
    # split dofs into two sets:
    #   1. free
    #   2. prescribed
    # Set up new system:
    #
    #  | K_ff  K_fp |  [ u_f ]   | F_f |
    #  | K_pf  K_pp |  [ u_p ]   | F_p |
    #
    # Eliminate prescribed dofs:
    #   K_ff.u_f = Ff - K_fp.u_p
    # `prescribed_dofs` and `prescribed_vals` are returned from
    # `apply_boundary_conditions` and used below for elimination.

    # Automatically constrain DOFs with zero stiffness (uncoupled DOFs)
    all_dofs = np.arange(num_dof)
    zero_rows = np.where(np.all(np.isclose(K, 0.0, atol=1e-12), axis=1))[0]
    for zr in zero_rows:
        if zr not in prescribed_dofs:
            prescribed_dofs.append(int(zr))
            prescribed_vals.append(0.0)

    free_dofs = np.setdiff1d(all_dofs, prescribed_dofs)
    Kff = K[np.ix_(free_dofs, free_dofs)]
    Kfp = K[np.ix_(free_dofs, prescribed_dofs)]
    Ff = F[free_dofs] - np.dot(Kfp, prescribed_vals)
    uf = np.linalg.solve(Kff, Ff)

    # solve the system
    dofs = np.zeros(num_dof, dtype=float)
    dofs[free_dofs] = uf
    dofs[prescribed_dofs] = prescribed_vals

    solution = {"dofs": dofs, "stiff": K, "force": F}

    # compute and return residual = external - internal for inspection
    F_int_full = assemble_internal_forces(
        dofs, coords, blocks, materials, block_elem_map, dof_per_node, integration=integration
    )
    F_ext_full, _, _ = assemble_external_forces(
        coords, blocks, bcs, dloads, materials, block_elem_map, dof_per_node, integration=integration
    )
    R_full = F_ext_full - F_int_full
    solution["residual"] = R_full
    return solution


def global_dof(node: int, local_dof: int, dof_per_node: int) -> int:
    r"""Return the global degree of freedom index for a given node and local dof

    NOTE: Assumes elements have uniform degrees of freedom across the mesh.

    Compute the global degree of freedom (DOF) index for a given node and local DOF.

    This helper function maps a local degree of freedom at a specific node
    to its corresponding global index in the assembled finite element system.
    The mapping assumes that each node has the same number of DOFs
    (`dof_per_node`) throughout the mesh. It is typically used during
    assembly of the global stiffness matrix and force vector.

    Parameters
    ----------
    node : int
        The node index (0-based) for which the global DOF index is desired.

    local_dof : int
        The local degree of freedom number at the given node.
        For 1D elements, this is almost always 0.

    dof_per_node : int
        Number of degrees of freedom associated with each node.
        For a simple 1D axial bar or truss problem, this value is 1.

    Returns
    -------
    global_index : int
        The global DOF index corresponding to the specified node and local DOF.

    Notes
    -----
    The relationship between local and global DOFs is defined as

    .. math::
                ext{global\_dof} = \text{node} \times \text{dof\_per\_node} + \text{local\_dof}

    This formula ensures unique indexing for all nodal DOFs in the mesh.

    Examples
    --------
    For a system with one DOF per node (e.g., 1D axial deformation):

    >>> global_dof(0, 0, 1)
    0
    >>> global_dof(1, 0, 1)
    1

    For a 2D problem with 2 DOFs per node (u_x, u_y):

    >>> global_dof(0, 0, 2)
    0  # x-displacement at node 0
    >>> global_dof(0, 1, 2)
    1  # y-displacement at node 0
    >>> global_dof(1, 0, 2)
    2  # x-displacement at node 1
    >>> global_dof(1, 1, 2)
    3  # y-displacement at node 1

    See Also
    --------
    first_fe_code : Uses this function to assemble global stiffness and force matrices.

    """
    return node * dof_per_node + local_dof


def apply_boundary_conditions(
    K: np.ndarray, F: np.ndarray, bcs: list[dict], dof_per_node: int
) -> tuple[list[int], list[float]]:
    """Apply boundary conditions to the system.

    This helper currently applies Neumann BCs (mutating `F` in-place) and
    collects Dirichlet prescribed DOFs and values. `K` is accepted so future
    boundary-condition types can modify the stiffness matrix in-place if needed.
    """
    prescribed_dofs: list[int] = []
    prescribed_vals: list[float] = []

    for bc in bcs:
        btype = bc.get("type")
        if btype == NEUMANN:
            for n in bc["nodes"]:
                I = global_dof(n, bc["local_dof"], dof_per_node)
                F[I] += bc["value"]
        elif btype == DIRICHLET:
            for n in bc["nodes"]:
                I = global_dof(n, bc["local_dof"], dof_per_node)
                prescribed_dofs.append(I)
                prescribed_vals.append(bc["value"])
        else:
            # Unknown BC types are currently ignored. We can raise here if
            # you prefer strict behavior.
            continue

    return prescribed_dofs, prescribed_vals


def apply_distributed_loads(
    F: np.ndarray,
    dloads: list[dict],
    coords: np.ndarray,
    blocks: list[dict],
    materials: dict[str, Any],
    block_elem_map: dict[int, tuple[int, int]],
    dof_per_node: int,
    integration: dict | None = None,
) -> None:
    """Assemble and add equivalent nodal forces for distributed loads.

    This function mutates `F` in-place. It supports the existing 1D dload
    types: "BX" (uniform body load per length) and "GRAV" (density*area*accel).

    Parameters mirror the inputs previously used inline in `first_fe_code`.
    """
    if integration is None:
        integration = {"internal": "analytic", "ngp": 2}
    ngp = int(integration.get("ngp", 2))
    for dload in dloads:
        dtype = dload["type"]
        direction = np.array(dload["direction"], dtype=float)
        if direction.size != 1:
            raise ValueError(f"1D problem expects one direction component, got {direction}")
        sign = np.sign(direction[0])
        if sign == 0.0:
            raise ValueError(f"dload direction must be ±1, got {direction[0]}")
        for eid in dload["elements"]:
            if eid not in block_elem_map:
                raise ValueError(
                    f"Element {eid} in distributed load "
                    f"{dload['name']} not found in any element block"
                )
            block_index, local_index = block_elem_map[eid]
            block = blocks[block_index]
            # merge block-specific integration settings
            block_integration = dict(integration or {})
            block_integration.update(block.get("integration", {}))
            ngp_elem = int(block_integration.get("ngp", 2))
            nodes = block["connect"][local_index]
            xe = coords[nodes]
            he = xe[1, 0] - xe[0, 0]
            A = block["element"]["properties"]["area"]
            # prepare q(x) evaluation: handle gravity and BX
            if dtype == "BX":
                q_spec = dload["value"]
            elif dtype == "GRAV":
                mat = materials[block["material"]]
                rho = mat["density"]
                # gravity uses acceleration value in dload["value"]
                g_spec = dload["value"]
                # represent q(x) = rho * A * g(x)
                if callable(g_spec):

                    def q_func(x, rho=rho, A=A, g_spec=g_spec):
                        return rho * A * g_spec(x)
                else:
                    q_const = float(rho * A * g_spec)
                    q_func = None
                    q_spec = q_const
            else:
                raise NotImplementedError(f"dload type {dtype!r} not supported for 1D")

            eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]

            # If user requested gauss or q is callable, perform numerical integration
            use_gauss = block_integration.get("internal", "analytic") == "gauss"
            is_callable_q = callable(dload["value"]) or (dtype == "GRAV" and callable(g_spec))
            if use_gauss or is_callable_q:
                # Gauss integration for equivalent nodal forces
                xi_pts, wts = _gauss_1d(ngp_elem)
                fe = np.zeros(2, dtype=float)
                x1 = xe[0, 0]
                x2 = xe[1, 0]
                J = he / 2.0
                for xi, w in zip(xi_pts, wts):
                    N = _shape_funcs(xi)
                    x = N[0] * x1 + N[1] * x2
                    if dtype == "BX":
                        qx = (
                            dload["value"](x)
                            if callable(dload["value"])
                            else float(dload["value"]) * sign
                        )
                    else:  # GRAV
                        if callable(g_spec):
                            qx = rho * A * g_spec(x) * sign
                        else:
                            qx = float(q_spec) * sign
                    fe += w * N * qx * J
                # Expand 2-entry elemental force into multi-DOF element vector
                # Place the equivalent nodal nodal force into the local axial
                # DOF (index 0) for each node; non-axial DOFs remain zero.
                nd = int(dof_per_node)
                fe_big = np.zeros(2 * nd, dtype=float)
                fe_big[0] = fe[0]
                fe_big[nd] = fe[1]
                F[eft] += fe_big
            else:
                # analytic constant q per element
                if dtype == "BX":
                    q = float(dload["value"]) * sign
                else:
                    # GRAV with scalar g_spec already reduced to q_spec
                    q = float(q_spec) * sign
                qe = q * he / 2 * np.ones(2)
                nd = int(dof_per_node)
                # Map analytic equivalent nodal loads into axial DOFs only
                qe_big = np.zeros(2 * nd, dtype=float)
                qe_big[0] = qe[0]
                qe_big[nd] = qe[1]
                F[eft] += qe_big


def element_stiffness(
    area: float, E: float, length: float, integration: str = "analytic", ngp: int = 2
) -> np.ndarray:
    """Return the 2x2 elemental stiffness matrix for a 1D linear bar element.

    Supports two integration modes:
      - "analytic" (default): uses the closed-form ke = (A*E/L) * [[1,-1],[-1,1]].
      - "gauss": for constant scalar A and E this returns the same analytic
        result (Gauss integration over the element reduces to the same value).

    Note: Gauss integration with spatially varying or callable `E`/`area` is
    not implemented here because the function does not receive nodal
    coordinates. If you need that functionality, extend this signature to
    accept element nodal coordinates and evaluate material/area callables at
    Gauss points.
    """
    if np.isclose(length, 0.0):
        raise ValueError("Zero-length element provided to element_stiffness")

    integration = (integration or "analytic").lower()
    # For the current API, only scalar/constant area and E are supported.
    if integration == "analytic":
        k = float(area) * float(E) / float(length)
        return k * np.array([[1.0, -1.0], [-1.0, 1.0]], dtype=float)
    elif integration == "gauss":
        # If A and E are scalars, Gauss integration yields the same stiffness.
        if not (isinstance(area, (int, float)) and isinstance(E, (int, float))):
            raise NotImplementedError(
                "Gauss integration with non-scalar area or E is not supported by this helper;"
                " extend element_stiffness to accept nodal coordinates to evaluate callables."
            )
        k = float(area) * float(E) / float(length)
        return k * np.array([[1.0, -1.0], [-1.0, 1.0]], dtype=float)
    else:
        raise ValueError(f"Unknown integration option {integration!r} for element_stiffness")


def _shape_funcs(xi: float) -> np.ndarray:
    """Linear shape functions N1,N2 on reference coordinate xi in [-1,1]."""
    N1 = 0.5 * (1.0 - xi)
    N2 = 0.5 * (1.0 + xi)
    return np.array([N1, N2], dtype=float)


def _gauss_1d(ngp: int):
    """Return Gauss-Legendre points and weights on [-1,1] for ngp=1..3."""
    if ngp == 1:
        return np.array([0.0], dtype=float), np.array([2.0], dtype=float)
    if ngp == 2:
        a = 1.0 / np.sqrt(3.0)
        return np.array([-a, a], dtype=float), np.array([1.0, 1.0], dtype=float)
    if ngp == 3:
        a = np.sqrt(3.0 / 5.0)
        return np.array([-a, 0.0, a], dtype=float), np.array(
            [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0], dtype=float
        )
    raise ValueError("unsupported ngp")


def element_internal_force(
    area: float,
    E: float,
    length: float,
    ue: np.ndarray,
    integration: str = "analytic",
    ngp: int = 2,
    dof_per_node: int = 1,
) -> np.ndarray:
    """Compute the elemental internal nodal force vector for a 2-node bar element.

    This helper supports uniform multi-DOF per node. If `ue` has length 2 it
    behaves as before. If `ue` has length 2*dof_per_node the function extracts
    axial DOFs (local index 0) from each node, computes the 2x2 axial internal
    force, and expands it into a 2*dof_per_node vector placing axial entries
    into the local axial DOF index (0) for each node. Non-axial DOFs receive
    zero internal contribution (uncoupled axial-only behaviour).
    """
    ue = np.asarray(ue, dtype=float).ravel()
    if ue.size == 2:
        ke = element_stiffness(area, E, length, integration=integration, ngp=ngp)
        return ke.dot(ue)

    expected = 2 * int(dof_per_node)
    if ue.size != expected:
        raise ValueError(f"element displacement vector `ue` must have length 2 or {expected}")

    # extract axial components (local axial DOF index 0)
    nd = int(dof_per_node)
    ue_ax = np.array([ue[0], ue[nd]], dtype=float)
    # axial ke and axial internal force
    ke_ax = element_stiffness(area, E, length, integration=integration, ngp=ngp)
    fe_ax = ke_ax.dot(ue_ax)
    # expand into full DOF vector: place fe_ax[0] at local index 0 of node0,
    # and fe_ax[1] at local index 0 of node1
    fe_big = np.zeros(expected, dtype=float)
    fe_big[0] = fe_ax[0]
    fe_big[nd] = fe_ax[1]
    return fe_big


def element_axial_force(area: float, E: float, length: float, ue: np.ndarray) -> float:
    """Return the scalar axial internal force for a 2-node bar element.

    This is the axial result (N) computed as k*(u2 - u1) where k = A*E/L.
    Positive results indicate tension.
    """
    ue = np.asarray(ue, dtype=float).ravel()
    if ue.size != 2:
        raise ValueError("element displacement vector `ue` must have length 2")
    k = area * E / length
    return float(k * (ue[1] - ue[0]))


def assemble_internal_forces(
    dofs: np.ndarray,
    coords: np.ndarray,
    blocks: list[dict],
    materials: dict[str, Any],
    block_elem_map: dict[int, tuple[int, int]],
    dof_per_node: int = 1,
    integration: dict | None = None,
) -> np.ndarray:
    """Assemble the global internal force vector from element internal forces.

    Parameters
    ----------
    dofs : ndarray
        Global displacement vector.
    coords : NDArray[float]
        Nodal coordinates array (num_nodes, ndim).
    blocks : list[dict]
        Element block definitions (same format used elsewhere).
    materials : dict
        Material database mapping names to property dicts.
    block_elem_map : dict
        Mapping from global element id to (block_index, local_element_index).
    dof_per_node : int
        Degrees of freedom per node (default 1).

    Returns
    -------
    F_int : ndarray
        Global internal force vector (same length as global DOF vector).
    """
    num_node = coords.shape[0]
    num_dof = int(num_node * dof_per_node)
    F_int = np.zeros(num_dof, dtype=float)

    # Iterate over all element entries via block_elem_map to ensure we match
    # the same global element numbering used elsewhere.
    if integration is None:
        integration = {"stiffness": "analytic", "internal": "analytic", "ngp": 2}
    ngp = int(integration.get("ngp", 2))

    for eid, (block_index, local_index) in block_elem_map.items():
        block = blocks[block_index]
        nodes = block["connect"][local_index]
        A = block["element"]["properties"]["area"]
        mat = materials[block["material"]]
        E = material_tangent_modulus(mat)
        # merge per-block integration settings
        block_integration = dict(integration or {})
        block_integration.update(block.get("integration", {}))
        ngp_elem = int(block_integration.get("ngp", 2))
        xe = coords[list(nodes)]
        L = float(xe[1, 0] - xe[0, 0])
        eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]
        ue = dofs[eft]
        # `ue` may be length 2 (standard) or 2*dof_per_node (expanded).
        # `element_internal_force` will extract axial DOFs (local index 0)
        # when dof_per_node > 1 and expand the axial internal force into
        # the larger element vector.
        fe_int = element_internal_force(
            A,
            E,
            L,
            ue,
            integration=block_integration.get("internal", "analytic"),
            ngp=ngp_elem,
            dof_per_node=dof_per_node,
        )
        F_int[eft] += fe_int

    return F_int


def assemble_external_forces(
    coords: np.ndarray,
    blocks: list[dict],
    bcs: list[dict],
    dloads: list[dict],
    materials: dict[str, Any],
    block_elem_map: dict[int, tuple[int, int]],
    dof_per_node: int = 1,
    integration: dict | None = None,
) -> tuple[np.ndarray, list[int], list[float]]:
    """Assemble the global external force vector from Neumann BCs and dloads.

    Returns a tuple (F_ext, prescribed_dofs, prescribed_vals) where
    - F_ext is the assembled global external force vector (numpy array),
    - prescribed_dofs is a list of global DOF indices prescribed by Dirichlet BCs,
    - prescribed_vals is the corresponding list of prescribed values.

    This mirrors the loads applied inside `first_fe_code`, but keeps assembly
    separate so callers can inspect or compare internal vs external forces.
    """
    num_node = int(coords.shape[0])
    num_dof = int(num_node * dof_per_node)
    F_ext = np.zeros(num_dof, dtype=float)

    # Apply distributed loads first (they may be element-based)
    apply_distributed_loads(
        F_ext, dloads, coords, blocks, materials, block_elem_map, dof_per_node, integration=integration
    )

    # Apply Neumann BCs and collect Dirichlet lists
    prescribed_dofs: list[int] = []
    prescribed_vals: list[float] = []
    for bc in bcs:
        btype = bc.get("type")
        if btype == NEUMANN:
            for n in bc["nodes"]:
                I = global_dof(n, bc["local_dof"], dof_per_node)
                F_ext[I] += bc["value"]
        elif btype == DIRICHLET:
            for n in bc["nodes"]:
                I = global_dof(n, bc["local_dof"], dof_per_node)
                prescribed_dofs.append(I)
                prescribed_vals.append(bc["value"])
        else:
            continue

    return F_ext, prescribed_dofs, prescribed_vals


def material_tangent_modulus(material: dict) -> float:
    """Return the tangent (Young's) modulus for a material definition.

    Parameters
    ----------
    material : dict
        Material dictionary expected to contain a nested "parameters" mapping
        with key "E" for Young's modulus.

    Returns
    -------
    E : float
        Young's modulus for the material.

    Raises
    ------
    ValueError
        If the material dictionary does not contain the expected key.
    """
    # Allow different material types (ELASTIC, NEO-HOOKE). For Neo-Hookean
    # materials we accept either an (E, nu) pair or Lame parameters (mu, lambda)
    # and compute an equivalent Young's modulus for assembly (linearized).
    params = material.get("parameters", {})
    mtype = str(material.get("type", "ELASTIC")).upper()
    # normalize variants like NEO_HOOKE or NEO-HOOKE
    mtype_norm = mtype.replace("-", "").replace("_", "")
    if mtype_norm == "ELASTIC":
        try:
            E = params["E"]
        except Exception as exc:
            raise ValueError("Material definition must include parameters['E'] for ELASTIC") from exc
        return float(E)
    if "NEO" in mtype_norm and "HOOK" in mtype_norm:
        # If user provided E, use it. Otherwise, compute E from mu and lambda.
        if "E" in params:
            return float(params["E"])
        if "mu" in params and "lambda" in params:
            mu = float(params["mu"])
            lam = float(params["lambda"])
            # Convert Lame parameters to Young's modulus E:
            # E = mu*(3*lambda + 2*mu)/(lambda + mu)
            if (lam + mu) == 0:
                raise ValueError("Invalid Lame parameters: lambda + mu must be non-zero")
            E = mu * (3.0 * lam + 2.0 * mu) / (lam + mu)
            return float(E)
        if "mu" in params and "kappa" in params:
            mu = float(params["mu"])
            kappa = float(params["kappa"])
            lam = kappa - 2.0 / 3.0 * mu
            if (lam + mu) == 0:
                raise ValueError("Invalid parameters: computed lambda + mu is zero")
            E = mu * (3.0 * lam + 2.0 * mu) / (lam + mu)
            return float(E)
        raise ValueError(
            "Neo-Hookean material requires either parameters['E'] or parameters['mu'] and parameters['lambda']"
        )
    # Fallback: try to read E and raise a helpful error if missing
    try:
        E = params["E"]
    except Exception as exc:
        raise ValueError("Material definition must include parameters['E']") from exc
    return float(E)


def material_constitutive(material: dict, ndim: int = 1) -> np.ndarray:
    """Return a small constitutive matrix for the material.

    For the current 1D solver this returns a 1x1 matrix containing Young's
    modulus E. For higher-dimensional use-cases this function can be
    extended to return plane-stress/plane-strain tensors.
    """
    E = material_tangent_modulus(material)
    if ndim == 1:
        return np.array([[E]], dtype=float)
    raise NotImplementedError("material_constitutive currently supports ndim==1 only")


def _neo_PK1_and_tangent(material: dict, F: float) -> tuple[float, float]:
    """Return first Piola-Kirchhoff stress P and tangent dP/dF for 1D

    Uses a compressible Neo-Hookean energy reduced to 1D:

        W = (mu/2)*(F**2 - 3) - mu*ln(F) + (kappa/2)*(ln(F))**2

    so that

        P = dW/dF = mu*(F - 1/F) + (kappa*ln(F))/F

    and

        dP/dF = mu*(1 + 1/F^2) + kappa*(1 - ln(F))/F^2

    The function accepts materials specified either by (E, nu) or (mu, lambda)
    (or mu and kappa). If E/nu are provided they are converted to mu/kappa.
    """
    if F <= 0.0:
        raise ValueError("Invalid deformation gradient F <= 0 for Neo-Hookean material")

    params = material.get("parameters", {})
    # Determine mu and kappa (bulk modulus)
    if "mu" in params and "lambda" in params:
        mu = float(params["mu"])
        lam = float(params["lambda"])
        # convert to bulk modulus
        kappa = lam + 2.0 / 3.0 * mu
    elif "mu" in params and "kappa" in params:
        mu = float(params["mu"])
        kappa = float(params["kappa"])
    elif "E" in params and "nu" in params:
        E = float(params["E"])
        nu = float(params["nu"])
        mu = E / (2.0 * (1.0 + nu))
        kappa = E / (3.0 * (1.0 - 2.0 * nu))
    else:
        raise ValueError(
            "Neo-Hookean material requires either (E,nu) or (mu,lambda) or (mu,kappa) in parameters"
        )

    # compute P and dP/dF
    lnF = float(np.log(F))
    P = mu * (F - 1.0 / F) + (kappa * lnF) / F
    dPdF = mu * (1.0 + 1.0 / (F * F)) + kappa * (1.0 - lnF) / (F * F)
    return float(P), float(dPdF)


def element_strain(length: float, ue: np.ndarray) -> float:
    """Compute the axial strain for a 2-node 1D bar element.

    Strain is assumed constant over the linear element and computed as

        eps = (u2 - u1) / L

    Parameters
    ----------
    length : float
        Element length L (must be non-zero).
    ue : ndarray
        Element nodal displacement vector [u1, u2].

    Returns
    -------
    eps : float
        Axial strain (dimensionless).
    """
    ue = np.asarray(ue, dtype=float).ravel()
    if ue.size != 2:
        raise ValueError("element displacement vector `ue` must have length 2")
    if np.isclose(length, 0.0):
        raise ValueError("Zero-length element provided to element_strain")
    return float((ue[1] - ue[0]) / length)


def element_stress_from_strain(material: dict, strain: float) -> float:
    """Return Cauchy axial stress (= E * strain) for a given material and strain.

    Parameters
    ----------
    material : dict
        Material dictionary with parameters['E'] available.
    strain : float
        Axial strain value.

    Returns
    -------
    sigma : float
        Axial stress (same sign convention as positive tension).
    """
    E = material_tangent_modulus(material)
    return float(E * strain)


def element_stress(area: float, material: dict, length: float, ue: np.ndarray) -> float:
    """Compute the scalar axial stress for a 2-node bar element from nodal displacements.

    This function computes the element strain and multiplies by Young's modulus
    to produce the axial stress. It is equivalent to element_axial_force / area
    (assuming small strains and linear elasticity).
    """
    eps = element_strain(length, ue)
    return element_stress_from_strain(material, eps)
