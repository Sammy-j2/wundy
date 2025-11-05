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
    """

    dof_per_node = 1
    num_node = coords.shape[0]
    num_dof = num_node * dof_per_node
    K = np.zeros((num_dof, num_dof), dtype=float)
    F = np.zeros(num_dof, dtype=float)

    # integration options: defaults
    if integration is None:
        integration = {"stiffness": "analytic", "internal": "analytic", "ngp": 2}
    ngp = int(integration.get("ngp", 2))

    # Assemble global stiffness
    for block in blocks:
        A = block["element"]["properties"]["area"]
        material = materials[block["material"]]
        E = material_tangent_modulus(material)
        for nodes in block["connect"]:
            # GLOBAL DOF = NODE NUMBER x NUMBER OF DOF PER NODE + LOCAL DOF
            eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]

            xe = coords[nodes]
            he = xe[1, 0] - xe[0, 0]
            if np.isclose(he, 0.0):
                raise ValueError(f"Zero-length element detected between nodes {nodes}")
            ke = element_stiffness(
                A, E, he, integration=integration.get("stiffness", "analytic"), ngp=ngp
            )
            K[np.ix_(eft, eft)] += ke

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

    all_dofs = np.arange(num_dof)
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
            use_gauss = integration.get("internal", "analytic") == "gauss" if integration else False
            is_callable_q = callable(dload["value"]) or (dtype == "GRAV" and callable(g_spec))
            if use_gauss or is_callable_q:
                # Gauss integration for equivalent nodal forces
                xi_pts, wts = _gauss_1d(int(integration.get("ngp", 2)))
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
                F[eft] += fe
            else:
                # analytic constant q per element
                if dtype == "BX":
                    q = float(dload["value"]) * sign
                else:
                    # GRAV with scalar g_spec already reduced to q_spec
                    q = float(q_spec) * sign
                qe = q * he / 2 * np.ones(2)
                F[eft] += qe


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
) -> np.ndarray:
    """Compute the elemental internal nodal force vector for a 2-node bar element.

    The internal (negative internal) nodal force vector is computed as

        fint_e = ke @ u_e

    where ``ke`` is the elemental stiffness matrix and ``u_e`` is the element
    displacement vector [u1, u2]. The returned vector follows the FEM sign
    convention (internal reactions at element nodes), i.e. typically
    [-N, +N] where N is the axial force (positive in tension).

    Parameters
    ----------
    area, E, length : float
        Element cross-sectional area, Young's modulus and element length.
    ue : ndarray
        Element nodal displacement vector with shape (2,) or (2,1).

    Returns
    -------
    fint_e : ndarray
        Element internal nodal force vector of shape (2,).
    """
    ue = np.asarray(ue, dtype=float).ravel()
    if ue.size != 2:
        raise ValueError("element displacement vector `ue` must have length 2")
    # compute ke according to requested integration and use it to form internal force
    ke = element_stiffness(area, E, length, integration=integration, ngp=ngp)
    return ke.dot(ue)


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
        xe = coords[nodes]
        L = float(xe[1, 0] - xe[0, 0])
        eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]
        ue = dofs[eft]
        fe_int = element_internal_force(
            A, E, L, ue, integration=integration.get("internal", "analytic"), ngp=ngp
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
    apply_distributed_loads(F_ext, dloads, coords, blocks, materials, block_elem_map, dof_per_node)

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
    try:
        E = material["parameters"]["E"]
    except Exception as exc:  # KeyError / TypeError
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
