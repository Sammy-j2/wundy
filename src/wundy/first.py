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
) -> dict[str, Any]:
    """
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

    # Assemble global stiffness
    for block in blocks:
        A = block["element"]["properties"]["area"]
        material = materials[block["material"]]
        E = material["parameters"]["E"]
        for nodes in block["connect"]:
            # GLOBAL DOF = NODE NUMBER x NUMBER OF DOF PER NODE + LOCAL DOF
            eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]

            xe = coords[nodes]
            he = xe[1, 0] - xe[0, 0]
            if np.isclose(he, 0.0):
                raise ValueError(f"Zero-length element detected between nodes {nodes}")
            ke = A * E / he * np.array([[1.0, -1.0], [-1.0, 1.0]])
            K[np.ix_(eft, eft)] += ke

    # Apply Neumann boundary conditions to force
    for bc in bcs:
        if bc["type"] == NEUMANN:
            for n in bc["nodes"]:
                I = global_dof(n, bc["local_dof"], dof_per_node)
                F[I] += bc["value"]

    # Apply distributed loads
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
            if dtype == "BX":
                q = dload["value"] * sign
            elif dtype == "GRAV":
                mat = materials[block["material"]]
                rho = mat["density"]
                q = rho * A * dload["value"] * sign
            else:
                raise NotImplementedError(f"dload type {dtype!r} not supported for 1D")
            eft = [global_dof(n, j, dof_per_node) for n in nodes for j in range(dof_per_node)]
            qe = q * he / 2 * np.ones(2)
            F[eft] += qe

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
    prescribed_dofs: list[int] = []
    prescribed_vals: list[float] = []
    for bc in bcs:
        if bc["type"] == DIRICHLET:
            for n in bc["nodes"]:
                I = global_dof(n, bc["local_dof"], dof_per_node)
                prescribed_dofs.append(I)
                prescribed_vals.append(bc["value"])

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
    """Return the global degree of freedom index for a given node and local dof

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
        \text{global\_dof} = \text{node} \times \text{dof\_per\_node} + \text{local\_dof}

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
