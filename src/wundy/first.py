"""First-order 1D finite element solver utilities.

This module provides a minimal linear 1D finite element assembly and solver
for two-node (linear) bar elements. It supports assembling the global
stiffness matrix, applying Neumann (point) loads, applying distributed
loads (types "BX" and "GRAV"), and eliminating Dirichlet (prescribed)
degrees of freedom via a symmetry-preserving elimination.

Key functions
- ``first_fe_code``: assemble and solve a 1D linear finite element problem.
- ``global_dof``: compute the global degree-of-freedom index for a node.

Data shapes and conventions
- ``coords``: numpy array of shape (num_nodes, ndim) where ndim == 1 for
    this solver (but stored as Nx1 arrays in the project). Node indices are
    assumed zero-based.
- ``blocks``: list of element blocks; each block contains element
    connectivity under the key ``connect`` and element properties under
    ``element`` -> ``properties``.
- ``bcs``: list of boundary condition dicts with types defined in
    ``wundy.schemas`` (``DIRICHLET`` and ``NEUMANN``).
- ``dloads``: list of distributed-load dicts; supported ``type`` values in
    this module are ``"BX"`` (body/axial load per unit length) and
    ``"GRAV"`` (gravity, requires material density).

Raises
- ``ValueError`` for zero-length elements or invalid distributed-load
    directions.
- ``NotImplementedError`` for unsupported distributed-load types.

Note: this file focuses on clarity and testability rather than performance.
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
) -> dict[str, Any]:
    """Assemble and solve a 1D linear finite element problem.

    Parameters
    ----------
    coords:
        Array of nodal coordinates (shape: (num_nodes, ndim)). Node indices
        are assumed zero-based.
    blocks:
        List of element blocks. Each block must include an ``element``
        mapping with ``properties`` (e.g. ``area``) and a ``connect`` list
        of node index pairs for each element.
    bcs:
        List of boundary condition dicts. Use values from
        ``wundy.schemas`` for ``DIRICHLET`` and ``NEUMANN`` types.
    dloads:
        List of distributed load dicts. Supported types: ``"BX"`` and
        ``"GRAV"`` (requires material ``density``).
    materials:
        Mapping of material name to material definition (parameters,
        density, etc.).
    block_elem_map:
        Mapping from global element id to a tuple ``(block_index,
        local_element_index)`` used to locate element connectivity for
        distributed loads.

    Returns
    -------
    dict
        A solution dictionary with keys:
        - ``dofs``: numpy array of nodal degrees of freedom (displacements).
        - ``stiff``: assembled global stiffness matrix K.
        - ``force``: global force vector F (before Dirichlet elimination).

    Raises
    ------
    ValueError
        On zero-length elements or invalid distributed-load directions.
    NotImplementedError
        If an unsupported distributed-load ``type`` is encountered.

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

    """
    return node * dof_per_node + local_dof
