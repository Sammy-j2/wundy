import io

import numpy as np

import wundy.first as first
import wundy.ui as ui


def test_element_stiffness():
    A = 1.0
    E = 10.0
    L = 2.0
    ke = first.element_stiffness(A, E, L)
    k = A * E / L
    assert np.allclose(ke, k * np.array([[1.0, -1.0], [-1.0, 1.0]]))
    # Gauss result should match analytic for constant A,E
    ke_g = first.element_stiffness(A, E, L, integration="gauss", ngp=2)
    assert np.allclose(ke_g, ke)


def test_element_internal_and_axial():
    A = 1.0
    E = 20.0
    L = 0.5
    ue = np.array([0.1, 0.3])
    fe_int = first.element_internal_force(A, E, L, ue)
    N = first.element_axial_force(A, E, L, ue)
    # internal vector should be [-N, +N]
    assert np.allclose(fe_int, np.array([-N, N]))


def test_material_helpers():
    mat = {"parameters": {"E": 15.0}}
    assert first.material_tangent_modulus(mat) == 15.0
    C = first.material_constitutive(mat, ndim=1)
    assert np.allclose(C, np.array([[15.0]]))


def _build_simple_input(with_dload: bool = False):
    file = io.StringIO()
    file.write("""\
wundy:
  nodes: [[1, 0], [2, 1], [3, 2], [4, 3], [5, 4]]
  elements: [[1, 1, 2], [2, 2, 3], [3, 3, 4], [4, 4, 5]]
  boundary conditions:
  - name: fix-nodes
    dof: x
    nodes: [1]
  concentrated loads:
  - name: cload-1
    nodes: [5]
    value: 2.0
  materials:
  - type: elastic
    name: mat-1
    parameters:
      E: 10.0
      nu: 0.3
  element blocks:
  - material: mat-1
    name: block-1
    elements: all
    element:
      type: t1d1
      properties:
        area: 1
""")
    if with_dload:
        # append distributed load block
        file.write("""\
  distributed loads:
  - name: dload-1
    elements: all
    type: BX
    value: 8
    direction: [1]
""")
    file.seek(0)
    data = ui.load(file)
    inp = ui.preprocess(data)
    return inp


def test_internal_external_equilibrium_no_dload():
    inp = _build_simple_input(with_dload=False)
    soln = first.first_fe_code(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
    )
    F_ext, pres_dofs, pres_vals = first.assemble_external_forces(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
    )
    F_int = first.assemble_internal_forces(
        soln["dofs"], inp["coords"], inp["blocks"], inp["materials"], inp["block_elem_map"]
    )
    # For prescribed DOFs the solver does not include reaction forces in the
    # assembled external vector; compare only free DOFs for equilibrium here.
    all_dofs = np.arange(inp["coords"].shape[0])
    free_dofs = np.setdiff1d(all_dofs, pres_dofs)
    assert np.allclose(F_ext[free_dofs], F_int[free_dofs])


def test_internal_external_equilibrium_with_dload():
    inp = _build_simple_input(with_dload=True)
    soln = first.first_fe_code(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
    )
    F_ext, pres_dofs, pres_vals = first.assemble_external_forces(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
    )
    F_int = first.assemble_internal_forces(
        soln["dofs"], inp["coords"], inp["blocks"], inp["materials"], inp["block_elem_map"]
    )
    all_dofs = np.arange(inp["coords"].shape[0])
    free_dofs = np.setdiff1d(all_dofs, pres_dofs)
    assert np.allclose(F_ext[free_dofs], F_int[free_dofs])
