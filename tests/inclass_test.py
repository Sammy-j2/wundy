import io
import wundy
import wundy.first

ss = """
wundy:
  # Request 2 DOFs per node (axial DOF is local index 0)
  dof_per_node: 2
  # Simple 1D bar split into 4 elements (5 nodes)
  nodes:
    - [1, 0.0]
    - [2, 0.25]
    - [3, 0.5]
    - [4, 0.75]
    - [5, 1.0]

  elements:
    - [1, 1, 2]
    - [2, 2, 3]
    - [3, 3, 4]
    - [4, 4, 5]

  # Use a Neo-Hookean material to exercise nonlinear assembly
  materials: {materials}
   

  # Request a nonlinear Newton-Raphson solution and use Gauss quadrature
  integration:
    nonlinear: nonlinear
    internal: gauss
    ngp: 2

  element blocks:
    - name: block
      material: rubber
      elements: [1, 2, 3, 4]
      element:
        type: T1D1
        properties:
          area: 1.0

  boundary conditions:
    - nodes: [1]
      dof: X
      type: DIRICHLET
      value: 0.0
    - nodes: [5]
      dof: X
      type: NEUMANN
      value: 1000.0
"""

def test_NL_elastic_inclass():
    materials = """
    - type: elastic
      name: rubber
      parameters:
        E: 1.0e3
        nu: .3
      density: 1200.0
    """
    f = io.StringIO(ss.format(materials=materials))
    data = wundy.ui.load(f)
    inp = wundy.ui.preprocess(data)
    print(inp)
    sol = wundy.first.first_fe_code(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
    )

    # Check that solution contains displacements
    assert "dofs" in sol
    
    disp = sol["dofs"]
    


    # Check that the displacement at the last node in the axial direction is positive
    assert abs( disp[-1]  - 1.0) < 1e-12, f"Expected positive axial displacement at node 5, got {disp[-1]}"
    
# def test_NL_neohook_inclass():
#     materials = f"""
#     - type: neo_hooke
#       name: rubber
#       parameters:
#         mu: {100/(2*(1+0.3))}
#         lambda: {100*0.3/(1+0.3)/(1-2*0.3)}
#       density: 1200.0
#     """
   
#     f = io.StringIO(ss.format(materials=materials))
#     data = wundy.ui.load(f)
#     inp = wundy.ui.preprocess(data)
#     print(inp)
#     sol = wundy.first.first_fe_code(
#       inp["coords"],
#       inp["blocks"],
#       inp["bcs"],
#       inp["dload"],
#       inp["materials"],
#       inp["block_elem_map"],
#       integration=inp.get("integration"),
#       dof_per_node=inp.get("dof_per_node", 1),
#     )

#     # Check that solution contains displacements
#     assert "dofs" in sol
    
#     disp = sol["dofs"]
    


#     # Check that the displacement at the last node in the axial direction is positive
#     assert abs( disp[-1]  - 1.0) < 1e-12, f"Expected value to be 1, got {disp[-1]}"