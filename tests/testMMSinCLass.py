import io

import numpy as np

import wundy
import wundy.first


def test_first_2():
    file = io.StringIO()
    L = 2
    num_elems = 50
    q = -np.pi**2 * np.sin(np.pi)  # so that exact solution is u = sin(pi*x)
    nodes = [[i+1,float(u)] for i,u in enumerate(np.linspace(0, L, num_elems + 1, dtype=float))]
    elements = np.array([[i+1, i+1, i+2] for i in range(num_elems)], dtype=int)
    file.write("""\
wundy:
  nodes: {nodes}
  elements: {elements}
  boundary conditions:
  - name: fix-nodes
    dof: x
    nodes: [1]
  - name: fix-right
    dof: x
    nodes: [51]
  
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

  distributed loads:
  - name: dload-1
    elements: all
    type: bx
    expression: 10*pi**2*sin(pi * x)
    direction: [1]
""".format(nodes=nodes, elements=elements.tolist(), q=q))
    file.seek(0)
    data = wundy.ui.load(file)
    inp = wundy.ui.preprocess(data)
    soln = wundy.first.first_fe_code(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
    )

    dofs = soln["dofs"]
    K = soln["stiff"]
    F = soln["force"]
    R = np.dot(K, dofs) - F
    
    # For sinusoidal load with zero net integral, reaction forces should sum to zero
    total_reaction = R[0] + R[-1]
    assert np.allclose(total_reaction, 0.0, atol=1e-10)
    
    # Verify computed displacements match manufactured solution u = sin(pi*x)
    u = dofs
    ue = np.sin( np.pi * np.linspace(0,L,num_elems+1) )
    norm = np.linalg.norm(ue-u)
    print("Norm of error:", norm)
    print(f"Computed u first 5: {u[:5]}")
    print(f"Expected ue first 5: {ue[:5]}")
    print(f"Computed u last 5: {u[-5:]}")
    print(f"Expected ue last 5: {ue[-5:]}")
    
    # Check displacement error is small
    assert np.allclose(u, ue, atol=1e-3, rtol=1e-2)
    
    
 

  
    