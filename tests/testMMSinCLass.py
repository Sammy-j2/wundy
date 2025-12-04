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
    expression: -pi**2*sin(pi * x)
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
    
    assert np.allclose( R[0] , -q * L )

  
    
