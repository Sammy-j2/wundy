import io

import numpy as np

import wundy
from wundy.schemas import DIRICHLET
from wundy.schemas import NEUMANN


def user_input() -> io.StringIO:
    file = io.StringIO()
    file.write("""\
wundy:
  nodes: [[1, 0], [2, 1], [3, 2]]
  elements: [[1, 1, 2], [2, 2, 3]]
  node sets:
  - name: nset-1
    nodes: [2]
  boundary conditions:
  - nodes: 1
    dof: 'x'
  - nodes: 'nset-1'
    dof: 'x'
    value: 1.0
  concentrated loads:
  - nodes: 3
    dof: x
    value: 2.0
  materials:
  - type: neo_hooke
    name: mat-1
    parameters:
      E: 10.0
      nu: 0.3
  distributed loads:
  - type: BX
    elements: all
    value: 2.0
    direction: [1.0]
  distributed loads:
  - type: GRAV
    elements: all
    direction: [-1.0]
    value: 2.0
  element blocks:
  - material: mat-1
    name: block-1
    elements: ALL
    element:
      type: T1D1
""")
    file.seek(0)
    return file


def test_validate_input():
  file = user_input()
  data = wundy.ui.load(file)
  inp = data["wundy"]

  # Input file validator converted lists to arrays
  assert inp["nodes"] == [[1, 0.0], [2, 1.0], [3, 2.0]]
  assert inp["elements"] == [[1, 1, 2], [2, 2, 3]]

  # And inserts default values
  assert inp["boundary conditions"][0] == {
    "nodes": [1],
    "type": DIRICHLET,
    "dof": 0,
    "value": 0.0,
  }
  assert inp["boundary conditions"][1] == {
    "nodes": "NSET-1",
    "type": DIRICHLET,
    "dof": 0,
    "value": 1.0,
  }

  nsets = inp["node sets"]
  assert isinstance(nsets, list)
  assert len(nsets) == 1
  assert np.allclose(nsets[0]["nodes"], [2])

  assert inp["concentrated loads"] == [{"nodes": [3], "dof": 0, "value": 2.0}]

  materials = inp["materials"]
  assert isinstance(materials, list)
  assert len(materials) == 1

  material = materials[0]
  # Accept either classic linear elastic or Neo-Hookean canonical tokens
  assert material["type"] in {"ELASTIC", "NEO_HOOKE", "NEOHOOKE"}
  assert material["name"] == "MAT-1"
  assert material["density"] == 0.0
  assert material["parameters"] == {"E": 10.0, "nu": 0.3}

  blocks = inp["element blocks"]
  assert isinstance(blocks, list)
  assert len(blocks) == 1

  block = blocks[0]
  assert block["element"] == {
    "type": "T1D1",
    "properties": {"area": 1.0},
  }
  assert block["name"] == "BLOCK-1"
  assert block["material"] == "MAT-1"
  assert block["elements"] == "ALL"


def test_preprocess_input():
  file = user_input()
  data = wundy.ui.load(file)
  d = wundy.ui.preprocess(data)

  assert isinstance(d["blocks"], list)
  assert len(d["blocks"]) == 1
  block = d["blocks"][0]
  assert block["name"] == "BLOCK-1"
  assert block["material"] == "MAT-1"
  assert np.allclose(block["connect"], np.array([[0, 1], [1, 2]], dtype=int))
  assert block["element"] == {
    "type": "T1D1",
    "properties": {
      "area": 1.0,
      "node_per_elem": 2,
      "freedom_table": [(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)],
    },
  }

  assert isinstance(d["bcs"], list)

  assert np.allclose(d["coords"], np.array([[0], [1], [2]]))
  assert np.allclose(d["coords"], np.array([[0], [1], [2]]))

  assert isinstance(d["nsets"], dict)
  assert d["nsets"]["NSET-1"] == [1]

  assert isinstance(d["bcs"], list)
  bc = d["bcs"][0]
  assert bc["nodes"] == [0]
  assert bc["local_dof"] == 0
  assert bc["name"] == "BOUNDARY-1"
  assert bc["value"] == 0.0
  assert bc["type"] == DIRICHLET
  bc = d["bcs"][1]
  assert bc["nodes"] == [1]
  assert bc["local_dof"] == 0
  assert bc["name"] == "BOUNDARY-2"
  assert bc["value"] == 1.0
  assert bc["type"] == DIRICHLET
  bc = d["bcs"][2]
  assert bc["nodes"] == [2]
  assert bc["local_dof"] == 0
  assert bc["value"] == 2.0
  assert bc["name"] == "CLOAD-1"
  assert bc["type"] == NEUMANN

  assert isinstance(d["materials"], dict)
  material = d["materials"]["MAT-1"]
  assert material["type"] in {"ELASTIC", "NEO_HOOKE", "NEOHOOKE"}
  assert material["parameters"] == {"E": 10.0, "nu": 0.3}
