import io
import pytest

from wundy import ui


YAML_VALID = """wundy:
  dof_per_node: 2
  nodes:
    - [1, 0.0]
    - [2, 1.0]

  elements:
    - [1, 1, 2]

  materials:
    - type: elastic
      name: steel
      parameters:
        E: 2.1e11

  element blocks:
    - name: block
      material: steel
      elements: [1]
      element:
        type: T1D1
        properties:
          area: 1.0

  boundary conditions:
    - nodes: [1]
      dof: X
      type: DIRICHLET
      value: 0.0
    - nodes: [2]
      dof: Y
      type: DIRICHLET
      value: 0.0
"""

YAML_INVALID = """wundy:
  dof_per_node: 1
  nodes:
    - [1, 0.0]
    - [2, 1.0]

  elements:
    - [1, 1, 2]

  materials:
    - type: elastic
      name: steel
      parameters:
        E: 2.1e11

  element blocks:
    - name: block
      material: steel
      elements: [1]
      element:
        type: T1D1
        properties:
          area: 1.0

  boundary conditions:
    - nodes: [1]
      dof: Y
      type: DIRICHLET
      value: 0.0
"""


def test_preprocess_accepts_dof_per_node():
    data = ui.load(io.StringIO(YAML_VALID))
    pre = ui.preprocess(data)
    assert pre.get("dof_per_node") == 2
    # boundary conditions local_dof should be integers within range
    bcs = pre.get("bcs", [])
    assert len(bcs) == 2
    assert bcs[0]["local_dof"] == 0  # X -> 0
    assert bcs[1]["local_dof"] == 1  # Y -> 1


def test_preprocess_rejects_invalid_dof():
    data = ui.load(io.StringIO(YAML_INVALID))
    with pytest.raises(ValueError):
        ui.preprocess(data)
