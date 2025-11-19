import io
import numpy as np

from wundy import ui, first


SIMPLE_NEO_YAML = """wundy:
  integration:
    nonlinear: linearize
    internal: analytic
    ngp: 2

  nodes:
    - [1, 0.0]
    - [2, 1.0]

  elements:
    - [1, 1, 2]

  materials:
    - type: neo_hooke
      name: rubber
      parameters:
        mu: 3.0e5
        lambda: 2.0e5
      density: 1200.0

  element blocks:
    - name: single
      material: rubber
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

  concentrated loads:
    - nodes: [2]
      dof: X
      value: 10000.0
"""


def test_integration_propagation_and_solver_runs():
    # Load YAML and preprocess
    f = io.StringIO(SIMPLE_NEO_YAML)
    data = ui.load(f)
    pre = ui.preprocess(data)

    # Global integration should be present and default to 'linearize'
    assert "integration" in pre
    assert pre["integration"]["nonlinear"] == "linearize"

    # Blocks list should exist and have merged integration settings
    blocks = pre["blocks"]
    assert isinstance(blocks, list) and len(blocks) == 1
    blk = blocks[0]
    assert "integration" in blk
    assert blk["integration"]["nonlinear"] == "linearize"

    # Prepare inputs for first_fe_code
    coords = pre["coords"]
    blocks_input = pre["blocks"]
    bcs = pre["bcs"]
    dloads = pre["dload"]
    materials = pre["materials"]
    block_elem_map = pre["block_elem_map"]

    # Run in linearized mode (default)
    sol_lin = first.first_fe_code(
      coords,
      blocks_input,
      bcs,
      dloads,
      materials,
      block_elem_map,
      integration=pre["integration"],
      dof_per_node=pre.get("dof_per_node", 1),
    )
    assert "dofs" in sol_lin and isinstance(sol_lin["dofs"], np.ndarray)
    u_lin = sol_lin["dofs"]
    assert u_lin.size == coords.shape[0]
    assert np.all(np.isfinite(u_lin))
    # residual should be present and finite; save for inspection
    assert "residual" in sol_lin
    r_lin = sol_lin["residual"]
    assert isinstance(r_lin, np.ndarray)
    assert r_lin.size == u_lin.size
    assert np.all(np.isfinite(r_lin))
    np.save("tests/last_residual_linear.npy", r_lin)

    # Now request full nonlinear NR globally and ensure solver runs
    pre_nl = dict(pre)
    pre_nl["integration"] = dict(pre_nl["integration"])
    pre_nl["integration"]["nonlinear"] = "nonlinear"

    sol_nl = first.first_fe_code(
      coords,
      blocks_input,
      bcs,
      dloads,
      materials,
      block_elem_map,
      integration=pre_nl["integration"],
      dof_per_node=pre.get("dof_per_node", 1),
    )
    assert "dofs" in sol_nl and isinstance(sol_nl["dofs"], np.ndarray)
    u_nl = sol_nl["dofs"]
    assert u_nl.size == coords.shape[0]
    assert np.all(np.isfinite(u_nl))
    assert "residual" in sol_nl
    r_nl = sol_nl["residual"]
    assert isinstance(r_nl, np.ndarray)
    assert r_nl.size == u_nl.size
    assert np.all(np.isfinite(r_nl))
    np.save("tests/last_residual_nonlinear.npy", r_nl)

    # For a soft Neo-Hooke and a reasonably large load we expect some difference
    # between linearized and full nonlinear responses (not necessarily large).
    assert not np.allclose(u_lin, u_nl)
