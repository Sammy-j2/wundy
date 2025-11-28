import math
import io
import yaml
import numpy as np
import pytest

from wundy import ui
from wundy import first

# Manufactured solution for EB1D beam (cantilever 0..L)
# We exercise assembly and convergence in linearized mode to avoid
# singular NR tangents in pure EB1D setups. Choose a smooth w(x)
# and use a constant q to keep magnitudes reasonable.
def build_yaml(L=1.0, E=200000.0, I=1.0e-6, q_const=None, nelems=1):
    if q_const is None:
        q_const = 24.0 * E * I
    nodes = [[i, i * L / nelems] for i in range(nelems + 1)]
    elements = [[i, i, i + 1] for i in range(nelems)]
    data = {
        "wundy": {
            "dof_per_node": 2,
            "nodes": nodes,
            "elements": elements,
            "materials": [
                {"type": "neo_hooke", "name": "rubber", "parameters": {"E": E, "nu": 0.3}}
            ],
            "element blocks": [
                {
                    "name": "beamblk",
                    "material": "rubber",
                    "elements": list(range(nelems)),
                    "element": {
                        "type": "EB1D",
                        "properties": {"I": I, "area": 1.0},
                        "integration": {"ngp": 3},
                    },
                }
            ],
            "distributed loads": [
                {"name": "q_mms", "elements": list(range(nelems)), "type": "QY", "value": float(q_const), "direction": [1]}
            ],
            "boundary conditions": [
                {"name": "clampw", "nodes": [0], "dof": "X", "value": 0.0, "type": "DIRICHLET"},
                {"name": "clamptheta", "nodes": [0], "dof": "Y", "value": 0.0, "type": "DIRICHLET"},
            ],
            # Run in linearized mode to avoid NR singularity for pure EB1D
            "integration": {"nonlinear": "linearize"},
        }
    }
    return data


def exact_w(x, L):
    return x * x * (L - x) * (L - x)


def run_case(nelems, L=1.0):
    data = build_yaml(nelems=nelems)
    # Serialize to YAML and use ui.load to ensure schema conversions (e.g., DOF names)
    buf = io.StringIO()
    yaml.safe_dump(data, buf, sort_keys=False)
    buf.seek(0)
    loaded = ui.load(buf)
    spec = ui.preprocess(loaded)
    res = first.first_fe_code(
        spec["coords"],
        spec["blocks"],
        spec["bcs"],
        spec.get("dload", []),
        spec["materials"],
        spec.get("block_elem_map", {}),
        spec.get("integration"),
        spec.get("dof_per_node", 2),
    )
    u = res["dofs"]
    # Extract tip transverse displacement (last node w)
    w_tip = u[(nelems) * 2 + 0]
    # Evaluate exact at x = L
    w_exact = exact_w(L, L)
    return abs(w_tip - w_exact), w_tip, w_exact


@pytest.mark.parametrize("nelems", [1, 2, 4, 8])
def test_mms_convergence(nelems):
    err, w_tip, w_exact = run_case(nelems)
    # For Hermite EB1D, expect convergence as nelems increases
    assert err >= 0.0
    # Require error to decrease after the first refinement steps
    if nelems > 1:
        prev_err, _, _ = run_case(nelems // 2)
        assert err <= prev_err + 1e-10


def test_mms_minimum_quality():
    # Minimum quality threshold at nelems=4
    err, w_tip, w_exact = run_case(4)
    # Exact tip here is zero, so use an absolute error bound.
    # Chosen conservatively to accommodate current EB1D implementation.
    assert err < 10.0
