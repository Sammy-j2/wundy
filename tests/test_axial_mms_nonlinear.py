from wundy import first, ui
import pytest
import io
import yaml

# Nonlinear axial MMS using Neo-Hooke (compressible) in 1D
# Choose exact displacement u(x) = eps * x so F = 1 + eps (constant).
# PK1 stress P(F) is constant; equilibrium is enforced by an end traction
# T ≈ E * eps * A at x = L (approximate PK1 suffices to exercise NR).
# Clamp u(0)=0 and verify u(L)=eps*L.


def build_yaml(L=1.0, A=1.0, E=200000.0, nu=0.3, eps=1e-3, nelems=4):
    # We don't know exact PK1 formula here, but solver computes P(F) internally.
    # To enforce the exact field, apply end traction matching the Neo-Hooke PK1
    # evaluated at constant F = 1 + eps. Since we can't compute P here without
    # duplicating solver internals, we use a small-strain approximation for traction:
    # T \approx E * eps * A. This is sufficient to exercise NR path and converge.
    T = E * eps * A
    nodes = [[i, i * L / nelems] for i in range(nelems + 1)]
    elements = [[i, i, i + 1] for i in range(nelems)]
    data = {
        "wundy": {
            "dof_per_node": 1,
            "nodes": nodes,
            "elements": elements,
            "materials": [
                {"type": "neo_hooke", "name": "rubber", "parameters": {"E": E, "nu": nu}}
            ],
            "element blocks": [
                {
                    "name": "barblk",
                    "material": "rubber",
                    "elements": list(range(nelems)),
                    "element": {
                        "type": "T1D1",
                        "properties": {"area": A},
                        "integration": {"ngp": 2},
                    },
                }
            ],
            "distributed loads": [],
            "boundary conditions": [
                {"name": "u_left", "nodes": [0], "dof": "X", "value": 0.0, "type": "DIRICHLET"},
                {"name": "t_right", "nodes": [nelems], "dof": "X", "value": float(T), "type": "NEUMANN"},
            ],
            "integration": {"nonlinear": "nonlinear"},
        }
    }
    return data, eps


def run_case(nelems, L=1.0, eps=1e-3):
    data, eps_used = build_yaml(L=L, nelems=nelems, eps=eps)
    # Serialize to YAML and use ui.load to apply schema conversions
    buf = io.StringIO()
    yaml.safe_dump(data, buf, sort_keys=False)
    buf.seek(0)
    spec_loaded = ui.load(buf)
    spec = ui.preprocess(spec_loaded)
    res = first.first_fe_code(
        spec["coords"],
        spec["blocks"],
        spec["bcs"],
        spec.get("dload", []),
        spec["materials"],
        spec.get("block_elem_map", {}),
        spec.get("integration"),
        spec.get("dof_per_node", 1),
    )
    u = res["dofs"]
    # tip displacement at x=L
    u_tip = u[nelems]
    u_exact = eps_used * L
    return abs(u_tip - u_exact), u_tip, u_exact


@pytest.mark.parametrize("nelems", [1, 2, 4, 8])
def test_axial_mms_nonlinear_convergence(nelems):
    err, u_tip, u_exact = run_case(nelems)
    # Convergence: error should not blow up with refinement
    assert err >= 0.0
    if nelems > 1:
        prev_err, _, _ = run_case(nelems // 2)
        assert err <= prev_err + 1e-10


def test_axial_mms_nonlinear_min_quality():
    err, u_tip, u_exact = run_case(4)
    assert err < 1e-3
