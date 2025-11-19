# Material models (wundy)

This page documents how to specify material models for use with the 1‑D
finite‑element solver (`first_fe_code`) in this repository.

Contents
- What the solver expects
- YAML example
- Programmatic (Python) example
- How materials are used by the solver
- Extending/custom materials

## What the solver expects

The solver accepts a material database as a Python mapping (or the
equivalent YAML section). Each material is a dictionary with at least the
following structure:

- `parameters`: a mapping containing at least `E` (Young's modulus).
- `density`: optional float, used by gravity (`GRAV`) distributed loads.

Required keys (summary):

- `parameters['E']` — numeric Young's modulus (Pa).
- `density` — optional (kg/m^3) for gravity loads.

The relevant helper functions in `src/wundy/first.py` expect this layout:

- `material_tangent_modulus(material: dict) -> float` — returns
  `material['parameters']['E']` (raises ValueError if missing).
- `material_constitutive(material: dict, ndim: int=1) -> np.ndarray` —
  returns a small constitutive matrix (1×1 array with `E` for 1‑D).

## YAML example

When using the YAML input format (see `docs/YAML_User_Input_Spec_README.md`),
materials are given like this:

```yaml
materials:
  - type: elastic
    name: steel_a
    parameters:
      E: 210e9      # Pa
      nu: 0.3       # Poisson's ratio (kept for compatibility)
    density: 7850   # kg/m^3 (optional; used with GRAV-type dloads)
```

Notes:

- `type` is currently informational and normalized to uppercase. The
  implementation presently supports linear elastic (`ELASTIC`) materials.
- `nu` is accepted in the YAML schema for forward compatibility but is not
  used by the 1‑D constitutive helper (it is ignored unless you extend
  `material_constitutive` for higher dimensions).

### Neo‑Hookean example

You can also specify compressible Neo‑Hookean materials. The loader accepts
either the classic `(E, nu)` pair or Lame-style parameters `(mu, lambda)` or
`(mu, kappa)`. By default the solver treats Neo‑Hookean materials in a
linear FEM context by using the small‑strain/tangent modulus (i.e. it
``linearize``s the material). This preserves the usual linear FEM workflow
while still allowing you to specify hyperelastic material parameters.

If you want the full nonlinear Newton–Raphson solution, set the `integration`
option `nonlinear: nonlinear` either globally or per-element-block (see
examples below). The `integration.nonlinear` option accepts three values:

- `linearize` (default): linear FEM using the small-strain tangent modulus.
- `nonlinear`: run Newton–Raphson when any Neo‑Hookean material is present.
- `auto`: legacy behavior — equivalent to `nonlinear` (keeps old behavior).

Example material entry (Neo‑Hookean parameters):

```yaml
materials:
  - type: neo_hooke
    name: rubber
    parameters:
      mu: 3.0e5
      lambda: 2.0e5
    density: 1200.0
```

Global control of nonlinear behavior (request full NR globally):

```yaml
wundy:
  integration:
    nonlinear: nonlinear    # 'linearize'|'nonlinear'|'auto' (default 'linearize')
    internal: analytic      # 'analytic'|'gauss' for internal force evaluation
    ngp: 2                  # Gauss points when using 'gauss'
```

Per-block (per-element) override example — request NR for a specific block:

```yaml
element blocks:
  - name: rubber_block
    material: rubber
    elements: [1]
    element:
      type: T1D1
      properties:
        area: 1.0
      integration:
        nonlinear: nonlinear
        internal: gauss
        ngp: 3
```

When `integration.nonlinear` is left as `linearize` (the default), Neo‑Hookean
materials are accepted but the solver runs the linear FEM using the small‑strain
tangent modulus computed from the provided parameters (E or mu/lambda). Use
the `nonlinear` option only when you need the full geometric/material nonlinearity
solution returned by the Newton–Raphson loop.

## Programmatic (Python) example

You can provide the same data directly to `first_fe_code` as a Python
dictionary. Example:

```python
import numpy as np
from wundy import first

coords = np.array([[0.0], [1.0]])
blocks = [
    {
        "element": {"properties": {"area": 1.0}},
        "material": "steel",
        "connect": [(0, 1)],
    }
]
materials = {
    "steel": {"parameters": {"E": 210e9, "nu": 0.3}, "density": 7850.0}
}
bcs = [
    {"type": first.DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},
    {"type": first.NEUMANN, "nodes": [1], "local_dof": 0, "value": 1000.0},
]
dloads = []
block_elem_map = {0: (0, 0)}

result = first.first_fe_code(coords, blocks, bcs, dloads, materials, block_elem_map)
u = result["dofs"]

# Compute element stresses (example for the single element):
area = blocks[0]["element"]["properties"]["area"]
L = coords[1, 0] - coords[0, 0]
ue = u[[0, 1]]
sigma = first.element_stress(area, materials["steel"], L, ue)
print("element stress (Pa):", sigma)
```

## How materials are used by the solver

- At assembly time the code looks up the block's material by name and
  calls `material_tangent_modulus(material)` to obtain `E` for use in the
  analytic element stiffness formula ke = (A*E/L) * [[1,-1],[-1,1]].
- For grav‑type distributed loads the solver uses `density` from the
  material definition to compute a per‑length body force `q = rho * A * g`.
- The helper `material_constitutive` returns a 1×1 constitutive matrix
  containing `E` (it can be extended to return tensors for plane stress/strain
  behaviour in higher-dimensional solvers).

## Extension points / custom materials

The code currently assumes small‑strain linear elasticity (constant E per
material). If you need nonlinear materials, spatially varying E, or
different constitutive laws, consider one of the following approaches:

1. Extend `material_tangent_modulus(material)` to handle additional
   material `type` values and compute secant/tangent moduli. Update
   `material_constitutive` accordingly.
2. Replace the `E` parameter with a callable (e.g. `E(x)` or
   `E(strain)`) and modify the element routines to accept and evaluate
   that callable at Gauss points (the solver already contains small
   Gauss helpers `_gauss_1d` and `_shape_funcs` used for dloads)
   — this requires code changes in `element_stiffness` and
   `element_internal_force` to evaluate stiffness/internals by integration.
3. Add a new material `type` (for instance, `HYPERELASTIC`) and
   implement the corresponding constitutive/tangent routines. Then
   update the assembly to call the proper tangent stiffness for each
   element when forming the global stiffness (needed for nonlinear
   Newton iterations).

## Quick tips

- Ensure `parameters.E` is present (float) or `first_fe_code` will raise a
  ValueError when trying to read material parameters.
- If you use `GRAV` dloads, include `density` in the material so the
  solver can compute body forces automatically.
- Keep names consistent: the block's `material` field must match a
  `materials` dictionary key (case‑insensitive in YAML; programmatic
  keys are matched as provided).

## Multi-DOF (`dof_per_node`) support and behavior

The solver accepts a top-level `dof_per_node` integer (via YAML or
programmatic input) which sets the number of degrees of freedom associated
with each node. This enables using the same mesh for problems with multiple
DOFs per node (for example, in-plane `X` and `Y` displacements) while
preserving the simple 1‑D axial element formulation used by `first_fe_code`.

Key points:

- **Uniform DOFs per node:** The implementation assumes a uniform number of
  DOFs per node across the mesh. Set `dof_per_node` in your YAML under the
  top-level `wundy:` mapping, for example `dof_per_node: 2`.
- **Axial DOF index:** For each node the *axial* degree of freedom is
  assumed to be the local DOF with index `0`. That is, when `dof_per_node == 2`
  we assume local DOFs `[0,1]` correspond to `[X, Y]` and axial coupling
  (the bar/truss axial stiffness) only connects the local index `0` DOFs.
- **Axial-only coupling (Option A):** The current multi-DOF support expands
  the standard 2×2 axial element stiffness and 2-entry internal force into
  a `2 * dof_per_node` sized element matrix/vector by placing axial entries
  at local index `0` for each node. Non-axial DOFs are intentionally
  uncoupled (receive zero internal contribution). This keeps the element
  formulation simple and backward compatible.
- **Automatic zero-prescription:** Because non-axial DOFs are uncoupled they
  may produce rows/columns of zero in the global stiffness matrix. To avoid
  singular systems the solver will automatically treat any *completely zero*
  global DOF as a prescribed DOF with value `0.0` (equivalent to applying a
  Dirichlet condition of zero). This is a conservative default: if you need
  free non-axial DOFs (e.g. full 2D mechanics) you should use a genuinely
  vectorial element formulation (Option B) which couples DOFs correctly.

Example YAML snippet (enable two DOFs per node):

```yaml
wundy:
  dof_per_node: 2
  nodes:
    - [1, 0.0]
    - [2, 1.0]
  # ... rest of the input follows the usual format
```

When `dof_per_node` is set to `2`, any boundary condition `dof` value such
as `X` or `Y` is mapped to local indices `0` and `1` respectively. If a
`dof` provided in the YAML is out of range for the declared `dof_per_node`
the preprocessor raises a `ValueError`.

If you want full vector coupling (i.e. stiffness that couples `X` and `Y`
DOFs), let me know and I can implement the Option B vectorial element
formulation — it requires changing element shape functions, constitutive
coupling, and assembly for a genuinely multi-dimensional element.

---

See also: `docs/YAML_User_Input_Spec_README.md` for the full YAML
format and `src/wundy/first.py` for the implementation details.

## Example input file

An example YAML input that demonstrates material definitions, element
blocks, distributed loads and boundary conditions is provided at
`docs/examples/material_input_example.yaml`. Use `wundy.ui.load` and
`wundy.ui.preprocess` to convert the YAML into the Python structures
expected by `first.first_fe_code()`.
