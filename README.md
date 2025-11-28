# wundy
One-dimension finite element program for my Theory of FEM course — now with axial bars and Euler–Bernoulli beams.

## Install

### Clone repository

```console
git clone git@github.com:<user>/wundy
```

where `user` is your user name if you forked the repo, `tjfulle` otherwise.

### Create virtual environment

```console
python3 -m venv venv
source activate venv/bin/activate
```

### Install in editable mode

```console
cd wundy
python3 -m pip install -e .
```

## Test

In the `wundy` directory, execute

```console
pytest
```

## Material models and nonlinear solver

See `docs/material_models.md` for details on how to specify materials in
YAML or programmatically, what keys are required, and examples showing how
the solver uses `E` and `density`.

The solver supports compressible Neo‑Hookean materials for axial deformation.
Set the material `type` to include both `neo` and `hook` (e.g., `neo_hooke` or
`neo-hooke`) and supply parameters either as `(E, nu)` or Lame-style `(mu, lambda)`
or `(mu, kappa)`. When Neo‑Hooke is present and `integration.nonlinear: nonlinear`
is selected, the solver runs a Newton–Raphson path using PK1 stress with a consistent
tangent for axial elements. Current limitation: EB1D bending is treated linearly
with tangent `E` (no geometric nonlinearity); this allows mixed axial+beam problems
to run NR stably while beams use linear bending.

## Beam elements (EB1D)

An Euler–Bernoulli beam element type `EB1D` is available for simple 1‑D bending
problems. Each node must have at least two DOFs (`dof_per_node: 2`) interpreted as
`X -> w` (transverse displacement) and `Y -> theta` (rotation). The element stiffness
matrix (Hermite cubic) is assembled with ordering `[w1, theta1, w2, theta2]`.

Element block example:

```yaml
element blocks:
	- name: beamblk
		material: rubber
		elements: [0]
		element:
			type: EB1D
			properties:
				I: 1.0e-6   # second moment of area
				area: 1.0   # optional, used by body/gravity loads
			integration:
				ngp: 3
```

Distributed transverse loads use types `QY` or `QBEAM` and can specify one of:

- `value: <constant>` (uniform load)
- `table: [[x,q], ...]` (piecewise linear interpolation)
- `expression: "100.0 * x / L"` (evaluated safely with variables `x` and `L`)
- `value: !!python/name:module.func` (callable in programmatic usage)

Examples are in tests (`tests/test_beam_yaml.py`) and docs (`docs/examples`).

Neo‑Hookean bending currently uses the linear tangent modulus `E` (no geometric
nonlinearity for rotations). If `integration.nonlinear: nonlinear` is set, axial elements
participate in NR while EB1D contributions remain linear.

## Example: run the material input

An example YAML input is included at `docs/examples/material_input_example.yaml`.
To run it from PowerShell using the project virtual environment (assumes
the venv is named ` .venvwundy` in the project root), activate the venv
and execute a short Python snippet to load and run the solver:

```powershell
# activate venv (PowerShell)
.\.venvwundy\Scripts\Activate.ps1

# run the example (this prints the global DOF vector)
python - <<'PY'
from pathlib import Path
from wundy import ui, first

p = Path('docs/examples/material_input_example.yaml')
with p.open() as fh:
	data = ui.load(fh)
inp = ui.preprocess(data)
sol = first.first_fe_code(
	inp['coords'], inp['blocks'], inp['bcs'], inp['dload'], inp['materials'], inp['block_elem_map']
)
print('dofs:', sol['dofs'])
PY
```

Alternatively, copy the snippet into a small script (for example
`bin/run_example.py`) and run it with `python bin/run_example.py`.

Alternatively, use the included runner script which accepts a YAML path
or opens a GUI file picker on Windows:

```powershell
# run by path
python .\bin\run_yaml.py docs\examples\material_input_example.yaml

# open file picker (Windows)
python .\bin\run_yaml.py --pick
```

## Example: multi-DOF input and pipeline

An example showing `dof_per_node: 2` and how to specify `X`/`Y` DOFs is
included at `docs/examples/multidof_example.yaml`. This demonstrates the
preprocessor mapping of symbolic DOF names to local indices and the
axial-only assembly behaviour documented in `docs/material_models.md`.

To run the multi-DOF example:

```powershell
# activate venv (PowerShell)
.\.venvwundy\Scripts\Activate.ps1

python - <<'PY'
from pathlib import Path
from wundy import ui, first

p = Path('docs/examples/multidof_example.yaml')
with p.open() as fh:
	data = ui.load(fh)
inp = ui.preprocess(data)
sol = first.first_fe_code(
	inp['coords'], inp['blocks'], inp['bcs'], inp['dload'], inp['materials'], inp['block_elem_map'], dof_per_node=inp['dof_per_node']
)
print('dofs:', sol['dofs'])
print('residual:', sol['residual'])
PY
```

Runner note: if you enable nonlinear Newton–Raphson (see docs/material_models.md)
you can control Newton behavior with the `integration` block in your YAML:

```yaml
wundy:
	integration:
		nonlinear: nonlinear    # 'linearize'|'nonlinear'|'auto'
		internal: gauss         # 'analytic'|'gauss'
		ngp: 2                 # number of Gauss points when using 'gauss'
		nr_tol: 1e-8           # residual tolerance for Newton
		nr_du_tol: 1e-10       # optional displacement-increment tolerance
		nr_max_it: 50         # maximum Newton iterations
```

See `docs/examples/newton_example.yaml` for a working Newton example.

Manufactured solutions (MMS):
- Beam MMS runs with `integration.nonlinear: linearize` to avoid singular NR tangents in pure EB1D cases; see `tests/test_beam_mms.py`.
- Axial MMS runs fully nonlinear with Neo‑Hooke; see `tests/test_axial_mms_nonlinear.py`.

## Design Rationale

- PK1 stress for axial Neo‑Hooke: Work-conjugate in the reference configuration; gives simple internal forces `f_e = A P [-1, 1]^T` and a consistent tangent for robust NR.
- EB1D linearization: Current beam bending uses tangent `E` without geometric nonlinearity. This keeps mixed axial+beam problems stable under NR while we focus nonlinearity on axial behavior.
- Automatic zero-prescription: With `dof_per_node > 1`, non-axial DOFs can be uncoupled for axial-only elements. We automatically prescribe any global DOF with a zero stiffness row to `0.0` to avoid singular systems. Use genuinely vectorial elements if you need coupled multi-DOF mechanics.

Beam Newton example (uniform load cantilever):

```yaml
wundy:
	dof_per_node: 2
	nodes:
		- [0, 0.0]
		- [1, 1.0]
	elements:
		- [0, 0, 1]
	materials:
		- name: rubber
			type: neo_hooke
			parameters:
				E: 200000.0
				nu: 0.3
	element blocks:
		- name: beamblk
			material: rubber
			elements: [0]
			element:
				type: EB1D
				properties:
					I: 1.0e-6
					area: 1.0
				integration:
					ngp: 3
	distributed loads:
		- name: q_uniform
			elements: [0]
			type: QY
			value: 100.0
			direction: [1]
	boundary conditions:
		- {name: clampw, nodes: [0], dof: X, value: 0.0, type: DIRICHLET}
		- {name: clamptheta, nodes: [0], dof: Y, value: 0.0, type: DIRICHLET}
	integration:
		nonlinear: nonlinear
```

Run with:

```powershell
python .\bin\run_yaml.py docs\examples\beam_uniform.yaml
```

Run the test suite with:

```powershell
pytest -q
```
