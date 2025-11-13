# wundy
One dimension finite element program

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

## Material models

See `docs/material_models.md` for details on how to specify materials in
YAML or programmatically, what keys are required, and examples showing how
the solver uses `E` and `density`.

Note: the solver supports compressible Neo‑Hookean materials. To use the
nonlinear Neo‑Hookean path set the material `type` to any token containing
`neo` and `hook` (for example `neo_hooke` or `neo-hooke`) and supply
parameters either as `(E, nu)` or as Lame-style `(mu, lambda)` or `(mu, kappa)`.
When Neo‑Hooke is present the solver performs a Newton–Raphson solve with
default tolerance 1e-8 and max iterations 25 and returns the nonlinear
nodal displacements. See `docs/material_models.md` for examples.

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
