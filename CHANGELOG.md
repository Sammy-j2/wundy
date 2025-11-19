# Changelog

All notable changes to this project are documented in this file.

## 2025-11-19 — Multi-DOF, materials, integration, docs, tests

Implemented features (user-requested during development session):

- Add Neo‑Hookean material support (1D compressible PK1 and consistent tangent).
  - Material input accepts `(E, nu)` or Lame-style `(mu, lambda)` / `(mu, kappa)`.
  - Solver can run in a full Newton–Raphson path when `integration.nonlinear: nonlinear`.

- YAML-driven integration options:
  - Global and per-block `integration` settings (`stiffness`, `internal`, `ngp`, `nonlinear`).
  - Gauss quadrature selectable via `internal: gauss` and per-block `ngp` overrides.

- Added `dof_per_node` (uniform DOFs per node) support and plumbing:
  - Solver and preprocessor accept `dof_per_node` (e.g. `2` for X/Y per node).
  - Implementation uses axial-only coupling (Option A): axial DOF is local index 0.
  - Non-axial DOFs are currently uncoupled (receive zero internal force/stiffness).
  - Solver auto-prescribes fully-zero global DOFs to `0.0` to avoid singular K.

- Tests and docs:
  - Added unit test `tests/test_multidof_solver.py` validating multi-DOF axial behavior.
  - Updated `docs/material_models.md` with multi-DOF notes and Neo‑Hooke examples.
  - Added `docs/examples/multidof_example.yaml` and updated `README.md` with usage snippets.
  - Added inline module docstrings and comments in `src/wundy/first.py`, `src/wundy/ui.py`, and `src/wundy/schemas.py`.

Notes / next steps (planned):

- Implement nonlinear Gauss integration for consistent tangent integration in NR path.
- Add robust Newton solver features: load-stepping and line-search.
- Optionally implement full vectorial element coupling (Option B) for true multi-DOF coupling.

For more details, see `docs/material_models.md` and the unit tests in `tests/`.
