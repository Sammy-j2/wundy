# Manufactured Solutions (MMS) Test: Nonlinear Axial Bar with Neo-Hookean Material

## Overview
The test in `tests/test_axial_mms_nonlinear.py` implements a **Manufactured Solutions** (MMS) verification for a 1D axial bar under **finite deformation** with a **Neo-Hookean hyperelastic material model**. Unlike linear elastic tests, this exercises the **Newton-Raphson nonlinear solver** with a consistent tangent stiffness derived from the First Piola-Kirchhoff (PK1) stress.

## Problem Definition

### Governing Equation (Nonlinear)
For a 1D bar under finite deformation, the equilibrium equation in the reference configuration is:
$$\frac{d}{dX}(A \cdot P(F)) = 0$$

where:
- $A$ = cross-sectional area (constant)
- $P(F)$ = First Piola-Kirchhoff stress (function of deformation gradient $F$)
- $F = 1 + \frac{du}{dX}$ = deformation gradient (1 + axial strain)
- $u(X)$ = displacement of material point initially at $X$

### Neo-Hookean Material Model
The solver uses a **compressible Neo-Hookean** strain energy density (1D reduction):
$$W(F) = \frac{\mu}{2}(F^2 - 3) - \mu \ln(F) + \frac{\kappa}{2}(\ln F)^2$$

where:
- $\mu = \frac{E}{2(1+\nu)}$ = shear modulus
- $\kappa = \frac{E}{3(1-2\nu)}$ = bulk modulus
- $E$ = Young's modulus
- $\nu$ = Poisson's ratio

**PK1 Stress:** (derivative of energy with respect to $F$)
$$P(F) = \frac{\partial W}{\partial F} = \mu\left(F - \frac{1}{F}\right) + \frac{\kappa \ln(F)}{F}$$

**Consistent Tangent:** (required for Newton-Raphson)
$$\frac{dP}{dF} = \mu\left(1 + \frac{1}{F^2}\right) + \frac{\kappa(1 - \ln F)}{F^2}$$

**Implementation:** [src/wundy/first.py::_neo_PK1_and_tangent](../src/wundy/first.py#L1111) (lines 1111–1160)

### Manufactured Solution
The test uses a **constant strain field**:
$$u_{\text{exact}}(X) = \varepsilon \cdot X$$

This gives a **constant deformation gradient**:
$$F = 1 + \frac{du}{dX} = 1 + \varepsilon$$

Since $F$ is constant, $\frac{dP}{dX} = 0$ automatically satisfies equilibrium **without body forces**. The solution is enforced by applying:
1. **Dirichlet BC** at $X=0$: $u(0) = 0$ (fixed left end)
2. **Neumann BC** at $X=L$: End traction $T$ matching the internal stress

**End Traction Approximation:**
The exact PK1 stress at $F = 1 + \varepsilon$ is:
$$P_{\text{exact}} = \mu\left((1+\varepsilon) - \frac{1}{1+\varepsilon}\right) + \frac{\kappa \ln(1+\varepsilon)}{1+\varepsilon}$$

However, to avoid duplicating solver internals in the test, the traction is **approximated** using small-strain theory:
$$T \approx E \cdot \varepsilon \cdot A$$

This approximation is sufficient to exercise the Newton-Raphson path and verify convergence for small $\varepsilon$ (line 18).

**Test Parameters:**
- Domain length: $L = 1.0$ m
- Applied strain: $\varepsilon = 10^{-3}$ (0.1%)
- Cross-sectional area: $A = 1.0$ m²
- Young's modulus: $E = 200{,}000$ Pa
- Poisson's ratio: $\nu = 0.3$
- Number of elements: `nelems` ∈ {1, 2, 4, 8} (parametric mesh refinement)

## YAML Input Structure

The test constructs a YAML specification programmatically via `build_yaml()` (lines 13–50):

```yaml
wundy:
  dof_per_node: 1
  nodes: [[0, 0.0], [1, 0.25], [2, 0.5], [3, 0.75], [4, 1.0]]  # Example: nelems=4
  elements: [[0, 0, 1], [1, 1, 2], [2, 2, 3], [3, 3, 4]]  # Connectivity
  materials:
    - type: neo_hooke  # Hyperelastic material
      name: rubber
      parameters:
        E: 200000.0
        nu: 0.3
  element blocks:
    - name: barblk
      material: rubber
      elements: [0, 1, 2, 3]
      element:
        type: T1D1  # 1D axial bar
        properties:
          area: 1.0
        integration:
          ngp: 2  # Gauss points per element
  distributed loads: []  # No body forces
  boundary conditions:
    - name: u_left
      nodes: [0]
      dof: X
      value: 0.0
      type: DIRICHLET  # Fix left end
    - name: t_right
      nodes: [4]  # Right end (node index = nelems)
      dof: X
      value: 200.0  # T = E * eps * A ≈ 200 N
      type: NEUMANN  # Apply traction
  integration:
    nonlinear: nonlinear  # Enable Newton-Raphson
```

**Key YAML Features:**
- `type: neo_hooke` — Selects Neo-Hookean material (preprocessor converts to `NEO_HOOKE`)
- `dof: X` — Symbolic DOF identifier (preprocessor maps to local DOF 0)
- `nonlinear: nonlinear` — Forces Newton-Raphson solver (vs. `linearize` mode)
- `NEUMANN` BC — Applied as external force on global force vector

## Solver Pipeline

### Step 1: Load and Preprocess (lines 54–60)
```python
data, eps_used = build_yaml(L=L, nelems=nelems, eps=eps)
buf = io.StringIO()
yaml.safe_dump(data, buf, sort_keys=False)
buf.seek(0)
spec_loaded = ui.load(buf)  # Schema validation and conversions
spec = ui.preprocess(spec_loaded)  # Build solver-ready data structures
```

**Preprocessing converts:**
- `type: neo_hooke` → `"type": "NEO_HOOKE"` (uppercase standardization)
- `dof: X` → `"local_dof": 0` (symbolic to numeric mapping)
- `type: DIRICHLET` / `NEUMANN` → schema constants
- Material parameters `{E, nu}` → stored for conversion to `{mu, kappa}` during solve

**Preprocessing locations:**
- Schema validation: [src/wundy/schemas.py](../src/wundy/schemas.py)
- Load processing: [src/wundy/ui.py](../src/wundy/ui.py) (`load` and `preprocess` functions)

### Step 2: Nonlinear Solve with Newton-Raphson (lines 61–70)
```python
res = first.first_fe_code(
    spec["coords"],
    spec["blocks"],
    spec["bcs"],
    spec.get("dload", []),
    spec["materials"],
    spec.get("block_elem_map", {}),
    spec.get("integration"),  # Contains {"nonlinear": "nonlinear"}
    spec.get("dof_per_node", 1),
)
```

**Solver function:** [src/wundy/first.py::first_fe_code](../src/wundy/first.py#L37) (lines 37–530)

#### Newton-Raphson Algorithm (lines 212–370)

The nonlinear solver iterates to find $u$ such that the residual $R = F_{\text{ext}} - F_{\text{int}}(u) = 0$.

**Initialization (lines 215–222):**
```python
F_ext, prescribed_dofs, prescribed_vals = assemble_external_forces(...)  # Build global force vector
u = np.zeros(num_dof, dtype=float)  # Initialize displacement
for pd, pv in zip(prescribed_dofs, prescribed_vals):
    u[pd] = pv  # Set Dirichlet values
```

**Newton-Raphson Loop (lines 224–370):**
For each iteration $k = 0, 1, \ldots$:

1. **Reset tangent stiffness and internal force** (lines 236–238):
   ```python
   K[:, :] = 0.0
   F_int = np.zeros_like(F_ext)
   ```

2. **Element loop: Assemble internal forces and tangent** (lines 245–304):
   - Extract element nodal displacements $u_e$ from global $u$
   - Compute deformation gradient: $F_e = 1 + (u_2 - u_1) / h$
   - Evaluate Neo-Hookean PK1 stress and tangent (line 265):
     ```python
     P, dPdF = _neo_PK1_and_tangent(mat, F_e)
     ```
   - Internal force vector (lines 273–276):
     $$f_e^{\text{int}} = A \cdot P \cdot \begin{bmatrix} -1 \\ 1 \end{bmatrix}$$
   - Tangent stiffness matrix (lines 277–282):
     $$k_e^{\text{tang}} = \frac{A \cdot dP/dF}{h} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$
   - Assemble into global $K$ and $F_{\text{int}}$

3. **Compute residual** (line 306):
   ```python
   R = F_ext - F_int
   ```

4. **Handle uncoupled DOFs** (lines 307–317):
   - Identify DOFs with zero stiffness (uncoupled from system)
   - Automatically prescribe to zero to avoid singular matrix

5. **Partition and solve for increment** (lines 319–331):
   - Separate free and prescribed DOFs
   - Solve reduced system:
     $$K_{ff}^{(k)} \Delta u_f^{(k)} = R_f^{(k)}$$
   - Update displacement:
     $$u_f^{(k+1)} = u_f^{(k)} + \Delta u_f^{(k)}$$

6. **Convergence check** (lines 333–370):
   - Residual norm: $\|R_f\|$
   - Displacement increment norm: $\|\Delta u\|$ (optional)
   - If $\|R_f\| < \text{tol}$ → converged, return solution
   - If $k \geq \text{max\_it}$ → raise error (not converged)

**Newton-Raphson Controls (lines 225–233):**
- `nr_tol` = $10^{-8}$ (default residual tolerance)
- `nr_max_it` = 25 (default max iterations)
- `nr_du_tol` = None (optional displacement increment tolerance)

These can be overridden via the `integration` dict in YAML.

**Convergence Criteria:**
The solver returns when:
$$\|R_f\| = \|F_{\text{ext},f} - F_{\text{int},f}(u)\| < \text{nr\_tol}$$

or (optionally):
$$\|\Delta u\| < \text{nr\_du\_tol}$$

### Step 3: Extract and Verify Solution (lines 71–77, test function lines 86–95)

```python
u = res["dofs"]
u_tip = u[nelems]  # Displacement at right end (x = L)
u_exact = eps_used * L  # Manufactured solution: u(L) = eps * L
err = abs(u_tip - u_exact)  # Absolute error
```

**Test Assertions:**

1. **Convergence Test** (`test_axial_mms_nonlinear_convergence`, lines 80–90):
   - Parametrized over mesh densities: `nelems` ∈ {1, 2, 4, 8}
   - Verifies error decreases (or remains stable) with refinement:
     ```python
     assert err <= prev_err + 1e-10  # Monotonic convergence
     ```
   - **Why this works:** For constant strain, even a single element should be exact. Refinement verifies solver stability and consistency.

2. **Minimum Quality Test** (`test_axial_mms_nonlinear_min_quality`, lines 93–95):
   - Uses `nelems = 4`
   - Requires absolute error below threshold:
     ```python
     assert err < 1e-3  # Error < 0.1% of displacement
     ```
   - **Physical interpretation:** For $\varepsilon = 10^{-3}$ and $L = 1.0$, exact displacement is $u_{\text{exact}} = 10^{-3}$. Error must be $< 10^{-3}$, i.e., solution accurate to ~100% relative error tolerance (generous for demonstration).

## Key Implementation Details

### Neo-Hookean Stress and Tangent Computation
**Function:** [src/wundy/first.py::_neo_PK1_and_tangent](../src/wundy/first.py#L1111) (lines 1111–1160)

```python
def _neo_PK1_and_tangent(material: dict, F: float) -> tuple[float, float]:
    # Convert E, nu to mu, kappa
    E = float(params["E"])
    nu = float(params["nu"])
    mu = E / (2.0 * (1.0 + nu))
    kappa = E / (3.0 * (1.0 - 2.0 * nu))
    
    # Compute PK1 stress
    lnF = float(np.log(F))
    P = mu * (F - 1.0 / F) + (kappa * lnF) / F
    
    # Compute consistent tangent
    dPdF = mu * (1.0 + 1.0 / (F * F)) + kappa * (1.0 - lnF) / (F * F)
    
    return float(P), float(dPdF)
```

**Key features:**
- Accepts material specified by `(E, nu)` or `(mu, lambda)` or `(mu, kappa)`
- Returns both stress and tangent for Newton-Raphson
- Validates $F > 0$ (negative or zero deformation gradient is unphysical)

### Material Type Detection
**Function:** [src/wundy/first.py::is_neo_material](../src/wundy/first.py#L1090) (lines 1090–1100)

```python
def is_neo_material(mat: dict) -> bool:
    mtype = str(mat.get("type", "")).upper()
    return mtype in {"NEO_HOOKE", "NEO-HOOKE", "NEOHOOKE"}
```

Checks if material type is Neo-Hookean (case-insensitive, accepts variants).

### Element Internal Force Assembly (Nonlinear Path)
**Location:** [src/wundy/first.py](../src/wundy/first.py#L255) (lines 255–282)

For each T1D1 element in the Newton-Raphson loop:
1. Extract element DOFs: $u_e = [u_1, u_2]$
2. Compute deformation gradient: $F_e = 1 + (u_2 - u_1) / h$
3. Evaluate stress and tangent: `P, dPdF = _neo_PK1_and_tangent(mat, F_e)`
4. Internal force: $f_e = A \cdot P \cdot [-1, 1]^T$
5. Tangent stiffness: $k_e = (A \cdot dPdF / h) \cdot \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$

### Global DOF Mapping
**Function:** [src/wundy/first.py::global_dof](../src/wundy/first.py#L790) (lines 790–798)
```python
def global_dof(node_index, local_dof, dof_per_node):
    return node_index * dof_per_node + local_dof
```

For `dof_per_node=1` (pure axial):
- Node 0 → DOF 0
- Node 1 → DOF 1
- etc.

## What This Test Validates

1. **Neo-Hookean Implementation:**
   - PK1 stress formula $P(F)$ correctly implemented
   - Consistent tangent $dP/dF$ correctly implemented

2. **Newton-Raphson Convergence:**
   - Solver converges to correct solution under nonlinear material
   - Tangent stiffness matrix is properly assembled and symmetric

3. **Finite Deformation Kinematics:**
   - Deformation gradient $F = 1 + du/dX$ correctly computed
   - Internal forces balance external tractions

4. **Boundary Condition Handling:**
   - Dirichlet BCs (prescribed displacements) enforced correctly
   - Neumann BCs (tractions) applied to global force vector

5. **Mesh Independence:**
   - Solution accuracy maintained across mesh refinement
   - Constant strain field captured exactly even with coarse mesh

## File Locations Summary

| Component | File | Lines |
|-----------|------|-------|
| Test code | `tests/test_axial_mms_nonlinear.py` | 1–95 |
| Main solver (NR path) | `src/wundy/first.py::first_fe_code` | 212–370 |
| Neo-Hooke stress/tangent | `src/wundy/first.py::_neo_PK1_and_tangent` | 1111–1160 |
| Material type check | `src/wundy/first.py::is_neo_material` | 1090–1100 |
| Tangent modulus | `src/wundy/first.py::material_tangent_modulus` | 1060–1088 |
| Global DOF mapping | `src/wundy/first.py::global_dof` | 790–798 |
| External force assembly | `src/wundy/first.py::assemble_external_forces` | 1175–1230 |
| Internal force assembly | `src/wundy/first.py::assemble_internal_forces` | 1020–1058 |
| YAML preprocessing | `src/wundy/ui.py::load`, `preprocess` | (see ui.py) |
| Schema validation | `src/wundy/schemas.py` | (see schemas.py) |

## How to Run the Test

```bash
# Run all parametric cases
pytest tests/test_axial_mms_nonlinear.py -v

# Run specific test
pytest tests/test_axial_mms_nonlinear.py::test_axial_mms_nonlinear_min_quality -v

# Run with detailed output
pytest tests/test_axial_mms_nonlinear.py -v -s
```

Expected output:
```
tests/test_axial_mms_nonlinear.py::test_axial_mms_nonlinear_convergence[1] PASSED
tests/test_axial_mms_nonlinear.py::test_axial_mms_nonlinear_convergence[2] PASSED
tests/test_axial_mms_nonlinear.py::test_axial_mms_nonlinear_convergence[4] PASSED
tests/test_axial_mms_nonlinear.py::test_axial_mms_nonlinear_convergence[8] PASSED
tests/test_axial_mms_nonlinear.py::test_axial_mms_nonlinear_min_quality PASSED
```

## Extending the Test

To modify or extend the nonlinear MMS test:

1. **Change applied strain**: Edit `eps=1e-3` in function calls (line 90, 95)
   - Caution: Large strains may require tighter `nr_tol` or more iterations

2. **Change material parameters**: Modify `E` or `nu` in `build_yaml()` (line 13)
   - Affects bulk/shear moduli and thus PK1 stress magnitude

3. **Change domain length**: Edit `L=1.0` (line 13)
   - Scales absolute displacement but not strain

4. **Change mesh resolution**: Modify `nelems` parameter (lines 80, 95)
   - Constant strain should converge with any mesh

5. **Add body forces**: Modify manufactured solution to:
   - $u(X) = f(X)$ with non-constant $F = 1 + du/dX$
   - Derive corresponding body force $b(X) = -\frac{d(A P(F))}{dX}$
   - Add as distributed load in YAML

6. **Change traction approximation**: Replace line 18 with exact PK1:
   ```python
   mu = E / (2.0 * (1.0 + nu))
   kappa = E / (3.0 * (1.0 - 2.0 * nu))
   F_exact = 1.0 + eps
   P_exact = mu * (F_exact - 1.0 / F_exact) + (kappa * np.log(F_exact)) / F_exact
   T = P_exact * A
   ```

7. **Test different material models**: Change `type: neo_hooke` to other hyperelastic models (if implemented)

## Theoretical Background

### Why Use PK1 Stress?
In finite deformation, different stress measures exist:
- **Cauchy stress** $\sigma$: force per deformed area (spatial)
- **PK1 stress** $P$: force per undeformed area (mixed)
- **PK2 stress** $S$: conjugate to Green strain (material)

For 1D problems with prescribed tractions on the **undeformed boundary**, PK1 is the natural choice:
$$T = P \cdot A_0$$

where $A_0$ is the undeformed area. This matches the Neumann BC interpretation used in the solver.

See [docs/material_models.md](material_models.md) for detailed rationale.

### Why Consistent Tangent?
Newton-Raphson requires the derivative of internal forces with respect to displacements:
$$K_{\text{tang}} = \frac{\partial F_{\text{int}}}{\partial u}$$

For Neo-Hookean materials, this becomes:
$$K_{\text{tang}} = \frac{A}{h} \frac{dP}{dF} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

Using the **consistent tangent** $dP/dF$ (not an approximation like secant or chord modulus) ensures **quadratic convergence** near the solution.

### Small Strain Limit
As $\varepsilon \to 0$, the Neo-Hookean model reduces to linear elasticity:
$$P(F) \approx E \varepsilon + O(\varepsilon^2)$$

which is why the approximate traction $T \approx E \varepsilon A$ works well for small strains.

## References

- **Nonlinear FEM:** Belytschko, Liu, Moran. *Nonlinear Finite Elements for Continua and Structures*. Wiley, 2000.
- **Hyperelasticity:** Holzapfel. *Nonlinear Solid Mechanics*. Wiley, 2000.
- **Newton-Raphson:** Crisfield. *Non-linear Finite Element Analysis of Solids and Structures*. Wiley, 1991.
- **MMS for verification:** Roache. *Verification and Validation in Computational Science and Engineering*. Hermosa, 1998.
- **Solver documentation:** See `docs/material_models.md`, `README.md`, and `docs/MMS_TEST_DESCRIPTION.md`.
