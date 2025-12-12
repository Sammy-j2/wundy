# Manufactured Solutions (MMS) Test: Axial Bar with Sinusoidal Distributed Load

## Overview
The test in `tests/testMMSinCLass.py` implements a **Manufactured Solutions** (MMS) verification for a 1D elastic axial bar subjected to a distributed load. This test validates that the solver correctly computes nodal displacements and internal forces under a prescribed analytical load distribution.

## Problem Definition

### Governing Equation
For a 1D axial bar under distributed load $q(x)$, the strong form is:
$$\frac{d}{dx}\left(A E \frac{du}{dx}\right) + q(x) = 0$$

with:
- $A$ = cross-sectional area
- $E$ = Young's modulus  
- $u(x)$ = axial displacement
- $q(x)$ = distributed load per unit length

### Manufactured Solution
The test uses a **sinusoidal manufactured solution**:
$$u_{\text{exact}}(x) = \sin(\pi x)$$

To generate the corresponding load, we substitute into the governing equation:
$$q(x) = -\frac{d}{dx}\left(A E \frac{du}{dx}\right) = -A E \pi^2 \sin(\pi x)$$

With material properties $A=1, E=10$:
$$q(x) = -10\pi^2 \sin(\pi x)$$

**Test Parameters:**
- Domain length: $L = 2$
- Number of elements: $\text{num\_elems} = 50$
- Cross-sectional area: $A = 1.0$
- Young's modulus: $E = 10.0$
- Boundary conditions: **Fixed at both ends** ($u(0)=0$ and $u(L)=0$, nodes 1 and 51)
- Distributed load: $q(x) = 10\pi^2 \sin(\pi x)$ (expression-based, includes E factor)

## YAML Input Structure

The test constructs a YAML specification programmatically (lines 15-47 in `testMMSinCLass.py`):

```yaml
wundy:
  nodes: [[1, 0.0], [2, 0.04], ..., [51, 2.0]]  # num_elems+1 nodes uniformly spaced
  elements: [[1, 1, 2], [2, 2, 3], ..., [50, 50, 51]]  # Element connectivity
  boundary conditions:
    - name: fix-nodes
      dof: x
      nodes: [1]  # Fixed at first node (x=0)
    - name: fix-right
      dof: x
      nodes: [51]  # Fixed at last node (x=2)
  materials:
    - type: elastic
      name: mat-1
      parameters:
        E: 10.0
        nu: 0.3
  element blocks:
    - material: mat-1
      name: block-1
      elements: all
      element:
        type: t1d1  # 1D bar element
        properties:
          area: 1
  distributed loads:
    - name: dload-1
      elements: all
      type: bx  # Body load (axial)
      expression: "10*pi**2*sin(pi * x)"  # q(x) = 10π²sin(πx) with E factor
      direction: [1]  # Positive x-direction
```

## Solver Pipeline

### Step 1: Load and Preprocess (lines 48-50)
```python
data = wundy.ui.load(file)          # Parse YAML using schemas.py
inp = wundy.ui.preprocess(data)     # Convert YAML to solver-ready format
```

**Preprocessing converts:**
- Node list → coordinate array `coords` (shape: `num_nodes × 1`)
- Element list → connectivity array `blocks[0]["connect"]`
- Boundary conditions → list of dicts with Dirichlet prescriptions
- Distributed loads → list of dicts with expression evaluation context
- Material specs → material database

**Key preprocessing locations:**
- Schema validation: [src/wundy/schemas.py](../src/wundy/schemas.py) (defines allowed YAML keys and types)
- Load processing: [src/wundy/ui.py](../src/wundy/ui.py) (loads YAML, applies schema, builds preprocessed dict)

### Step 2: Solve Using `first_fe_code` (lines 51-57)
```python
soln = wundy.first.first_fe_code(
    inp["coords"],           # Node coordinates
    inp["blocks"],           # Element blocks with connectivity
    inp["bcs"],              # Boundary conditions
    inp["dload"],            # Distributed loads
    inp["materials"],        # Material properties
    inp["block_elem_map"],   # Maps global element IDs to (block_idx, local_idx)
)
```

**Solver function:** [src/wundy/first.py::first_fe_code](../src/wundy/first.py#L37) (lines 37–530)

#### Solver Algorithm Overview

The `first_fe_code` function performs the following:

1. **Initialize global stiffness matrix K and force vector F** (lines 160–180)
   - `num_dof = num_nodes * dof_per_node`
   - `K` shape: `(num_dof, num_dof)`, initialized to zeros
   - `F` shape: `(num_dof,)`, initialized to zeros

2. **Assemble element stiffness matrices** (lines 370–430)
   - For each element block, compute element stiffness $k_e$ and add to global $K$
   - For T1D1 (axial bar) elements:
     $$k_e = \frac{A E}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$
   - **Element stiffness function:** [src/wundy/first.py::element_stiffness](../src/wundy/first.py#L860) (lines 860–900)
   - **Assembly:** Direct stiffness method, place $k_e$ into global $K$ at DOF indices (line 428)

3. **Apply Boundary Conditions** (lines 432–434)
   - Identify which DOFs are prescribed (Dirichlet) and their values
   - Store in `prescribed_dofs` and `prescribed_vals`
   - **BC processing function:** [src/wundy/first.py::apply_boundary_conditions](../src/wundy/first.py#L500) (lines 500–545)

4. **Assemble Distributed Loads** (lines 436–438)
   - Add equivalent nodal forces from distributed load $q(x)$ to global $F$
   - **Load assembly function:** [src/wundy/first.py::apply_distributed_loads](../src/wundy/first.py#L590) (lines 590–760)
   
   #### Distributed Load Assembly Details
   
   For a distributed load of type "BX" (body load) with expression $q(x) = -\pi^2 \sin(\pi x)$:
   
   **Process for each element:**
   - Extract element nodes, coordinates, and length
   - If load is expression-based (line 626):
     ```python
     dload["value"] = lambda x: eval(dload["expression"], {"sin": np.sin, "pi": np.pi, "x": x})
     ```
   - Evaluate load at Gauss quadrature points (line 742):
     ```python
     xi_pts, wts = _gauss_1d(ngp_elem)  # Get Gauss points and weights
     ```
   - For each quadrature point, compute:
     - Physical coordinate: $x = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2$
     - Load value: $q(x)$
     - Shape functions: $N_1 = \frac{x_2 - x}{h}$, $N_2 = \frac{x - x_1}{h}$
     - Element equivalent loads: $f_e = \sum_{qp} w_{qp} N(x_{qp}) q(x_{qp}) J$
   - **Gauss integration function:** [src/wundy/first.py::_gauss_1d](../src/wundy/first.py#L830) (lines 830–844)
   - **Shape functions:** [src/wundy/first.py::_shape_funcs](../src/wundy/first.py#L800) (lines 800–810)

5. **Handle Uncoupled DOFs** (lines 440–450)
   - Identify DOFs with all-zero stiffness rows (uncoupled from rest of system)
   - Automatically prescribe such DOFs to zero to avoid singular system
   - **Zero-row detection:** Lines 450–456

6. **Solve Reduced System** (lines 458–480)
   - Partition DOFs into free and prescribed
   - Build reduced system for free DOFs: $K_{ff} u_f = F_f - K_{fp} u_p$
   - Solve using `np.linalg.solve(Kff, Rf)` (line 472)
   - Back-substitute to get full displacement vector

7. **Return Solution** (lines 482–490)
   - `dofs`: Solved nodal displacement vector $u$
   - `stiff`: Global stiffness matrix $K$
   - `force`: Global force vector $F$

### Step 3: Verification (lines 58–63)

The test verifies solution correctness by checking the **residual**:

```python
dofs = soln["dofs"]
K = soln["stiff"]
F = soln["force"]
R = np.dot(K, dofs) - F  # Residual: should be zero at free DOFs

# Verify residual at prescribed DOF boundary
assert np.allclose(R[0], -q * L)
```

**What this tests:**
- $K u = F$ is satisfied at the algebraic level
- The reaction force at the fixed node equals the integral of the distributed load: $\int_0^L q(x) \, dx$
- For sinusoidal load $q(x) = -\pi^2 \sin(\pi x)$: $\int_0^2 q(x) \, dx = -\pi^2 \sin(\pi) = 0$ (but test adjusts expectation)

## Key Implementation Details

### Global Degree of Freedom Indexing
The solver maps local element DOFs to global DOF numbers using:

**Function:** [src/wundy/first.py::global_dof](../src/wundy/first.py#L790) (lines 790–798)
```python
def global_dof(node_index, local_dof, dof_per_node):
    """Map (node, local_DOF) to global DOF index."""
    return node_index * dof_per_node + local_dof
```

For `dof_per_node=1` (1D axial):
- Node 0 → global DOF 0
- Node 1 → global DOF 1
- etc.

For `dof_per_node>1` (e.g., beams with DOF=[w, theta]):
- Node 0: global DOFs [0, 1]
- Node 1: global DOFs [2, 3]
- etc.

### Material Properties
Material database mapping is accessed by block reference:

**Line 368:**
```python
mat = materials[block["material"]]  # e.g., "MAT-1"
E = float(mat["parameters"]["E"])  # Young's modulus from params
```

### Integration Controls
Quadrature order and solver behavior are controlled via optional `integration` dict:

```python
integration = {
    "ngp": 2,                  # Number of Gauss points (default 2)
    "internal": "analytic",    # Method for internal forces (analytic|gauss)
    "stiffness": "analytic",   # Method for stiffness assembly (analytic|gauss)
    "nonlinear": "linearize"   # Nonlinear solver mode (linearize|nonlinear|auto)
}
```

**Default values set at lines 197–200.**

## Test Assertions

The test checks two key properties:

**Assertion 1 (Line 67-68):**
```python
total_reaction = R[0] + R[-1]
assert np.allclose(total_reaction, 0.0, atol=1e-10)
```

For the sinusoidal load $q(x) = 10\pi^2 \sin(\pi x)$, the integral over $[0, 2]$ is:
$$\int_0^2 10\pi^2 \sin(\pi x) \, dx = 10\pi^2 \left[-\frac{1}{\pi}\cos(\pi x)\right]_0^2 = 0$$

Therefore, reaction forces at both ends must sum to zero. This validates:
- **Load integration**: Distributed load expressions correctly evaluated and assembled
- **Force equilibrium**: Reaction forces preserve force balance

**Assertion 2 (Line 78):**
```python
assert np.allclose(u, ue, atol=1e-3, rtol=1e-2)
```

Where:
- `u` = computed nodal displacements
- `ue` = exact manufactured solution $u_{\text{exact}}(x) = \sin(\pi x)$

This validates:
1. **Stiffness assembly**: Global $K$ matrix correctly combines element matrices
2. **Boundary condition enforcement**: Fixed-fixed BCs correctly applied
3. **System solution**: Linear solver produces displacement field matching exact solution

4. **Boundary conditions**: Dirichlet constraints are properly enforced

## File Locations Summary

| Component | File | Lines |
|-----------|------|-------|
| Test code | `tests/testMMSinCLass.py` | 1–68 |
| Main solver | `src/wundy/first.py::first_fe_code` | 37–530 |
| Stiffness assembly | `src/wundy/first.py::element_stiffness` | 860–900 |
| Distributed loads | `src/wundy/first.py::apply_distributed_loads` | 590–760 |
| Gauss quadrature | `src/wundy/first.py::_gauss_1d` | 830–844 |
| Shape functions | `src/wundy/first.py::_shape_funcs` | 800–810 |
| BC application | `src/wundy/first.py::apply_boundary_conditions` | 500–545 |
| Global DOF mapping | `src/wundy/first.py::global_dof` | 790–798 |
| YAML preprocessing | `src/wundy/ui.py::load` | (see ui.py) |
| Schema validation | `src/wundy/schemas.py` | (see schemas.py) |

## How to Run the Test

```bash
pytest tests/testMMSinCLass.py -v
```

Expected output:
```
tests/testMMSinCLass.py::test_first_2 PASSED
```

## Extending the Test

To modify or extend the MMS test:

1. **Change domain length**: Edit `L = 2` (line 13)
2. **Change mesh resolution**: Edit `num_elems = 50` (line 14)
3. **Change manufactured solution**: Modify the expression at line 40 and corresponding load at line 12
4. **Change material properties**: Modify `E: 10.0` (line 33)
5. **Change element type**: Replace `type: t1d1` with another type (e.g., `eb1d` for beams)

## References

- **FEM theory**: Hughes, T.J.R. *The Finite Element Method*. Dover, 2000.
- **Manufactured solutions**: Oberkampf & Roy. *Verification and Validation in Scientific Computing*. Cambridge, 2010.
- **Solver documentation**: See `docs/material_models.md` for material model details and `README.md` for usage examples.
