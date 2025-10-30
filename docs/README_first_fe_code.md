# 1D Finite Element Solver — `first_fe_code`

This module implements a simple **1D finite element (FE)** solver for bar or truss-type problems using linear two-node elements. It constructs the global stiffness matrix and force vector, applies boundary conditions and loads, and solves for nodal displacements.

The function `first_fe_code()` performs the entire finite element assembly and solution sequence, while `global_dof()` provides consistent global degree of freedom (DOF) indexing across the mesh.

---

## Overview

In the finite element method, the physical domain (a bar, rod, etc.) is divided into smaller subdomains called **elements**. Each element has **nodes** where the displacement field is interpolated.

### Governing Equation

\[
\frac{d}{dx}(EAu') + q(x) = 0
\]

Applying the **Galerkin method** and assembling the system yields:

\[
K u = F
\]

where:

- \( K \) is the global stiffness matrix  
- \( u \) is the nodal displacement vector  
- \( F \) is the global load vector

---

## Function: `first_fe_code`

```python
first_fe_code(
    coords,
    blocks,
    bcs,
    dloads,
    materials,
    block_elem_map
)
```

### Purpose

Builds and solves a 1D FEM system using linear bar elements. The process includes:

1. **Assembling element stiffness matrices**
2. **Applying distributed and body loads**
3. **Applying Neumann and Dirichlet boundary conditions**
4. **Solving for nodal displacements**

### Typical Use Case

This function is useful for analyzing **axially loaded bars** or **simple truss members**, where each node carries a single translational degree of freedom.

### Data Flow

```
       Element Data (A, E)
             │
             ▼
   Local Stiffness Matrix ke
             │
             ▼
      Assembled into Global K
             │
             ▼
   Apply Loads and BCs → Solve → Displacements (u)
```

### Inputs

| Parameter | Type | Description |
|------------|------|-------------|
| `coords` | `ndarray (n_nodes, 1)` | Nodal coordinates along the x-axis |
| `blocks` | `list[dict]` | Element connectivity and material grouping |
| `bcs` | `list[dict]` | Boundary conditions (`DIRICHLET` or `NEUMANN`) |
| `dloads` | `list[dict]` | Distributed or gravitational loads |
| `materials` | `dict[str, Any]` | Material definitions (E, density, etc.) |
| `block_elem_map` | `dict[int, tuple[int, int]]` | Maps global element IDs to their block and local indices |

### Outputs

| Key | Description |
|-----|--------------|
| `"dofs"` | Computed nodal displacements |
| `"stiff"` | Assembled global stiffness matrix \( K \) |
| `"force"` | Global force vector \( F \) after loads and BCs |

---

## Function: `global_dof`

```python
global_dof(node, local_dof, dof_per_node)
```

### Purpose

Converts a **local** degree of freedom at a node into its **global** index.  
Used internally when assembling global stiffness and force matrices.

### Formula

\[
\text{global\_dof} = \text{node} \times \text{dof\_per\_node} + \text{local\_dof}
\]

### Example: DOF Mapping

For a 1D bar with 3 nodes (each having 1 DOF):

```
 Node:   0     1     2
 DOF:   [0]   [1]   [2]
```

For a 2D problem with 2 DOFs per node (uₓ, uᵧ):

```
 Node:     0        1        2
 DOFs:   [0, 1]   [2, 3]   [4, 5]
```

---

## Example Problem

A steel bar of length 2.0 m is fixed at the left end and loaded by 1000 N at the right end.

| Parameter | Symbol | Value |
|------------|---------|-------|
| Cross-section | \( A \) | 0.01 m² |
| Young’s modulus | \( E \) | 210 GPa |
| Length | \( L \) | 2.0 m |
| Load | \( F \) | 1000 N |

```python
import numpy as np
from myfem import first_fe_code, global_dof, DIRICHLET, NEUMANN

coords = np.array([[0.0], [1.0], [2.0]])

blocks = [{
    "element": {"properties": {"area": 0.01}},
    "material": "steel",
    "connect": [[0, 1], [1, 2]]
}]

bcs = [
    {"type": DIRICHLET, "nodes": [0], "local_dof": 0, "value": 0.0},  # fixed
    {"type": NEUMANN, "nodes": [2], "local_dof": 0, "value": 1000.0}  # traction
]

materials = {"steel": {"parameters": {"E": 210e9}, "density": 7800}}
dloads = []
block_elem_map = {0: (0, 0), 1: (0, 1)}

solution = first_fe_code(coords, blocks, bcs, dloads, materials, block_elem_map)
print("Displacements:", solution["dofs"])
```

Expected displacement at node 2:

\[
u_2 = \frac{F L}{A E} = \frac{1000 \times 2}{0.01 \times 210\times10^9} \approx 9.52 \times 10^{-6}\ \text{m}
\]

---

## Visualization

A schematic of the 1D element assembly:

```
     Node 0               Node 1               Node 2
   (x=0)●──── Element 0 ───●──── Element 1 ───●(x=2)
        │                  │                  │
      u₀=0                u₁=?               u₂=?
```

Global system after assembly:

\[
\begin{bmatrix}
K_{11} & K_{12} & 0 \\
K_{21} & K_{22}+K_{33} & K_{34} \\
0 & K_{43} & K_{44}
\end{bmatrix}
\begin{Bmatrix}
u_0 \\ u_1 \\ u_2
\end{Bmatrix}
=
\begin{Bmatrix}
0 \\ 0 \\ 1000
\end{Bmatrix}
\]

---

## Notes for Developers

- The solver assumes **linear elasticity** and **small deformations**.
- Only **1D, two-node, linear bar elements** are currently supported.
- The code is modular and can be extended to include:
  - Higher-order elements
  - Thermal or dynamic effects
  - Multi-dimensional formulations

---

## References

- T.J.R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*, Dover Publications.
- A. Fish & T. Belytschko, *A First Course in Finite Elements*, Wiley, 2007.
- T. Fullem, *fem-with-python* (GitHub): https://github.com/tjfulle/fem-with-python

---
