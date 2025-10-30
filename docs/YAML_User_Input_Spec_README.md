# YAML User Input Spec for the 1‑D FE Solver (`wundy`)

This README explains how to write a **YAML** input file that passes the validation schemas shown in the code. It defines every section, the allowed values and types, defaults, and gives worked examples.

The top-level structure is:

```yaml
wundy:
  nodes:                # required
  elements:             # required
  boundary conditions:  # required (list; may be empty)
  materials:            # required
  element blocks:       # required
  node sets:            # optional
  element sets:         # optional
  concentrated loads:   # optional
  distributed loads:    # optional
```

> **Scope:** This version targets a **1‑D truss/bar** formulation with element type `T1D1`. Degrees of freedom are translational in the global \(X\) direction only. Extensions to 2‑D/3‑D are noted where relevant.

---

## 1) Nodes

**Schema:** a list of lists. Each inner list is `[node_id, x]` for 1‑D.

- `node_id`: integer (label/index of the node)
- `x`: float (the node coordinate along the global \(X\) axis)

```yaml
nodes:
  - [1, 0.0]
  - [2, 1.0]
  - [3, 2.0]
```

**Validation rules:**

- First entry of each inner list must be an `int`.
- Remaining entries are numeric; they are coerced to `float`.
- In 1‑D, supply exactly one coordinate per node.

---

## 2) Elements (Connectivity)

**Schema:** a list of integer lists. Each inner list is the node connectivity of one element.

- For `T1D1`, each element connects **two** node IDs: `[n1, n2]`.

```yaml
elements:
  - [1, 2]
  - [2, 3]
```

**Validation rules:**

- Every inner list must contain integers only.
- Node IDs must refer to existing nodes.

---

## 3) Materials

**Schema:** a list of material objects.

- `type` (required): only `"ELASTIC"` is supported.
- `name` (required): case-insensitive label; normalized to uppercase.
- `parameters` (required): must include:
  - `E` \(> 0\): Young’s modulus.
  - `nu` \(\in [-1, 0.5)\): Poisson’s ratio (kept for forward compatibility; not used in 1‑D axial).
- `density` (optional, default = `0.0`, must be \(>0\) if present).

```yaml
materials:
  - type: elastic
    name: steel_a
    parameters:
      E: 210e9         # Pa
      nu: 0.3
    density: 7850      # kg/m^3  (optional)
```

**Validation rules:**

- `type` and `name` are normalized to uppercase.
- `E` must be positive; `nu` must satisfy \(-1 \le \nu < 0.5\).

---

## 4) Element Blocks

Element blocks bind:
1) a **material**,  
2) a set of **elements** (by set name or explicit list), and  
3) the **element definition** (type and properties).

**Schema fields:**

- `name` (required): block label (normalized to uppercase).
- `material` (required): must match a material `name`.
- `elements` (required): either an element-set name **or** a list of element IDs.
- `element` (required):
  - `type`: must be `T1D1`.
  - `properties`:
    - `area` (optional; default = `1.0`; must be numeric and \(>0\)).

```yaml
element blocks:
  - name: blk1
    material: steel_a
    elements: [1, 2]             # or an elset name like "ESET_LEFT"
    element:
      type: T1D1
      properties:
        area: 1.25e-4            # m^2
```

**Validation rules:**

- `type` is checked against the allowed set `{T1D1}`.
- `area` is validated and defaulted if missing.

---

## 5) Node Sets (Optional)

Named groups of nodes.

```yaml
node sets:
  - name: left
    nodes: [1]
  - name: right
    nodes: [3]
```

- `name`: normalized to uppercase.
- `nodes`: list of integers.

---

## 6) Element Sets (Optional)

Named groups of elements.

```yaml
element sets:
  - name: midspan
    elements: [2]
```

- `name`: normalized to uppercase.
- `elements`: list of integers.

---

## 7) Boundary Conditions (Essential/Dirichlet or Natural/Neumann)

Each entry describes constraints or prescribed tractions on **nodes**.

**Fields:**

- `nodes` (required): one of
  - a node-set name (string), or
  - a single node id (int), or
  - a list of node ids (list[int]).
- `dof` (optional; default = `0`): degree of freedom, specified as a **string** DOF id that is mapped to an enum:
  - Allowed strings: `"X"` (1‑D).  
    Internally, `"X" → 0`.
- `type` (optional; default = `DIRICHLET`): `"DIRICHLET"` or `"NEUMANN"` (case-insensitive).
- `value` (optional; default = `0.0`): numeric value.  
  - For \( \text{DIRICHLET} \): displacement \(u\) in meters.  
  - For \( \text{NEUMANN} \): nodal force \(F\) in newtons (if you choose to represent natural BCs this way).
- `name` (optional): label (normalized to uppercase).

```yaml
boundary conditions:
  # Fix the left end: u(1) = 0
  - nodes: left
    dof: X                # optional; defaults to X (enum 0)
    type: DIRICHLET
    value: 0.0

  # Apply a nodal force at node 3: F = +1000 N
  - nodes: 3
    dof: X
    type: NEUMANN
    value: 1000.0
```

**Validation rules:**

- `nodes` may be a set name, a single int (coerced to a list), or a list of ints.
- `dof` string must be valid (`X`), and is mapped to enum `0`.
- `type` is validated and mapped to internal enums: `DIRICHLET → 1`, `NEUMANN → 0`.

---

## 8) Concentrated Loads (Optional)

Alternative place to define **nodal** loads, separate from boundary conditions. Same addressing of `nodes` and `dof`.

**Fields:**

- `nodes`: node-set name, single node id, or list of node ids.
- `dof` (optional; default = `0`, i.e., `"X"`).
- `value` (optional; default = `0.0`): nodal force magnitude.
- `name` (optional).

```yaml
concentrated loads:
  - nodes: [2, 3]
    dof: X
    value: -500.0
```

---

## 9) Distributed Loads (Optional; Body/Element Loads)

Loads applied **per element** (e.g., body force). In 1‑D, the direction vector has **length 1**.

**Fields:**

- `elements` (required): element-set name, single element id, or list of element ids.
- `type` (required): one of `"BX"` or `"GRAV"` (case-insensitive).
- `value` (required): `float` magnitude (e.g., acceleration \(g\), or body-force density).
- `direction` (required): list of length 1 in 1‑D. Each entry is `float`.
  - Example: `[1.0]` to act in +\(X\); `[-1.0]` to act in \(-X\).
- `name` (optional).

```yaml
distributed loads:
  - elements: midspan
    type: grav
    value: 9.81             # m/s^2
    direction: [ -1.0 ]     # act in -X
```

**Validation rules:**

- `type` normalized and validated (`BX`, `GRAV`).
- `direction` length must be exactly 1 in 1‑D; entries coerced to float.

---

## 10) Full Minimal Working Example

```yaml
wundy:
  nodes:
    - [1, 0.0]
    - [2, 1.0]
    - [3, 2.0]

  elements:
    - [1, 2]
    - [2, 3]

  node sets:
    - name: left
      nodes: [1]
    - name: right
      nodes: [3]

  element sets:
    - name: midspan
      elements: [2]

  materials:
    - type: elastic
      name: steel_a
      parameters:
        E: 210e9
        nu: 0.3

  element blocks:
    - name: blk1
      material: steel_a
      elements: [1, 2]        # or "midspan", etc.
      element:
        type: T1D1
        properties:
          area: 1.0e-4

  boundary conditions:
    - nodes: left
      dof: X
      type: dirichlet
      value: 0.0
    - nodes: right
      dof: X
      type: neumann
      value: 1000.0

  # Optional alternative to put nodal loads here:
  # concentrated loads:
  #   - nodes: right
  #     dof: X
  #     value: 1000.0

  # Optional body/element loads:
  # distributed loads:
  #   - elements: midspan
  #     type: grav
  #     value: 9.81
  #     direction: [-1.0]
```

---

## 11) Degrees of Freedom (DOF) Mapping

- DOF identifiers are given as strings and normalized to uppercase.
- In 1‑D:
  - `"X" → 0`.

Formally, the mapping implemented is:
\[
\text{dof\_id\_to\_enum}(\text{"X"}) = 0.
\]

---

## 12) Case Handling and Normalization

- Many string fields are **case-insensitive** and are normalized to uppercase internally:
  - `type`, `name`, set names, `DIRICHLET`/`NEUMANN`, `BX`/`GRAV`, etc.
- This means `"elastic"`, `"ELASTIC"`, and `"Elastic"` are all accepted as `"ELASTIC"`.

---

## 13) Units and Physical Meaning

There is no internal unit system. Be consistent:

- Coordinates \(x\): meters (m).
- Area \(A\): \( \mathrm{m}^2 \).
- Young’s modulus \(E\): pascals (Pa \(=\mathrm{N/m^2}\)).
- Density \(\rho\): \( \mathrm{kg/m^3} \) (if used).
- Displacement \(u\): meters (m).
- Nodal forces \(F\): newtons (N).
- Body load magnitude (e.g., gravity): \( \mathrm{m/s^2} \) if interpreted as acceleration; or \( \mathrm{N/kg} \) if you encode it that way—keep it consistent with how your solver converts `value × density × area` into element forces.

---

## 14) Common Validation Errors and How to Fix Them

1. **Unknown element type**  
   - Only `T1D1` is allowed:
     ```yaml
     element:
       type: T1D1
     ```

2. **Missing or invalid material parameters**  
   - Provide both `E` \(>0\) and `nu` \(\in [-1, 0.5)\):
     ```yaml
     parameters:
       E: 210e9
       nu: 0.3
     ```

3. **`area` not positive or missing**  
   - `area` defaults to `1.0` if omitted. If provided, it must be \(>0\).

4. **`direction` length incorrect in `distributed loads`**  
   - In 1‑D the list must have **exactly one** entry.

5. **Mixed or wrong case names**  
   - Internally normalized, but for readability prefer consistent casing.

6. **`nodes` or `elements` referencing non-existent ids or sets**  
   - Ensure all referenced ids and set names are defined.

---

## 15) Extending to 2‑D / 3‑D (Roadmap Notes)

The code contains hooks for higher dimensions:

- `valid_dof_id` would admit `"X"`, `"Y"` (and `"Z"` in 3‑D).
- `direction` in `distributed loads` would have length 2 (2‑D) or 3 (3‑D).
- `node_freedom_table` would return tuples indicating active DOFs per node for each element type.
- Additional element types (e.g., `T2D2`, `T3D2`, etc.) would be added to `element_types` and validated in `validate_element`.

---

## 16) Reference: Allowed Enums and Mappings

- **Boundary condition types**
  - User strings: `DIRICHLET`, `NEUMANN` (case-insensitive).
  - Internal enums:
    \[
    \text{bc\_type\_to\_enum}(\text{"DIRICHLET"}) = 1,\quad
    \text{bc\_type\_to\_enum}(\text{"NEUMANN"}) = 0.
    \]

- **Distributed load types**
  - User strings: `BX`, `GRAV` (case-insensitive).

- **Element types**
  - `T1D1` only (in this version).

---

## 17) Validation Checklist (Pre‑Run)

- [ ] All required top-level sections under `wundy` are present.  
- [ ] Nodes: each is `[int, float]` for 1‑D.  
- [ ] Elements: each is a list of integers; connectivity is valid.  
- [ ] Materials: at least one `ELASTIC` with `E>0`, `nu∈[-1,0.5)`.  
- [ ] Element blocks: `type=T1D1`; `area>0` (or omitted); `material` name matches.  
- [ ] Boundary conditions: valid `nodes` selector; `dof="X"`; `type` in {DIRICHLET, NEUMANN}; numeric `value`.  
- [ ] Optional loads: concentrated loads use node selectors; distributed loads use element selectors, valid `type`, numeric `value`, and `direction` of length 1.  
- [ ] All referenced set names exist.

---

## 18) Starter Template

Copy, edit, and save as `input.yaml`:

```yaml
wundy:
  nodes: []
  elements: []

  materials:
    - type: ELASTIC
      name: MAT1
      parameters: { E: 1.0e7, nu: 0.3 }

  element blocks:
    - name: BLK1
      material: MAT1
      elements: []
      element:
        type: T1D1
        properties: { area: 1.0 }

  node sets: []
  element sets: []

  boundary conditions: []
  concentrated loads: []
  distributed loads: []
```

---

### Notes on Math and Sign Conventions

- Dirichlet BCs prescribe displacement \( u \) (e.g., \( u = 0 \) at a fixed end).
- Neumann BCs prescribe force \( F \) at nodes (e.g., \( F = +1000\,\mathrm{N} \) in \(+X\)).
- Distributed loads interpreted as body/element loads typically induce equivalent nodal forces in assembly; ensure the combination of `value`, `direction`, `density`, and `area` is consistent so the resulting forces have units of newtons.

---

That’s the complete specification for preparing a valid YAML input for the 1‑D solver as written. To expand capability (2‑D/3‑D elements, additional DOFs, more load types), extend the enumerations and validation helpers in parallel with your solver implementation.
