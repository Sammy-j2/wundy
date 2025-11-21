import logging
import re
from typing import Any, IO

import numpy as np
import yaml

from .schemas import NEUMANN, input_schema

logger = logging.getLogger(__name__)


"""YAML loader and preprocessor for wundy.

This module exposes two user-facing helpers:

- `load(file)`: read YAML from a file-like object, coerce some common
    numeric-looking string tokens into numbers, and validate against the
    input schema defined in `src/wundy/schemas.py`.
- `preprocess(data)`: transform the validated YAML into the compact
    Python structures used by the assembler/solver (`coords`, `blocks`,
    `bcs`, `dload`, `materials`, `block_elem_map`, etc.).

The preprocessor also handles these important behaviors:

- `dof_per_node`: reads the user-requested number of DOFs per node and
    stores it in the preprocessed output. This enables the solver's
    uniform multi-DOF mode where the axial element is expanded into a
    `2 * dof_per_node` sized element contribution (axial DOF is local
    index 0).
- Integration options: global integration defaults are merged with any
    user-provided global `integration` block and with per-element-block
    `integration` overrides. The resulting dict is stored on each block as
    `block['integration']` for use during assembly.
 - Numeric coercion: the loader is forgiving of numeric-looking strings
     in material parameters and `integration` settings (for example `"1e-8"`).
     These are coerced to `int`/`float` before schema validation so YAML
     produced by different editors still validates.
"""


def load(file: IO[Any]) -> dict[str, dict[str, Any]]:
    """Load YAML from a file-like object, coerce numeric-like strings, then validate.

    Some YAML parsers or editor workflows may leave numeric-looking tokens as
    strings (e.g. "2.1e11"). To be forgiving, coerce any material parameters
    that look like floats into actual floats before schema validation.
    """
    text = file.read()
    data = yaml.safe_load(text)

    # Coerce numeric-looking strings in material parameters to floats so the
    # schema's numeric checks succeed even if PyYAML returned strings.
    float_re = re.compile(r"^[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?$")
    try:
        mats = data.get("wundy", {}).get("materials", [])
        for m in mats:
            params = m.get("parameters", {})
            for k, v in list(params.items()):
                if isinstance(v, str) and float_re.match(v.strip()):
                    try:
                        params[k] = float(v)
                    except Exception:
                        # leave as-is if conversion fails
                        pass
    except Exception:
        # be conservative: if data shape unexpected, skip coercion and let
        # schema validation handle errors
        pass

    # Coerce numeric-looking strings in 'integration' dictionaries so values
    # like "1e-8" or "2" are accepted by the schema numeric validators.
    try:
        wundy = data.get("wundy", {})
        integ = wundy.get("integration")
        if isinstance(integ, dict):
            for k, v in list(integ.items()):
                if isinstance(v, str) and float_re.match(v.strip()):
                    sval = v.strip()
                    if re.match(r"^[+-]?\d+$", sval):
                        try:
                            integ[k] = int(sval)
                            continue
                        except Exception:
                            pass
                    try:
                        integ[k] = float(sval)
                    except Exception:
                        pass

        # Per-element-block integration entries (if present)
        for eb in wundy.get("element blocks", []) or []:
            belem = eb.get("element", {})
            binteg = belem.get("integration")
            if isinstance(binteg, dict):
                for k, v in list(binteg.items()):
                    if isinstance(v, str) and float_re.match(v.strip()):
                        sval = v.strip()
                        if re.match(r"^[+-]?\d+$", sval):
                            try:
                                binteg[k] = int(sval)
                                continue
                            except Exception:
                                pass
                        try:
                            binteg[k] = float(sval)
                        except Exception:
                            pass
    except Exception:
        # If anything unexpected happens here, leave data as-is and allow
        # the schema validation to catch/describe the problem.
        pass

    return input_schema.validate(data)


def set_element_defaults(elem: dict[str, Any]) -> bool:
    if elem["type"].upper() == "T1D1":
        nft = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        props = {"node_per_elem": 2, "freedom_table": [nft, nft]}
        elem["properties"].update(props)
    else:
        raise ValueError(f"Unknown element type {elem['type']!r}")
    return True


def unique_name(named_items: list[dict], stem: str) -> str:
    names = [item.get("name") for item in named_items]
    i = 1
    while True:
        name = f"{stem.upper()}-{i}"
        if name not in names:
            return name
        i += 1


def preprocess(data: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Preprocess and transform user input.

    Assumptions: User input was loaded and validated by ``load``

    """
    errors: int = 0

    inp = data["wundy"]

    preprocessed: dict[str, Any] = {}

    num_node: int = len(inp["nodes"])
    max_dim: int = max(len(n[1:]) for n in inp["nodes"])
    node_map: dict[int, int] = preprocessed.setdefault("node_map", {})
    coords = preprocessed["coords"] = np.zeros((num_node, max_dim))
    for i, node in enumerate(inp["nodes"]):
        nid, *xc = node
        node_map[nid] = i
        coords[i, : len(xc)] = xc

    # Degrees of freedom per node: allow user to request 1/2/3 DOFs per node
    try:
        dof_per_node = int(inp.get("dof_per_node", 1))
    except Exception:
        dof_per_node = 1
    preprocessed["dof_per_node"] = dof_per_node
    # Note: `dof_per_node` influences validation of BC/load local DOF
    # indices later in this function. The solver assumes the axial DOF is
    # the local index 0 for each node when expanding element contributions
    # to the multi-DOF element vector/matrix.

    num_elem: int = len(inp["elements"])
    elem_map: dict[int, int] = preprocessed.setdefault("elem_map", {})
    # Elements expected format: [eid, n1, n2, ...]. Validate shape early so
    # downstream code (assembly/first_fe_code) does not hit obscure IndexErrors.
    for i, element in enumerate(inp["elements"]):
        if not isinstance(element, list) or len(element) < 3:
            raise ValueError(
                "Each element entry must be a list: [elem_id, node1, node2, ...] with at least two nodes"
            )
        elem_map[element[0]] = i

    # Put node sets in dictionary for easier look up
    nsets: dict[str, Any] = preprocessed.setdefault("nsets", {})
    nsets["all"] = list(range(num_node))
    for ns in inp.get("node sets", []):
        name = ns["name"]
        if name in nsets:
            errors += 1
            logger.error(f"Duplicate node set {name!r}")
        else:
            nodes: list[int] = []
            for n in ns["nodes"]:
                if n not in node_map:
                    errors += 1
                    logger.error(f"Node {n} in node set {name} is not defined")
                else:
                    nodes.append(node_map[n])
            nsets[name] = nodes

    # Put element sets in dictionary for easier look up
    elsets: dict[str, Any] = preprocessed.setdefault("elsets", {})
    elsets["ALL"] = list(range(num_elem))
    for es in inp.get("element sets", []):
        name = es["name"]
        if name in elsets:
            errors += 1
            logger.error(f"Duplicate element set {name!r}")
        else:
            elems: list[int] = []
            for e in es["elements"]:
                if e not in elem_map:
                    errors += 1
                    logger.error(f"Element {e} in element set {name} is not defined")
                else:
                    elems.append(elem_map[e])
            elsets[name] = elems

    # Put materials in dictionary for easier look up
    materials: dict[str, Any] = preprocessed.setdefault("materials", {})
    for material in inp["materials"]:
        name = material["name"]
        if name in materials:
            errors += 1
            logger.error(f"Duplicate material {name!r}")
        else:
            materials[name] = {"type": material["type"], "parameters": material["parameters"]}
    
        # Integration options: global defaults merged with user-provided top-level settings
        default_integration = {
            "stiffness": "analytic",
            "internal": "analytic",
            "ngp": 2,
            "nonlinear": "linearize",
        }
        user_integration = inp.get("integration") or {}
        # ensure ngp is an int when supplied as string
        try:
            if isinstance(user_integration.get("ngp", None), str):
                user_integration["ngp"] = int(user_integration["ngp"])
        except Exception:
            pass
        global_integration = {**default_integration, **user_integration}
        preprocessed["integration"] = global_integration

    # Put element blocks in dictionary for easier look up
    blocks: list[Any] = preprocessed.setdefault("blocks", [])
    for eb in inp["element blocks"]:
        name = eb["name"]
        if name in blocks:
            errors += 1
            logger.error(f"Duplicate element block {name!r}")
            continue
        if eb["material"] not in materials:
            errors += 1
            defined = sorted(materials.keys())
            logger.error(
                "material %r required by element block %r not defined; defined materials: %s",
                eb["material"],
                name,
                defined,
            )
            continue
        block: dict[str, Any] = {}
        block["name"] = name
        block["element"] = eb["element"]
        set_element_defaults(block["element"])
        block["material"] = eb["material"]
        elems: list[int] = []
        if isinstance(eb["elements"], str):
            # elements given as set name
            if eb["elements"] not in elsets:
                errors += 1
                logger.error(
                    f"element set {eb['elements']!r}, required by element block {name}, not defined"
                )
                continue
            elems.extend(elsets[eb["elements"]])
        else:
            for e in eb["elements"]:
                if e not in elem_map:
                    errors += 1
                    logger.error(f"Element {e}, required for element block {name}, is not defined")
                else:
                    elems.append(elem_map[e])
        if not elems:
            errors += 1
            logger.error(f"No elements defined for element block {name}")
        else:
            connect: list[list[int]] = []
            for e in elems:
                eid, *nodes = inp["elements"][e]
                if connect and len(nodes) != len(connect[0]):
                    errors += 1
                    logger.error(
                        f"Inconsistent element connectivity in element block {name}. "
                        "(All elements must have the same number of nodes)"
                    )
                    break
                row: list[int] = []
                for n in nodes:
                    if n not in node_map:
                        errors += 1
                        logger.error(
                            f"Node {n} of element {eid} in element block {name} is not in node map"
                        )
                    else:
                        row.append(node_map[n])
                connect.append(row)
            else:
                block["connect"] = np.array(connect, dtype=int)
                # Merge per-block integration options with global defaults
                block_user_integ = eb.get("integration") or {}
                try:
                    if isinstance(block_user_integ.get("ngp", None), str):
                        block_user_integ["ngp"] = int(block_user_integ["ngp"])
                except Exception:
                    pass
                block["integration"] = {**preprocessed.get("integration", {}), **block_user_integ}
                # Map from global index to local index
                block["elem_map"] = dict(zip(elems, range(len(elems))))
                blocks.append(block)

    # Convert boundary conditions to tags/vals that can be used by the assembler
    boundary: list[Any] = preprocessed.setdefault("bcs", [])
    for i, bc in enumerate(inp["boundary conditions"]):
        if "name" in bc:
            name = bc["name"]
        else:
            name = unique_name(inp["boundary conditions"], "BOUNDARY")
            bc["name"] = name
        nodes: list[int] = []
        if isinstance(bc["nodes"], str):
            if bc["nodes"] not in nsets:
                errors += 1
                logger.error(
                    f"Nodeset {bc['nodes']}, required by boundary condition {i + 1}, is not defined"
                )
            else:
                nodes.extend(nsets[bc["nodes"]])
        else:
            for n in bc["nodes"]:
                if n not in node_map:
                    errors += 1
                    logger.error(
                        f"Node {n}, required by boundary condition {i + 1}, is not defined"
                    )
                else:
                    nodes.append(node_map[n])
        # Validate requested local_dof against configured dof_per_node. The
        # input schema converts symbolic DOF names like 'X'/'Y' to numeric
        # local indices; here we ensure the index is within the declared
        # range [0, dof_per_node).
        local_dof = bc["dof"]
        if not (0 <= int(local_dof) < preprocessed["dof_per_node"]):
            errors += 1
            logger.error(
                f"Boundary condition {name} requests local DOF {local_dof}, but dof_per_node={preprocessed['dof_per_node']}"
            )
        boundary.append(
            {
                "name": name,
                "local_dof": int(local_dof),
                "type": bc["type"],
                "nodes": nodes,
                "value": bc["value"],
            }
        )

    # Convert concentrated loads to tags/vals that can be used by the assembler
    for i, cl in enumerate(inp.get("concentrated loads", [])):
        if "name" in cl:
            name = cl["name"]
        else:
            name = unique_name(inp["concentrated loads"], "CLOAD")
            cl["name"] = name
        nodes: list[int] = []
        if isinstance(cl["nodes"], str):
            if cl["nodes"] not in nsets:
                errors += 1
                logger.error(
                    f"Nodeset {cl['nodes']}, required by concentrated load {i + 1}, is not defined"
                )
            else:
                nodes.extend(nsets[cl["nodes"]])
        else:
            for n in cl["nodes"]:
                if n not in node_map:
                    errors += 1
                    logger.error(f"Node {n}, required by concentrated load {i + 1}, is not defined")
                else:
                    nodes.append(node_map[n])
        local_dof = cl["dof"]
        if not (0 <= int(local_dof) < preprocessed["dof_per_node"]):
            errors += 1
            logger.error(
                f"Concentrated load {name} requests local DOF {local_dof}, but dof_per_node={preprocessed['dof_per_node']}"
            )
        boundary.append(
            {
                "name": name,
                "local_dof": int(local_dof),
                "type": NEUMANN,
                "nodes": nodes,
                "value": cl["value"],
            }
        )

    # Process distributed load
    dload: list[Any] = preprocessed.setdefault("dload", [])
    for i, dl in enumerate(inp.get("distributed loads", [])):
        if "name" in dl:
            name = dl["name"]
        else:
            name = unique_name(inp["distributed loads"], "DLOAD")
            dl["name"] = name
        elems: list[int] = []
        if isinstance(dl["elements"], str):
            # elements given as set name
            if dl["elements"] not in elsets:
                errors += 1
                logger.error(
                    f"Element set {dl['elements']!r}, required by distributed load {i + 1}, not defined"
                )
            else:
                elems.extend(elsets[eb["elements"]])
        else:
            for e in dl["elements"]:
                if e not in elem_map:
                    errors += 1
                    logger.error(
                        f"Element {e}, required by distributed load {i + 1}, is not defined"
                    )
                else:
                    elems.append(elem_map[e])
        dload.append(
            {
                "name": name,
                "elements": elems,
                "type": dl["type"],
                "value": dl["value"],
                "direction": dl["direction"],
            }
        )

    # Create a mapping from global element index to block index, local elem index (within the block)
    block_elem_map: dict[int, tuple[int, int]] = preprocessed.setdefault("block_elem_map", {})
    for ib, block in enumerate(blocks):
        for global_elem_index, local_elem_index in block["elem_map"].items():
            if global_elem_index in block_elem_map:
                errors += 1
                logger.error(f"Duplicate element ID {e} found in multiple blocks")
            block_elem_map[global_elem_index] = (ib, local_elem_index)

    # Check if all elements are assigned to an element block
    if unassigned := set(range(num_elem)).difference(block_elem_map.keys()):
        errors += 1
        # Provide original element IDs to help the user correlate YAML ids
        unassigned_ids = [inp["elements"][e][0] for e in sorted(unassigned)]
        logger.error(
            "The following element indices (0-based) are unassigned: %s; corresponding YAML element IDs: %s",
            sorted(unassigned),
            unassigned_ids,
        )
        for e in sorted(unassigned):
            logger.error(f"Element {e} is not assigned to any element blocks")

    if errors:
        raise ValueError("Stopping due to previous errors")

    return preprocessed
