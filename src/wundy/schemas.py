from typing import Any

from schema import And
from schema import Optional
from schema import Or
from schema import Schema
from schema import Use

NEUMANN = 0
DIRICHLET = 1


element_types = {"T1D1"}
bc_types = {"DIRICHLET", "NEUMANN"}


def node_freedom_table(elem_type: str) -> tuple[int, ...]:
    if normalize_case(elem_type) == "T1D1":
        return (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    raise ValueError(f"Unknown element type {elem_type!r}")


def valid_element_type(name: str) -> bool:
    return normalize_case(name) in element_types


def isnumeric(x) -> bool:
    return isinstance(x, (int, float))


def ispositive(arg: float | int) -> bool:
    return arg > 0


def list_of_type(sequence: list, type) -> bool:
    return all(isinstance(n, type) for n in sequence)


def list_of_numeric(sequence: list) -> bool:
    return all(isinstance(x, (float, int)) for x in sequence)


def list_of_int(sequence) -> bool:
    return list_of_type(sequence, int)


def list_of_list(sequence) -> bool:
    return list_of_type(sequence, list)


def normalize_case(string: str) -> str:
    return string.upper()


def dof_id_to_enum(dof: str) -> int:
    return {"X": 0, "Y": 1, "Z": 2}[normalize_case(dof)]


def valid_bc_type(arg: str) -> bool:
    return normalize_case(arg) in bc_types


def bc_type_to_enum(bc_type: str) -> int:
    return {"DIRICHLET": DIRICHLET, "NEUMANN": NEUMANN}[normalize_case(bc_type)]


def valid_dof_id(dof: str):
    # extension to 2/3D: allow dof to be xyz
    return normalize_case(dof) in {"X"}


def valid_dload_type(arg: str):
    # extension to 2/3D: allow other DLOADs
    return normalize_case(arg) in {"BX", "GRAV"}


def validate_element(elem: dict[str, Any]) -> bool:
    if normalize_case(elem["type"]) == "T1D1":
        schema = Schema({Optional("area", default=1.0): And(isnumeric, ispositive)})
        v = schema.validate(elem["properties"])
        elem["properties"].update(v)
    else:
        raise ValueError(f"Unknown element type {elem['type']!r}")
    return True


def validate_material_parameters(material: dict[str, dict[str, Any]]) -> bool:
    elastic = Schema(
        {
            "E": And(isnumeric, ispositive, error="E must be > 0"),
            "nu": And(isnumeric, lambda x: -1.0 <= x < 0.5, error="nu must be between -1 and .5"),
        }
    )
    # Support for Neo-Hookean materials: accept either the same E/nu pair
    # or Lame/mu style parameters (mu and lambda). The solver currently
    # linearizes Neo-Hookean materials using an equivalent Young's modulus
    # for assembly, so we allow both forms for user convenience.
    mt = normalize_case(material["type"])
    if mt == "ELASTIC":
        elastic.validate(material["parameters"])
    elif "NEO" in mt and "HOOK" in mt:
        # Neo-Hookean can be specified either as (E, nu) or (mu, lambda)
        neo_lame = Schema({"mu": And(isnumeric, ispositive), "lambda": And(isnumeric)})
        try:
            # Prefer E/nu if provided
            elastic.validate(material["parameters"])
        except Exception:
            # Fall back to mu/lambda
            neo_lame.validate(material["parameters"])
    else:
        raise ValueError(f"Unknown material {material['type']!r}")
    return True


nodes_schema = Schema(
    And(
        list,
        list_of_list,
        lambda outer: all(isinstance(inner[0], int) for inner in outer),  # node label
        lambda outer: all(isinstance(f, (int, float)) for inner in outer for f in inner[1:]),
        Use(lambda outer: [[int(inner[0]), *[float(_) for _ in inner[1:]]] for inner in outer]),
    )
)

elements_schema = Schema(
    And(
        list,
        list_of_list,
        lambda outer: all(list_of_int(inner) for inner in outer),
    )
)

nset_schema = Schema(
    {
        "name": And(str, Use(normalize_case)),
        "nodes": And(list, list_of_int),
    },
)

elset_schema = Schema(
    {
        "name": And(str, Use(normalize_case)),
        "elements": And(list, list_of_int),
    },
)

boundary_schema = Schema(
    And(
        {
            "nodes": Or(
                And(str, Use(normalize_case)),  # node set name
                And(int, Use(lambda n: [n])),  # single node
                And(list, list_of_int),  # list of nodes
            ),
            Optional("dof", default=0): And(str, valid_dof_id, Use(dof_id_to_enum)),
            Optional("name"): And(str, Use(normalize_case)),
            Optional("value", default=0.0): And(isnumeric, Use(float)),
            Optional("type", default=DIRICHLET): And(str, valid_bc_type, Use(bc_type_to_enum)),
        },
    )
)

cload_schema = Schema(
    And(
        {
            "nodes": Or(
                And(str, Use(normalize_case)),  # node set name
                And(int, Use(lambda n: [n])),  # single node
                And(list, list_of_int),  # list of nodes
            ),
            Optional("dof", default=0): And(str, valid_dof_id, Use(dof_id_to_enum)),
            Optional("name"): And(str, Use(normalize_case)),
            Optional("value", default=0.0): Use(float),
        },
    )
)

dload_schema = Schema(
    And(
        {
            "elements": Or(
                And(str, Use(normalize_case)),  # element set name
                And(int, Use(lambda e: [e])),  # single element
                And(list, list_of_int),  # list of elements
            ),
            "type": And(str, valid_dload_type, Use(normalize_case)),
            "value": Use(float),
            "direction": And(
                list,
                list_of_numeric,
                lambda sequence: len(sequence) == 1,  # change to <= 2/3 for 2D/3D
                Use(lambda sequence: [float(x) for x in sequence]),
            ),
            Optional("name"): And(str, Use(normalize_case)),
        },
    )
)

material_schema = Schema(
    And(
        {
            "type": And(str, Use(normalize_case)),
            "name": And(str, Use(normalize_case)),
            "parameters": {str: object},
            Optional("density", default=0.0): And(isnumeric, ispositive),
        },
        lambda d: validate_material_parameters(d),
    )
)

block_schema = Schema(
    And(
        {
            "name": And(str, Use(normalize_case)),
            "material": And(str, Use(normalize_case)),
            "elements": Or(
                And(str, Use(normalize_case)),
                And(list, list_of_int),
            ),
            "element": {
                "type": And(str, valid_element_type, Use(normalize_case)),
                Optional("properties", default=dict()): {str: object},
            },
        },
        lambda d: validate_element(d["element"]),
    )
)

input_schema = Schema(
    {
        "wundy": {
            "nodes": nodes_schema,
            "elements": elements_schema,
            "boundary conditions": [boundary_schema],
            "materials": [material_schema],
            "element blocks": [block_schema],
            Optional("node sets"): [nset_schema],
            Optional("element sets"): [elset_schema],
            Optional("concentrated loads"): [cload_schema],
            Optional("distributed loads"): [dload_schema],
        }
    }
)
