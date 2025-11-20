#!/usr/bin/env python3
"""Run a wundy YAML input file and print the solution DOFs.

Usage:
    python bin/run_yaml.py path/to/input.yaml

If no path is given the script uses:
    docs/examples/material_input_example.yaml
"""

import argparse
import sys
try:
    import tkinter as tk
    from tkinter import filedialog
except Exception:
    tk = None
    filedialog = None
from pathlib import Path

from schema import SchemaError
import yaml

from wundy import first, ui


def run_yaml(path: str, debug: bool = False) -> None:
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"YAML file not found: {p}")
    try:
        with p.open() as fh:
            data = ui.load(fh)
    except SchemaError as exc:
        # In non-debug mode print a short message and re-raise. In debug
        # mode print the full diagnostics we used while developing.
        print("Schema validation failed:", exc)
        print("YAML file absolute path:", p.resolve())
        if debug:
            # Re-open raw YAML and show detailed diagnostics
            with p.open("rb") as fh2:
                raw_bytes = fh2.read()
            try:
                raw_text = raw_bytes.decode("utf-8")
            except Exception:
                raw_text = raw_bytes.decode("latin-1", errors="replace")
            try:
                parsed = yaml.safe_load(raw_text)
            except Exception as e:
                parsed = f"(yaml.safe_load failed: {e})"

            print("YAML raw text (repr, first 800 chars):")
            print(repr(raw_text[:800]))
            print("YAML raw bytes (first 512 bytes repr):")
            print(repr(raw_bytes[:512]))
            print("PyYAML parsed object (materials entry):")
            mats = None
            if isinstance(parsed, dict):
                mats = parsed.get("wundy", {}).get("materials", [])
            print(repr(mats))
            try:
                root_node = yaml.compose(raw_text)

                def find_node(node, path):
                    if not path:
                        return node
                    head, *tail = path
                    if isinstance(node, yaml.nodes.MappingNode):
                        for k_node, v_node in node.value:
                            if isinstance(k_node, yaml.nodes.ScalarNode) and k_node.value == head:
                                return find_node(v_node, tail)
                        return None
                    if isinstance(node, yaml.nodes.SequenceNode):
                        try:
                            idx = int(head)
                        except Exception:
                            return None
                        if 0 <= idx < len(node.value):
                            return find_node(node.value[idx], tail)
                        return None
                    return None

                target_path = ["wundy", "materials", "0", "parameters", "E"]
                node = find_node(root_node, target_path)
                if node is None:
                    print("Could not locate YAML node for wundy->materials[0]->parameters->E")
                else:
                    print("YAML node for E:")
                    print("  node.tag =", getattr(node, "tag", None))
                    print("  node.value =", getattr(node, "value", None))
                    print("  node.style =", getattr(node, "style", None))
            except Exception as e:
                print("YAML node inspection failed:", e)
        # Re-raise so calling process sees the SchemaError
        raise
    inp = ui.preprocess(data)
    # Pass through the user's requested `dof_per_node` (if present in preprocessed input)
    dof_per_node = int(inp.get("dof_per_node", 1))
    sol = first.first_fe_code(
        inp["coords"],
        inp["blocks"],
        inp["bcs"],
        inp["dload"],
        inp["materials"],
        inp["block_elem_map"],
        dof_per_node=dof_per_node,
    )
    print("dofs:", sol["dofs"])


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run a wundy YAML input file")
    parser.add_argument(
        "yaml",
        nargs="?",
        default=None,
        help="Path to YAML input file. If omitted a file dialog will open",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print detailed YAML diagnostics on schema/parse errors",
    )
    parser.add_argument(
        "--pick",
        action="store_true",
        help="Open a file dialog to pick a YAML input file (Windows GUI)",
    )
    args = parser.parse_args(argv)
    # If the user didn't pass a YAML path or requested a GUI pick, show a file dialog
    yaml_path = args.yaml
    if args.pick or yaml_path is None:
        if filedialog is None:
            print("tkinter not available; cannot open file picker")
            sys.exit(1)
        # create a hidden root and open file dialog
        root = tk.Tk()
        root.withdraw()
        p = filedialog.askopenfilename(
            title="Select a wundy YAML input file",
            filetypes=[("YAML files", "*.yaml;*.yml"), ("All files", "*.*")],
        )
        root.destroy()
        if not p:
            print("No file selected; exiting.")
            sys.exit(1)
        yaml_path = p

    run_yaml(yaml_path, debug=args.debug)


if __name__ == "__main__":
    main()
