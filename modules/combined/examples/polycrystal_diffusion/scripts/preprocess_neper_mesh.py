#!/usr/bin/env python3
"""Collapse Neper's per-entity Gmsh groups into MOOSE-ready physical groups."""

import argparse
import json
from pathlib import Path


POINT = 15
EDGE2 = 1
TRI3 = 2


def read_section(lines, name):
    start = lines.index(f"${name}\n") + 1
    end = lines.index(f"$End{name}\n", start)
    return lines[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--width", type=float, default=200.0, help="x extent in mesh units")
    parser.add_argument("--height", type=float, default=50.0, help="y extent in mesh units")
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--stats", type=Path)
    args = parser.parse_args()

    lines = args.input.read_text(encoding="ascii").splitlines(keepends=True)
    node_section = read_section(lines, "Nodes")
    element_section = read_section(lines, "Elements")

    node_count = int(node_section[0])
    nodes = {}
    for line in node_section[1:]:
        fields = line.split()
        nodes[int(fields[0])] = tuple(float(value) for value in fields[1:4])
    if len(nodes) != node_count:
        raise ValueError("Neper node count does not match the $Nodes section")

    grouped_elements = []
    counts = {"bulk_tri3": 0, "internal_gb_edge2": 0, "left": 0, "right": 0,
              "bottom": 0, "top": 0}
    grain_ids = set()
    bulk_node_ids = set()
    gb_node_ids = set()
    tol = args.tolerance

    for line in element_section[1:]:
        fields = [int(value) for value in line.split()]
        element_id, element_type, number_of_tags = fields[:3]
        tags = fields[3:3 + number_of_tags]
        node_ids = fields[3 + number_of_tags:]
        geometrical_id = tags[0]

        if element_type == POINT:
            continue
        if element_type == TRI3:
            if len(node_ids) != 3:
                raise ValueError(f"Element {element_id} is not TRI3")
            grain_id = tags[1] if tags[0] == 1000 else tags[0]
            geometrical_id = grain_id
            grain_ids.add(grain_id)
            physical_id = 1000
            counts["bulk_tri3"] += 1
            bulk_node_ids.update(node_ids)
        elif element_type == EDGE2:
            if len(node_ids) != 2:
                raise ValueError(f"Element {element_id} is not EDGE2")
            coordinates = [nodes[node_id] for node_id in node_ids]
            if all(abs(point[0]) <= tol for point in coordinates):
                physical_id, key = 1, "left"
            elif all(abs(point[0] - args.width) <= tol for point in coordinates):
                physical_id, key = 2, "right"
            elif all(abs(point[1]) <= tol for point in coordinates):
                physical_id, key = 3, "bottom"
            elif all(abs(point[1] - args.height) <= tol for point in coordinates):
                physical_id, key = 4, "top"
            else:
                physical_id, key = 10, "internal_gb_edge2"
                gb_node_ids.update(node_ids)
            counts[key] += 1
        else:
            raise ValueError(f"Unsupported Neper element type {element_type}")

        grouped_elements.append(
            f"{element_id} {element_type} 2 {physical_id} {geometrical_id} "
            + " ".join(str(node_id) for node_id in node_ids)
        )

    if len(grain_ids) != 100:
        raise ValueError(f"Expected 100 grain IDs, found {len(grain_ids)}")

    physical_names = [
        (1, 1, "left"),
        (1, 2, "right"),
        (1, 3, "bottom"),
        (1, 4, "top"),
        (1, 10, "grain_boundaries"),
        (2, 1000, "bulk"),
    ]
    output_lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$PhysicalNames",
                    str(len(physical_names))]
    output_lines.extend(f'{dim} {identifier} "{name}"' for dim, identifier, name in physical_names)
    output_lines.extend(["$EndPhysicalNames", "$Nodes", str(node_count)])
    output_lines.extend(line.rstrip("\n") for line in node_section[1:])
    output_lines.extend(["$EndNodes", "$Elements", str(len(grouped_elements))])
    output_lines.extend(grouped_elements)
    output_lines.append("$EndElements")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(output_lines) + "\n", encoding="ascii")

    stats = {"grains": len(grain_ids), "nodes": node_count, **counts}
    stats["msh_elements"] = len(grouped_elements)
    stats["gb_nodes"] = len(gb_node_ids)
    stats["gb_nodes_shared_with_bulk"] = gb_node_ids <= bulk_node_ids
    if args.stats:
        args.stats.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps(stats, sort_keys=True))


if __name__ == "__main__":
    main()
