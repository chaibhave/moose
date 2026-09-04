#!/usr/bin/env python3
"""Collapse a 3D Neper Gmsh mesh into bulk, exterior, and internal-GB groups."""

import argparse
import json
from pathlib import Path


POINT = 15
EDGE2 = 1
TRI3 = 2
TET4 = 4


def read_section(lines, name):
    start = lines.index(f"${name}\n") + 1
    end = lines.index(f"$End{name}\n", start)
    return lines[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--length", type=float, default=200.0)
    parser.add_argument("--width", type=float, default=50.0)
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

    counts = {"bulk_tet4": 0, "internal_gb_tri3": 0, "left": 0, "right": 0,
              "y_min": 0, "y_max": 0, "z_min": 0, "z_max": 0}
    grain_ids = set()
    bulk_node_ids = set()
    gb_node_ids = set()
    grouped_elements = []
    tol = args.tolerance

    for line in element_section[1:]:
        fields = [int(value) for value in line.split()]
        element_id, element_type, number_of_tags = fields[:3]
        tags = fields[3:3 + number_of_tags]
        node_ids = fields[3 + number_of_tags:]
        geometrical_id = tags[0]

        if element_type in (POINT, EDGE2):
            continue
        if element_type == TET4:
            grain_id = tags[1] if tags[0] == 1000 else tags[0]
            grain_ids.add(grain_id)
            physical_id = 1000
            geometrical_id = grain_id
            counts["bulk_tet4"] += 1
            bulk_node_ids.update(node_ids)
        elif element_type == TRI3:
            coordinates = [nodes[node_id] for node_id in node_ids]
            planes = (
                (0, 0.0, 1, "left"),
                (0, args.length, 2, "right"),
                (1, 0.0, 3, "y_min"),
                (1, args.width, 4, "y_max"),
                (2, 0.0, 5, "z_min"),
                (2, args.width, 6, "z_max"),
            )
            match = next(((identifier, key) for axis, value, identifier, key in planes
                          if all(abs(point[axis] - value) <= tol for point in coordinates)), None)
            if match:
                physical_id, key = match
            else:
                physical_id, key = 10, "internal_gb_tri3"
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
        (2, 1, "left"), (2, 2, "right"), (2, 3, "y_min"), (2, 4, "y_max"),
        (2, 5, "z_min"), (2, 6, "z_max"), (2, 10, "grain_boundaries"),
        (3, 1000, "bulk"),
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
