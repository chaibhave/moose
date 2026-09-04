#!/usr/bin/env bash
set -euo pipefail

# Neper meshes in micrometers; the MOOSE input scales the final mesh to SI.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
EXAMPLE_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
REPO_DIR=$(cd -- "${EXAMPLE_DIR}/../../../.." && pwd)
MESH_DIR="${EXAMPLE_DIR}/mesh"
LOCAL_ENV="${REPO_DIR}/.neper-env"

NEPER=${NEPER:-$(command -v neper || true)}
GMSH=${GMSH:-$(command -v gmsh || true)}
[[ -n "${NEPER}" ]] || NEPER="${LOCAL_ENV}/bin/neper"
[[ -n "${GMSH}" ]] || GMSH="${LOCAL_ENV}/bin/gmsh"
if [[ ! -x "${NEPER}" || ! -x "${GMSH}" ]]; then
  echo "Neper and Gmsh are required. Set NEPER and GMSH to their executable paths." >&2
  exit 1
fi

CHARACTERISTIC_LENGTH=${1:-5}
NAME=${2:-neper_polycrystal_3d}
mkdir -p "${MESH_DIR}"

echo "Using $("${NEPER}" --version) from ${NEPER}"
echo "Using Gmsh $("${GMSH}" --version) from ${GMSH}"

OMP_NUM_THREADS=1 "${NEPER}" -T \
  -n 100 -dim 3 -domain 'cube(200,50,50)' -morpho voronoi \
  -id 1 -regularization 1 -o "${MESH_DIR}/${NAME}"

OMP_NUM_THREADS=1 "${NEPER}" -M "${MESH_DIR}/${NAME}.tess" \
  -gmsh "${GMSH}" -tmp "${MESH_DIR}" -elttype tet -order 1 \
  -cl "${CHARACTERISTIC_LENGTH}" -format msh:ascii -o "${MESH_DIR}/${NAME}_raw"

python3 "${SCRIPT_DIR}/preprocess_neper_mesh_3d.py" \
  "${MESH_DIR}/${NAME}_raw.msh" "${MESH_DIR}/${NAME}.msh" \
  --stats "${MESH_DIR}/${NAME}_stats.json"
