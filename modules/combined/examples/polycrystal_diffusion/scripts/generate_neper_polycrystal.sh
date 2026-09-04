#!/usr/bin/env bash
set -euo pipefail

# Neper works most robustly at micrometer-scale coordinates. MOOSE scales the
# imported mesh by 1e-6 so that the calculation itself uses SI coordinates.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
EXAMPLE_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
REPO_DIR=$(cd -- "${EXAMPLE_DIR}/../../../.." && pwd)
MESH_DIR="${EXAMPLE_DIR}/mesh"
LOCAL_ENV="${REPO_DIR}/.neper-env"

NEPER=${NEPER:-$(command -v neper || true)}
GMSH=${GMSH:-$(command -v gmsh || true)}
if [[ -z "${NEPER}" && -x "${LOCAL_ENV}/bin/neper" ]]; then
  NEPER="${LOCAL_ENV}/bin/neper"
fi
if [[ -z "${GMSH}" && -x "${LOCAL_ENV}/bin/gmsh" ]]; then
  GMSH="${LOCAL_ENV}/bin/gmsh"
fi
if [[ -z "${NEPER}" || -z "${GMSH}" ]]; then
  echo "Neper and Gmsh are required. Set NEPER and GMSH to their executable paths." >&2
  exit 1
fi

CHARACTERISTIC_LENGTH=${1:-2.5}
NAME=${2:-neper_polycrystal}
mkdir -p "${MESH_DIR}"

echo "Using $(${NEPER} --version) from ${NEPER}"
echo "Using Gmsh $(${GMSH} --version) from ${GMSH}"

# 200 x 50 um gives 100 um^2 per grain, i.e. a 10 um characteristic
# area length. Seed id 1 makes the Voronoi realization reproducible.
OMP_NUM_THREADS=1 "${NEPER}" -T \
  -n 100 \
  -dim 2 \
  -domain 'square(200,50)' \
  -morpho voronoi \
  -id 1 \
  -regularization 1 \
  -o "${MESH_DIR}/${NAME}"

OMP_NUM_THREADS=1 "${NEPER}" -M "${MESH_DIR}/${NAME}.tess" \
  -gmsh "${GMSH}" \
  -tmp "${MESH_DIR}" \
  -elttype tri \
  -order 1 \
  -cl "${CHARACTERISTIC_LENGTH}" \
  -format msh:ascii \
  -o "${MESH_DIR}/${NAME}_raw"

python3 "${SCRIPT_DIR}/preprocess_neper_mesh.py" \
  "${MESH_DIR}/${NAME}_raw.msh" \
  "${MESH_DIR}/${NAME}.msh" \
  --stats "${MESH_DIR}/${NAME}_stats.json"
