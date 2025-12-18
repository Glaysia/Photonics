#!/usr/bin/env bash

set -euo pipefail
# Default to 12 ranks unless overridden.
NP=40
# Default to 1 OpenMP thread per rank to avoid oversubscription.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# Ensure Meep built under external/meep/install is found at runtime.
export LD_LIBRARY_PATH="/home/harry/CLionProjects/Photo/Photonics/external/meep/install/lib:${LD_LIBRARY_PATH:-}"

exec mpirun --oversubscribe -np "${NP}" ./build/Photonics "$@"

# exec ./build/Photonics "$@"
