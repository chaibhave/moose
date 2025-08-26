# How to work in this repo (Codex)

- Use the pre-activated conda/mamba environment.Run "micromamba load moose" before trying to build any code or run tests/input files.
- Clone MOOSE into $HOME/projects/moose if missing.
- Build a framework or a specific module before attempting to run its tests.
- Build only the requested or modified module(s) with `make -j$(nproc)` (METHOD=opt). 
- Do NOT run scripts/update_and_rebuild_*; we use the conda toolchain.
- If a dependency is missing, STOP and tell me the exact line to add to the setup script.
