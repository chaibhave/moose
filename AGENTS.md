# How to work in this repo (Codex)

- Use the pre-activated conda/mamba environment.
- Clone MOOSE into $HOME/projects/moose if missing.
- Build only the requested or modified module(s) with `make -j$(nproc)` (METHOD=opt). 
- Do NOT run scripts/update_and_rebuild_*; we use the conda toolchain.
- If a dependency is missing, STOP and tell me the exact line to add to the setup script.
