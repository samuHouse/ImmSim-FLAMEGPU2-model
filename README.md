# ImmSim-FLAMEGPU2-model
Immunitary System simulation model running in CUDA/C++ using the FLAMEGPU2 framework.
This model is based on the ImmSim model from Celada and Seiden published in 1992.

CMake instructions to build and run (from root dir):

## CMake setup
`cmake --fresh -S . -B build/`

## Build
`cmake --build build/ --target ImmSim` to speedup with multithreading you can also add `-j$(nproc)`

## Execute
`./build/bin/Release/ImmSim`
