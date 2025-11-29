# ImmSim-FLAMEGPU2-model

Immunitary System simulation model running in CUDA/C++ using the FLAMEGPU2 framework.
This model is based on the ImmSim model from Celada and Seiden published in 1992.

## Requirements

Building FLAME GPU from source has the following requirements. There are also optional dependencies which are required for some components, such as Documentation or Python bindings, however these are not strictly required, and are not required for this standalone example.

Building FLAME GPU has the following requirements. There are also optional dependencies which are required for some components, such as Documentation or Python bindings.

+ [CMake](https://cmake.org/download/) `>= 3.18`
  + `>= 3.20` if building python bindings using a multi-config generator (Visual Studio, Eclipse or Ninja Multi-Config)
+ [CUDA](https://developer.nvidia.com/cuda-downloads) `>= 11.0` and a [Compute Capability](https://developer.nvidia.com/cuda-gpus) `>= 3.5` NVIDIA GPU.
+ C++17 capable C++ compiler (host), compatible with the installed CUDA version
  + [Microsoft Visual Studio 2019 or 2022](https://visualstudio.microsoft.com/) (Windows)
    + *Note:* Visual Studio must be installed before the CUDA toolkit is installed. See the [CUDA installation guide for Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html) for more information.
  + [make](https://www.gnu.org/software/make/) and [GCC](https://gcc.gnu.org/) `>= 8.1` (Linux)
+ [git](https://git-scm.com/)

Optionally:

+ [cpplint](https://github.com/cpplint/cpplint) for linting code
+ [Doxygen](http://www.doxygen.nl/) to build the documentation
+ [Python](https://www.python.org/) `>= 3.7` for python integration
  + With `setuptools`, `wheel`, `build` and optionally `venv` python packages installed
+ [swig](http://www.swig.org/) `>= 4.0.2` for python integration
  + Swig `4.x` will be automatically downloaded by CMake if not provided (if possible).
+ [FLAMEGPU2-visualiser](https://github.com/FLAMEGPU/FLAMEGPU2-visualiser) dependencies (fetched if possible)
  + [SDL](https://www.libsdl.org/)
  + [GLM](http://glm.g-truc.net/) *(consistent C++/GLSL vector maths functionality)*
  + [GLEW](http://glew.sourceforge.net/) *(GL extension loader)*
  + [FreeType](http://www.freetype.org/)  *(font loading)*
  + [DevIL](http://openil.sourceforge.net/)  *(image loading)*
  + [Fontconfig](https://www.fontconfig.org/)  *(Linux only, font detection)*


## Building with CMake

Building via CMake is a three step process, with slight differences depending on your platform.

1. Create a build directory for an out-of tree build
2. Configure CMake into the build directory
    + Using the CMake GUI or CLI tools
    + Specifying build options such as the CUDA Compute Capabilities to target, the inclusion of Visualisation or Python components, or performance impacting features such as `FLAMEGPU_SEATBELTS`. See [CMake Configuration Options](#CMake-Configuration-Options) for details of the available configuration options
    + CMake will automatically find and select compilers
3. Build compilation targets using the configured build system

### Linux

CMake instructions to build and run (from root dir):

## CMake setup
`cmake --fresh -S . -B build/` to plot entity graphs you can also add `-DPLOT=ON`

## Build
`cmake --build build/ --target ImmSim` to speedup with multithreading you can also add `-j$(nproc)`

## Execute
`./build/bin/Release/ImmSim`


## License

FLAME GPU is distributed under the [MIT Licence](https://github.com/FLAMEGPU/FLAMEGPU2/blob/master/LICENSE.md).
