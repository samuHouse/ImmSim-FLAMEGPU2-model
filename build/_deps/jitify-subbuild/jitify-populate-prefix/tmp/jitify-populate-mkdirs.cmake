# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-src/jitify")
  file(MAKE_DIRECTORY "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-src/jitify")
endif()
file(MAKE_DIRECTORY
  "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-build"
  "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix"
  "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix/tmp"
  "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix/src/jitify-populate-stamp"
  "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix/src"
  "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix/src/jitify-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix/src/jitify-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/local/STUDENTI/samuele.casadei3/TesiImmSim/build/_deps/jitify-subbuild/jitify-populate-prefix/src/jitify-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
