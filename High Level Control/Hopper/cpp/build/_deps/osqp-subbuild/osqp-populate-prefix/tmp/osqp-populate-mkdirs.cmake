# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-src"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-build"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix/tmp"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix/src"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
