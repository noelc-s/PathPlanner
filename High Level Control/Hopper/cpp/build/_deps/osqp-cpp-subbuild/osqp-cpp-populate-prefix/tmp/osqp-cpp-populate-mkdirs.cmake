# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-src"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix/tmp"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix/src/osqp-cpp-populate-stamp"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix/src"
  "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix/src/osqp-cpp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix/src/osqp-cpp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-subbuild/osqp-cpp-populate-prefix/src/osqp-cpp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
