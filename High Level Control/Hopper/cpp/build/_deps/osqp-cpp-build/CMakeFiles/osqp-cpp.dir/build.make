# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build"

# Include any dependencies generated for this target.
include _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/flags.make

_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o: _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/flags.make
_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o: _deps/osqp-cpp-src/src/osqp++.cc
_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o: _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o -MF CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o.d -o CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o -c "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-src/src/osqp++.cc"

_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/osqp-cpp.dir/src/osqp++.cc.i"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-src/src/osqp++.cc" > CMakeFiles/osqp-cpp.dir/src/osqp++.cc.i

_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/osqp-cpp.dir/src/osqp++.cc.s"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-src/src/osqp++.cc" -o CMakeFiles/osqp-cpp.dir/src/osqp++.cc.s

# Object files for target osqp-cpp
osqp__cpp_OBJECTS = \
"CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o"

# External object files for target osqp-cpp
osqp__cpp_EXTERNAL_OBJECTS =

_deps/osqp-cpp-build/libosqp-cpp.a: _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/src/osqp++.cc.o
_deps/osqp-cpp-build/libosqp-cpp.a: _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/build.make
_deps/osqp-cpp-build/libosqp-cpp.a: _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libosqp-cpp.a"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" && $(CMAKE_COMMAND) -P CMakeFiles/osqp-cpp.dir/cmake_clean_target.cmake
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/osqp-cpp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/build: _deps/osqp-cpp-build/libosqp-cpp.a
.PHONY : _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/build

_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/clean:
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" && $(CMAKE_COMMAND) -P CMakeFiles/osqp-cpp.dir/cmake_clean.cmake
.PHONY : _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/clean

_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/depend:
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-src" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : _deps/osqp-cpp-build/CMakeFiles/osqp-cpp.dir/depend

