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
include _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/flags.make

_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o: _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/flags.make
_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o: _deps/abseil-cpp-src/absl/strings/internal/cordz_handle.cc
_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o: _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o -MF CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o.d -o CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o -c "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/strings/internal/cordz_handle.cc"

_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.i"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/strings/internal/cordz_handle.cc" > CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.i

_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.s"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/strings/internal/cordz_handle.cc" -o CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.s

# Object files for target absl_cordz_handle
absl_cordz_handle_OBJECTS = \
"CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o"

# External object files for target absl_cordz_handle
absl_cordz_handle_EXTERNAL_OBJECTS =

_deps/abseil-cpp-build/absl/strings/libabsl_cordz_handle.a: _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/internal/cordz_handle.cc.o
_deps/abseil-cpp-build/absl/strings/libabsl_cordz_handle.a: _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/build.make
_deps/abseil-cpp-build/absl/strings/libabsl_cordz_handle.a: _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libabsl_cordz_handle.a"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" && $(CMAKE_COMMAND) -P CMakeFiles/absl_cordz_handle.dir/cmake_clean_target.cmake
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/absl_cordz_handle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/build: _deps/abseil-cpp-build/absl/strings/libabsl_cordz_handle.a
.PHONY : _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/build

_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/clean:
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" && $(CMAKE_COMMAND) -P CMakeFiles/absl_cordz_handle.dir/cmake_clean.cmake
.PHONY : _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/clean

_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/depend:
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/strings" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : _deps/abseil-cpp-build/absl/strings/CMakeFiles/absl_cordz_handle.dir/depend

