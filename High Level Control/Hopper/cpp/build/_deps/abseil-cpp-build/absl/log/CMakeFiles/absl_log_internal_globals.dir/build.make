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
include _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/flags.make

_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o: _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/flags.make
_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o: _deps/abseil-cpp-src/absl/log/internal/globals.cc
_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o: _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o -MF CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o.d -o CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o -c "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/log/internal/globals.cc"

_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.i"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/log/internal/globals.cc" > CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.i

_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.s"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/log/internal/globals.cc" -o CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.s

# Object files for target absl_log_internal_globals
absl_log_internal_globals_OBJECTS = \
"CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o"

# External object files for target absl_log_internal_globals
absl_log_internal_globals_EXTERNAL_OBJECTS =

_deps/abseil-cpp-build/absl/log/libabsl_log_internal_globals.a: _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/internal/globals.cc.o
_deps/abseil-cpp-build/absl/log/libabsl_log_internal_globals.a: _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/build.make
_deps/abseil-cpp-build/absl/log/libabsl_log_internal_globals.a: _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libabsl_log_internal_globals.a"
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" && $(CMAKE_COMMAND) -P CMakeFiles/absl_log_internal_globals.dir/cmake_clean_target.cmake
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/absl_log_internal_globals.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/build: _deps/abseil-cpp-build/absl/log/libabsl_log_internal_globals.a
.PHONY : _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/build

_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/clean:
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" && $(CMAKE_COMMAND) -P CMakeFiles/absl_log_internal_globals.dir/cmake_clean.cmake
.PHONY : _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/clean

_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/depend:
	cd "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-src/absl/log" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log" "/home/noelcs/repos/ThinkingThroughThings/High Level Control/Hopper/cpp/build/_deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : _deps/abseil-cpp-build/absl/log/CMakeFiles/absl_log_internal_globals.dir/depend

