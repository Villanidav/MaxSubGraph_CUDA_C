# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.29

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.29.0/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.29.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build

# Include any dependencies generated for this target.
include CMakeFiles/test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/test.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test.dir/flags.make

CMakeFiles/test.dir/src/test.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/test.cpp.o: /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/test.cpp
CMakeFiles/test.dir/src/test.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test.dir/src/test.cpp.o"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/test.cpp.o -MF CMakeFiles/test.dir/src/test.cpp.o.d -o CMakeFiles/test.dir/src/test.cpp.o -c /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/test.cpp

CMakeFiles/test.dir/src/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/test.cpp.i"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/test.cpp > CMakeFiles/test.dir/src/test.cpp.i

CMakeFiles/test.dir/src/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/test.cpp.s"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/test.cpp -o CMakeFiles/test.dir/src/test.cpp.s

CMakeFiles/test.dir/src/gen_rotations.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/gen_rotations.cpp.o: /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/gen_rotations.cpp
CMakeFiles/test.dir/src/gen_rotations.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test.dir/src/gen_rotations.cpp.o"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/gen_rotations.cpp.o -MF CMakeFiles/test.dir/src/gen_rotations.cpp.o.d -o CMakeFiles/test.dir/src/gen_rotations.cpp.o -c /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/gen_rotations.cpp

CMakeFiles/test.dir/src/gen_rotations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/gen_rotations.cpp.i"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/gen_rotations.cpp > CMakeFiles/test.dir/src/gen_rotations.cpp.i

CMakeFiles/test.dir/src/gen_rotations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/gen_rotations.cpp.s"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/gen_rotations.cpp -o CMakeFiles/test.dir/src/gen_rotations.cpp.s

CMakeFiles/test.dir/src/select_vertex.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/select_vertex.cpp.o: /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/select_vertex.cpp
CMakeFiles/test.dir/src/select_vertex.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/test.dir/src/select_vertex.cpp.o"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/select_vertex.cpp.o -MF CMakeFiles/test.dir/src/select_vertex.cpp.o.d -o CMakeFiles/test.dir/src/select_vertex.cpp.o -c /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/select_vertex.cpp

CMakeFiles/test.dir/src/select_vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/select_vertex.cpp.i"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/select_vertex.cpp > CMakeFiles/test.dir/src/select_vertex.cpp.i

CMakeFiles/test.dir/src/select_vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/select_vertex.cpp.s"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/select_vertex.cpp -o CMakeFiles/test.dir/src/select_vertex.cpp.s

CMakeFiles/test.dir/src/hood.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/hood.cpp.o: /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/hood.cpp
CMakeFiles/test.dir/src/hood.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/test.dir/src/hood.cpp.o"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/hood.cpp.o -MF CMakeFiles/test.dir/src/hood.cpp.o.d -o CMakeFiles/test.dir/src/hood.cpp.o -c /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/hood.cpp

CMakeFiles/test.dir/src/hood.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/hood.cpp.i"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/hood.cpp > CMakeFiles/test.dir/src/hood.cpp.i

CMakeFiles/test.dir/src/hood.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/hood.cpp.s"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/hood.cpp -o CMakeFiles/test.dir/src/hood.cpp.s

CMakeFiles/test.dir/src/Label.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/Label.cpp.o: /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/Label.cpp
CMakeFiles/test.dir/src/Label.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/test.dir/src/Label.cpp.o"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/Label.cpp.o -MF CMakeFiles/test.dir/src/Label.cpp.o.d -o CMakeFiles/test.dir/src/Label.cpp.o -c /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/Label.cpp

CMakeFiles/test.dir/src/Label.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/Label.cpp.i"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/Label.cpp > CMakeFiles/test.dir/src/Label.cpp.i

CMakeFiles/test.dir/src/Label.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/Label.cpp.s"
	/opt/homebrew/bin/aarch64-apple-darwin23-g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/src/Label.cpp -o CMakeFiles/test.dir/src/Label.cpp.s

# Object files for target test
test_OBJECTS = \
"CMakeFiles/test.dir/src/test.cpp.o" \
"CMakeFiles/test.dir/src/gen_rotations.cpp.o" \
"CMakeFiles/test.dir/src/select_vertex.cpp.o" \
"CMakeFiles/test.dir/src/hood.cpp.o" \
"CMakeFiles/test.dir/src/Label.cpp.o"

# External object files for target test
test_EXTERNAL_OBJECTS =

test: CMakeFiles/test.dir/src/test.cpp.o
test: CMakeFiles/test.dir/src/gen_rotations.cpp.o
test: CMakeFiles/test.dir/src/select_vertex.cpp.o
test: CMakeFiles/test.dir/src/hood.cpp.o
test: CMakeFiles/test.dir/src/Label.cpp.o
test: CMakeFiles/test.dir/build.make
test: CMakeFiles/test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test.dir/build: test
.PHONY : CMakeFiles/test.dir/build

CMakeFiles/test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test.dir/clean

CMakeFiles/test.dir/depend:
	cd /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build /Users/francescovirgulti/Desktop/MaxSubGraph_CUDA_C/build/CMakeFiles/test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/test.dir/depend

