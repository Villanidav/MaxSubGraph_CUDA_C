# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /opt/install/cmake/3.27/bin/cmake

# The command to remove a file.
RM = /opt/install/cmake/3.27/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dvillani/subgraph/MaxSubGraph_CUDA_C

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dvillani/subgraph/MaxSubGraph_CUDA_C/build

# Include any dependencies generated for this target.
include CMakeFiles/test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/test.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test.dir/flags.make

CMakeFiles/test.dir/src/test.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/test.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/test.cpp
CMakeFiles/test.dir/src/test.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test.dir/src/test.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/test.cpp.o -MF CMakeFiles/test.dir/src/test.cpp.o.d -o CMakeFiles/test.dir/src/test.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/test.cpp

CMakeFiles/test.dir/src/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/test.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/test.cpp > CMakeFiles/test.dir/src/test.cpp.i

CMakeFiles/test.dir/src/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/test.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/test.cpp -o CMakeFiles/test.dir/src/test.cpp.s

CMakeFiles/test.dir/src/gen_bond_labels.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/gen_bond_labels.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp
CMakeFiles/test.dir/src/gen_bond_labels.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test.dir/src/gen_bond_labels.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/gen_bond_labels.cpp.o -MF CMakeFiles/test.dir/src/gen_bond_labels.cpp.o.d -o CMakeFiles/test.dir/src/gen_bond_labels.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp

CMakeFiles/test.dir/src/gen_bond_labels.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/gen_bond_labels.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp > CMakeFiles/test.dir/src/gen_bond_labels.cpp.i

CMakeFiles/test.dir/src/gen_bond_labels.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/gen_bond_labels.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp -o CMakeFiles/test.dir/src/gen_bond_labels.cpp.s

CMakeFiles/test.dir/src/gen_rotations.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/gen_rotations.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_rotations.cpp
CMakeFiles/test.dir/src/gen_rotations.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/test.dir/src/gen_rotations.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/gen_rotations.cpp.o -MF CMakeFiles/test.dir/src/gen_rotations.cpp.o.d -o CMakeFiles/test.dir/src/gen_rotations.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_rotations.cpp

CMakeFiles/test.dir/src/gen_rotations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/gen_rotations.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_rotations.cpp > CMakeFiles/test.dir/src/gen_rotations.cpp.i

CMakeFiles/test.dir/src/gen_rotations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/gen_rotations.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_rotations.cpp -o CMakeFiles/test.dir/src/gen_rotations.cpp.s

CMakeFiles/test.dir/src/select_vertex.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/select_vertex.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_vertex.cpp
CMakeFiles/test.dir/src/select_vertex.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/test.dir/src/select_vertex.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/select_vertex.cpp.o -MF CMakeFiles/test.dir/src/select_vertex.cpp.o.d -o CMakeFiles/test.dir/src/select_vertex.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_vertex.cpp

CMakeFiles/test.dir/src/select_vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/select_vertex.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_vertex.cpp > CMakeFiles/test.dir/src/select_vertex.cpp.i

CMakeFiles/test.dir/src/select_vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/select_vertex.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_vertex.cpp -o CMakeFiles/test.dir/src/select_vertex.cpp.s

CMakeFiles/test.dir/src/hood.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/hood.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/hood.cpp
CMakeFiles/test.dir/src/hood.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/test.dir/src/hood.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/hood.cpp.o -MF CMakeFiles/test.dir/src/hood.cpp.o.d -o CMakeFiles/test.dir/src/hood.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/hood.cpp

CMakeFiles/test.dir/src/hood.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/hood.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/hood.cpp > CMakeFiles/test.dir/src/hood.cpp.i

CMakeFiles/test.dir/src/hood.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/hood.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/hood.cpp -o CMakeFiles/test.dir/src/hood.cpp.s

CMakeFiles/test.dir/src/Label.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/Label.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/Label.cpp
CMakeFiles/test.dir/src/Label.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/test.dir/src/Label.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/Label.cpp.o -MF CMakeFiles/test.dir/src/Label.cpp.o.d -o CMakeFiles/test.dir/src/Label.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/Label.cpp

CMakeFiles/test.dir/src/Label.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/Label.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/Label.cpp > CMakeFiles/test.dir/src/Label.cpp.i

CMakeFiles/test.dir/src/Label.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/Label.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/Label.cpp -o CMakeFiles/test.dir/src/Label.cpp.s

CMakeFiles/test.dir/src/smiles_mcs.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/smiles_mcs.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp
CMakeFiles/test.dir/src/smiles_mcs.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/test.dir/src/smiles_mcs.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/smiles_mcs.cpp.o -MF CMakeFiles/test.dir/src/smiles_mcs.cpp.o.d -o CMakeFiles/test.dir/src/smiles_mcs.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp

CMakeFiles/test.dir/src/smiles_mcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/smiles_mcs.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp > CMakeFiles/test.dir/src/smiles_mcs.cpp.i

CMakeFiles/test.dir/src/smiles_mcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/smiles_mcs.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp -o CMakeFiles/test.dir/src/smiles_mcs.cpp.s

CMakeFiles/test.dir/src/calc_bound.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/calc_bound.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/calc_bound.cpp
CMakeFiles/test.dir/src/calc_bound.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/test.dir/src/calc_bound.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/calc_bound.cpp.o -MF CMakeFiles/test.dir/src/calc_bound.cpp.o.d -o CMakeFiles/test.dir/src/calc_bound.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/calc_bound.cpp

CMakeFiles/test.dir/src/calc_bound.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/calc_bound.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/calc_bound.cpp > CMakeFiles/test.dir/src/calc_bound.cpp.i

CMakeFiles/test.dir/src/calc_bound.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/calc_bound.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/calc_bound.cpp -o CMakeFiles/test.dir/src/calc_bound.cpp.s

CMakeFiles/test.dir/src/select_label.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/select_label.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_label.cpp
CMakeFiles/test.dir/src/select_label.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/test.dir/src/select_label.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/select_label.cpp.o -MF CMakeFiles/test.dir/src/select_label.cpp.o.d -o CMakeFiles/test.dir/src/select_label.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_label.cpp

CMakeFiles/test.dir/src/select_label.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/select_label.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_label.cpp > CMakeFiles/test.dir/src/select_label.cpp.i

CMakeFiles/test.dir/src/select_label.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/select_label.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/select_label.cpp -o CMakeFiles/test.dir/src/select_label.cpp.s

CMakeFiles/test.dir/src/gen_initial_labels.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/gen_initial_labels.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp
CMakeFiles/test.dir/src/gen_initial_labels.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/test.dir/src/gen_initial_labels.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/gen_initial_labels.cpp.o -MF CMakeFiles/test.dir/src/gen_initial_labels.cpp.o.d -o CMakeFiles/test.dir/src/gen_initial_labels.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp

CMakeFiles/test.dir/src/gen_initial_labels.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/gen_initial_labels.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp > CMakeFiles/test.dir/src/gen_initial_labels.cpp.i

CMakeFiles/test.dir/src/gen_initial_labels.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/gen_initial_labels.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp -o CMakeFiles/test.dir/src/gen_initial_labels.cpp.s

CMakeFiles/test.dir/src/gen_ring_classes.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/gen_ring_classes.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp
CMakeFiles/test.dir/src/gen_ring_classes.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/test.dir/src/gen_ring_classes.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/gen_ring_classes.cpp.o -MF CMakeFiles/test.dir/src/gen_ring_classes.cpp.o.d -o CMakeFiles/test.dir/src/gen_ring_classes.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp

CMakeFiles/test.dir/src/gen_ring_classes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/gen_ring_classes.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp > CMakeFiles/test.dir/src/gen_ring_classes.cpp.i

CMakeFiles/test.dir/src/gen_ring_classes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/gen_ring_classes.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp -o CMakeFiles/test.dir/src/gen_ring_classes.cpp.s

CMakeFiles/test.dir/src/search_mcs.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/search_mcs.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/search_mcs.cpp
CMakeFiles/test.dir/src/search_mcs.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/test.dir/src/search_mcs.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/search_mcs.cpp.o -MF CMakeFiles/test.dir/src/search_mcs.cpp.o.d -o CMakeFiles/test.dir/src/search_mcs.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/search_mcs.cpp

CMakeFiles/test.dir/src/search_mcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/search_mcs.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/search_mcs.cpp > CMakeFiles/test.dir/src/search_mcs.cpp.i

CMakeFiles/test.dir/src/search_mcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/search_mcs.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/search_mcs.cpp -o CMakeFiles/test.dir/src/search_mcs.cpp.s

CMakeFiles/test.dir/src/mol_mcs.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/mol_mcs.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mol_mcs.cpp
CMakeFiles/test.dir/src/mol_mcs.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/test.dir/src/mol_mcs.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/mol_mcs.cpp.o -MF CMakeFiles/test.dir/src/mol_mcs.cpp.o.d -o CMakeFiles/test.dir/src/mol_mcs.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mol_mcs.cpp

CMakeFiles/test.dir/src/mol_mcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/mol_mcs.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mol_mcs.cpp > CMakeFiles/test.dir/src/mol_mcs.cpp.i

CMakeFiles/test.dir/src/mol_mcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/mol_mcs.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mol_mcs.cpp -o CMakeFiles/test.dir/src/mol_mcs.cpp.s

CMakeFiles/test.dir/src/mc_split.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/mc_split.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mc_split.cpp
CMakeFiles/test.dir/src/mc_split.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/test.dir/src/mc_split.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/mc_split.cpp.o -MF CMakeFiles/test.dir/src/mc_split.cpp.o.d -o CMakeFiles/test.dir/src/mc_split.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mc_split.cpp

CMakeFiles/test.dir/src/mc_split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/mc_split.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mc_split.cpp > CMakeFiles/test.dir/src/mc_split.cpp.i

CMakeFiles/test.dir/src/mc_split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/mc_split.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/mc_split.cpp -o CMakeFiles/test.dir/src/mc_split.cpp.s

CMakeFiles/test.dir/src/g2mol.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/g2mol.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/g2mol.cpp
CMakeFiles/test.dir/src/g2mol.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/test.dir/src/g2mol.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/g2mol.cpp.o -MF CMakeFiles/test.dir/src/g2mol.cpp.o.d -o CMakeFiles/test.dir/src/g2mol.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/g2mol.cpp

CMakeFiles/test.dir/src/g2mol.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/g2mol.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/g2mol.cpp > CMakeFiles/test.dir/src/g2mol.cpp.i

CMakeFiles/test.dir/src/g2mol.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/g2mol.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/g2mol.cpp -o CMakeFiles/test.dir/src/g2mol.cpp.s

CMakeFiles/test.dir/src/gpu_mc_split.cu.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/gpu_mc_split.cu.o: CMakeFiles/test.dir/includes_CUDA.rsp
CMakeFiles/test.dir/src/gpu_mc_split.cu.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gpu_mc_split.cu
CMakeFiles/test.dir/src/gpu_mc_split.cu.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CUDA object CMakeFiles/test.dir/src/gpu_mc_split.cu.o"
	/opt/install/cuda/11.7/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/test.dir/src/gpu_mc_split.cu.o -MF CMakeFiles/test.dir/src/gpu_mc_split.cu.o.d -x cu -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/gpu_mc_split.cu -o CMakeFiles/test.dir/src/gpu_mc_split.cu.o

CMakeFiles/test.dir/src/gpu_mc_split.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CUDA source to CMakeFiles/test.dir/src/gpu_mc_split.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/test.dir/src/gpu_mc_split.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CUDA source to assembly CMakeFiles/test.dir/src/gpu_mc_split.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/test.dir/src/pair_vertex.cpp.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/pair_vertex.cpp.o: /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/pair_vertex.cpp
CMakeFiles/test.dir/src/pair_vertex.cpp.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object CMakeFiles/test.dir/src/pair_vertex.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/pair_vertex.cpp.o -MF CMakeFiles/test.dir/src/pair_vertex.cpp.o.d -o CMakeFiles/test.dir/src/pair_vertex.cpp.o -c /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/pair_vertex.cpp

CMakeFiles/test.dir/src/pair_vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test.dir/src/pair_vertex.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/pair_vertex.cpp > CMakeFiles/test.dir/src/pair_vertex.cpp.i

CMakeFiles/test.dir/src/pair_vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/pair_vertex.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dvillani/subgraph/MaxSubGraph_CUDA_C/src/pair_vertex.cpp -o CMakeFiles/test.dir/src/pair_vertex.cpp.s

# Object files for target test
test_OBJECTS = \
"CMakeFiles/test.dir/src/test.cpp.o" \
"CMakeFiles/test.dir/src/gen_bond_labels.cpp.o" \
"CMakeFiles/test.dir/src/gen_rotations.cpp.o" \
"CMakeFiles/test.dir/src/select_vertex.cpp.o" \
"CMakeFiles/test.dir/src/hood.cpp.o" \
"CMakeFiles/test.dir/src/Label.cpp.o" \
"CMakeFiles/test.dir/src/smiles_mcs.cpp.o" \
"CMakeFiles/test.dir/src/calc_bound.cpp.o" \
"CMakeFiles/test.dir/src/select_label.cpp.o" \
"CMakeFiles/test.dir/src/gen_initial_labels.cpp.o" \
"CMakeFiles/test.dir/src/gen_ring_classes.cpp.o" \
"CMakeFiles/test.dir/src/search_mcs.cpp.o" \
"CMakeFiles/test.dir/src/mol_mcs.cpp.o" \
"CMakeFiles/test.dir/src/mc_split.cpp.o" \
"CMakeFiles/test.dir/src/g2mol.cpp.o" \
"CMakeFiles/test.dir/src/gpu_mc_split.cu.o" \
"CMakeFiles/test.dir/src/pair_vertex.cpp.o"

# External object files for target test
test_EXTERNAL_OBJECTS =

test: CMakeFiles/test.dir/src/test.cpp.o
test: CMakeFiles/test.dir/src/gen_bond_labels.cpp.o
test: CMakeFiles/test.dir/src/gen_rotations.cpp.o
test: CMakeFiles/test.dir/src/select_vertex.cpp.o
test: CMakeFiles/test.dir/src/hood.cpp.o
test: CMakeFiles/test.dir/src/Label.cpp.o
test: CMakeFiles/test.dir/src/smiles_mcs.cpp.o
test: CMakeFiles/test.dir/src/calc_bound.cpp.o
test: CMakeFiles/test.dir/src/select_label.cpp.o
test: CMakeFiles/test.dir/src/gen_initial_labels.cpp.o
test: CMakeFiles/test.dir/src/gen_ring_classes.cpp.o
test: CMakeFiles/test.dir/src/search_mcs.cpp.o
test: CMakeFiles/test.dir/src/mol_mcs.cpp.o
test: CMakeFiles/test.dir/src/mc_split.cpp.o
test: CMakeFiles/test.dir/src/g2mol.cpp.o
test: CMakeFiles/test.dir/src/gpu_mc_split.cu.o
test: CMakeFiles/test.dir/src/pair_vertex.cpp.o
test: CMakeFiles/test.dir/build.make
test: CMakeFiles/test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Linking CXX executable test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test.dir/build: test
.PHONY : CMakeFiles/test.dir/build

CMakeFiles/test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test.dir/clean

CMakeFiles/test.dir/depend:
	cd /home/dvillani/subgraph/MaxSubGraph_CUDA_C/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dvillani/subgraph/MaxSubGraph_CUDA_C /home/dvillani/subgraph/MaxSubGraph_CUDA_C /home/dvillani/subgraph/MaxSubGraph_CUDA_C/build /home/dvillani/subgraph/MaxSubGraph_CUDA_C/build /home/dvillani/subgraph/MaxSubGraph_CUDA_C/build/CMakeFiles/test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/test.dir/depend

