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
CMAKE_SOURCE_DIR = /home/fvirgulti/MaxSubGraph_CUDA_C

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/fvirgulti/MaxSubGraph_CUDA_C/build

# Include any dependencies generated for this target.
include CMakeFiles/testFra.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/testFra.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/testFra.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testFra.dir/flags.make

CMakeFiles/testFra.dir/src/testFra.cu.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/testFra.cu.o: CMakeFiles/testFra.dir/includes_CUDA.rsp
CMakeFiles/testFra.dir/src/testFra.cu.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/testFra.cu
CMakeFiles/testFra.dir/src/testFra.cu.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/testFra.dir/src/testFra.cu.o"
	/opt/install/cuda/11.7/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/testFra.cu.o -MF CMakeFiles/testFra.dir/src/testFra.cu.o.d -x cu -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/testFra.cu -o CMakeFiles/testFra.dir/src/testFra.cu.o

CMakeFiles/testFra.dir/src/testFra.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CUDA source to CMakeFiles/testFra.dir/src/testFra.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/testFra.dir/src/testFra.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CUDA source to assembly CMakeFiles/testFra.dir/src/testFra.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp
CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o -MF CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o.d -o CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp

CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp > CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.i

CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_bond_labels.cpp -o CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.s

CMakeFiles/testFra.dir/src/gen_rotations.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/gen_rotations.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_rotations.cpp
CMakeFiles/testFra.dir/src/gen_rotations.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/testFra.dir/src/gen_rotations.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/gen_rotations.cpp.o -MF CMakeFiles/testFra.dir/src/gen_rotations.cpp.o.d -o CMakeFiles/testFra.dir/src/gen_rotations.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_rotations.cpp

CMakeFiles/testFra.dir/src/gen_rotations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/gen_rotations.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_rotations.cpp > CMakeFiles/testFra.dir/src/gen_rotations.cpp.i

CMakeFiles/testFra.dir/src/gen_rotations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/gen_rotations.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_rotations.cpp -o CMakeFiles/testFra.dir/src/gen_rotations.cpp.s

CMakeFiles/testFra.dir/src/select_vertex.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/select_vertex.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_vertex.cpp
CMakeFiles/testFra.dir/src/select_vertex.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/testFra.dir/src/select_vertex.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/select_vertex.cpp.o -MF CMakeFiles/testFra.dir/src/select_vertex.cpp.o.d -o CMakeFiles/testFra.dir/src/select_vertex.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_vertex.cpp

CMakeFiles/testFra.dir/src/select_vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/select_vertex.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_vertex.cpp > CMakeFiles/testFra.dir/src/select_vertex.cpp.i

CMakeFiles/testFra.dir/src/select_vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/select_vertex.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_vertex.cpp -o CMakeFiles/testFra.dir/src/select_vertex.cpp.s

CMakeFiles/testFra.dir/src/hood.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/hood.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/hood.cpp
CMakeFiles/testFra.dir/src/hood.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/testFra.dir/src/hood.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/hood.cpp.o -MF CMakeFiles/testFra.dir/src/hood.cpp.o.d -o CMakeFiles/testFra.dir/src/hood.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/hood.cpp

CMakeFiles/testFra.dir/src/hood.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/hood.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/hood.cpp > CMakeFiles/testFra.dir/src/hood.cpp.i

CMakeFiles/testFra.dir/src/hood.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/hood.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/hood.cpp -o CMakeFiles/testFra.dir/src/hood.cpp.s

CMakeFiles/testFra.dir/src/Label.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/Label.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/Label.cpp
CMakeFiles/testFra.dir/src/Label.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/testFra.dir/src/Label.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/Label.cpp.o -MF CMakeFiles/testFra.dir/src/Label.cpp.o.d -o CMakeFiles/testFra.dir/src/Label.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/Label.cpp

CMakeFiles/testFra.dir/src/Label.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/Label.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/Label.cpp > CMakeFiles/testFra.dir/src/Label.cpp.i

CMakeFiles/testFra.dir/src/Label.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/Label.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/Label.cpp -o CMakeFiles/testFra.dir/src/Label.cpp.s

CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp
CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o -MF CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o.d -o CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp

CMakeFiles/testFra.dir/src/smiles_mcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/smiles_mcs.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp > CMakeFiles/testFra.dir/src/smiles_mcs.cpp.i

CMakeFiles/testFra.dir/src/smiles_mcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/smiles_mcs.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/smiles_mcs.cpp -o CMakeFiles/testFra.dir/src/smiles_mcs.cpp.s

CMakeFiles/testFra.dir/src/calc_bound.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/calc_bound.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/calc_bound.cpp
CMakeFiles/testFra.dir/src/calc_bound.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/testFra.dir/src/calc_bound.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/calc_bound.cpp.o -MF CMakeFiles/testFra.dir/src/calc_bound.cpp.o.d -o CMakeFiles/testFra.dir/src/calc_bound.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/calc_bound.cpp

CMakeFiles/testFra.dir/src/calc_bound.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/calc_bound.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/calc_bound.cpp > CMakeFiles/testFra.dir/src/calc_bound.cpp.i

CMakeFiles/testFra.dir/src/calc_bound.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/calc_bound.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/calc_bound.cpp -o CMakeFiles/testFra.dir/src/calc_bound.cpp.s

CMakeFiles/testFra.dir/src/select_label.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/select_label.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_label.cpp
CMakeFiles/testFra.dir/src/select_label.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/testFra.dir/src/select_label.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/select_label.cpp.o -MF CMakeFiles/testFra.dir/src/select_label.cpp.o.d -o CMakeFiles/testFra.dir/src/select_label.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_label.cpp

CMakeFiles/testFra.dir/src/select_label.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/select_label.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_label.cpp > CMakeFiles/testFra.dir/src/select_label.cpp.i

CMakeFiles/testFra.dir/src/select_label.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/select_label.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/select_label.cpp -o CMakeFiles/testFra.dir/src/select_label.cpp.s

CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp
CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o -MF CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o.d -o CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp

CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp > CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.i

CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_initial_labels.cpp -o CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.s

CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp
CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o -MF CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o.d -o CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp

CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp > CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.i

CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/gen_ring_classes.cpp -o CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.s

CMakeFiles/testFra.dir/src/search_mcs.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/search_mcs.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/search_mcs.cpp
CMakeFiles/testFra.dir/src/search_mcs.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/testFra.dir/src/search_mcs.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/search_mcs.cpp.o -MF CMakeFiles/testFra.dir/src/search_mcs.cpp.o.d -o CMakeFiles/testFra.dir/src/search_mcs.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/search_mcs.cpp

CMakeFiles/testFra.dir/src/search_mcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/search_mcs.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/search_mcs.cpp > CMakeFiles/testFra.dir/src/search_mcs.cpp.i

CMakeFiles/testFra.dir/src/search_mcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/search_mcs.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/search_mcs.cpp -o CMakeFiles/testFra.dir/src/search_mcs.cpp.s

CMakeFiles/testFra.dir/src/mol_mcs.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/mol_mcs.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/mol_mcs.cpp
CMakeFiles/testFra.dir/src/mol_mcs.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/testFra.dir/src/mol_mcs.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/mol_mcs.cpp.o -MF CMakeFiles/testFra.dir/src/mol_mcs.cpp.o.d -o CMakeFiles/testFra.dir/src/mol_mcs.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/mol_mcs.cpp

CMakeFiles/testFra.dir/src/mol_mcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/mol_mcs.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/mol_mcs.cpp > CMakeFiles/testFra.dir/src/mol_mcs.cpp.i

CMakeFiles/testFra.dir/src/mol_mcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/mol_mcs.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/mol_mcs.cpp -o CMakeFiles/testFra.dir/src/mol_mcs.cpp.s

CMakeFiles/testFra.dir/src/mc_split.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/mc_split.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/mc_split.cpp
CMakeFiles/testFra.dir/src/mc_split.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/testFra.dir/src/mc_split.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/mc_split.cpp.o -MF CMakeFiles/testFra.dir/src/mc_split.cpp.o.d -o CMakeFiles/testFra.dir/src/mc_split.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/mc_split.cpp

CMakeFiles/testFra.dir/src/mc_split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/mc_split.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/mc_split.cpp > CMakeFiles/testFra.dir/src/mc_split.cpp.i

CMakeFiles/testFra.dir/src/mc_split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/mc_split.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/mc_split.cpp -o CMakeFiles/testFra.dir/src/mc_split.cpp.s

CMakeFiles/testFra.dir/src/g2mol.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/g2mol.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/g2mol.cpp
CMakeFiles/testFra.dir/src/g2mol.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/testFra.dir/src/g2mol.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/g2mol.cpp.o -MF CMakeFiles/testFra.dir/src/g2mol.cpp.o.d -o CMakeFiles/testFra.dir/src/g2mol.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/g2mol.cpp

CMakeFiles/testFra.dir/src/g2mol.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/g2mol.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/g2mol.cpp > CMakeFiles/testFra.dir/src/g2mol.cpp.i

CMakeFiles/testFra.dir/src/g2mol.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/g2mol.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/g2mol.cpp -o CMakeFiles/testFra.dir/src/g2mol.cpp.s

CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o: CMakeFiles/testFra.dir/includes_CUDA.rsp
CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/gpu_mc_split.cu
CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CUDA object CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o"
	/opt/install/cuda/11.7/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o -MF CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o.d -x cu -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/gpu_mc_split.cu -o CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o

CMakeFiles/testFra.dir/src/gpu_mc_split.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CUDA source to CMakeFiles/testFra.dir/src/gpu_mc_split.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/testFra.dir/src/gpu_mc_split.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CUDA source to assembly CMakeFiles/testFra.dir/src/gpu_mc_split.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/testFra.dir/src/pair_vertex.cpp.o: CMakeFiles/testFra.dir/flags.make
CMakeFiles/testFra.dir/src/pair_vertex.cpp.o: /home/fvirgulti/MaxSubGraph_CUDA_C/src/pair_vertex.cpp
CMakeFiles/testFra.dir/src/pair_vertex.cpp.o: CMakeFiles/testFra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object CMakeFiles/testFra.dir/src/pair_vertex.cpp.o"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testFra.dir/src/pair_vertex.cpp.o -MF CMakeFiles/testFra.dir/src/pair_vertex.cpp.o.d -o CMakeFiles/testFra.dir/src/pair_vertex.cpp.o -c /home/fvirgulti/MaxSubGraph_CUDA_C/src/pair_vertex.cpp

CMakeFiles/testFra.dir/src/pair_vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testFra.dir/src/pair_vertex.cpp.i"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fvirgulti/MaxSubGraph_CUDA_C/src/pair_vertex.cpp > CMakeFiles/testFra.dir/src/pair_vertex.cpp.i

CMakeFiles/testFra.dir/src/pair_vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testFra.dir/src/pair_vertex.cpp.s"
	/opt/install/gcc/11.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fvirgulti/MaxSubGraph_CUDA_C/src/pair_vertex.cpp -o CMakeFiles/testFra.dir/src/pair_vertex.cpp.s

# Object files for target testFra
testFra_OBJECTS = \
"CMakeFiles/testFra.dir/src/testFra.cu.o" \
"CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o" \
"CMakeFiles/testFra.dir/src/gen_rotations.cpp.o" \
"CMakeFiles/testFra.dir/src/select_vertex.cpp.o" \
"CMakeFiles/testFra.dir/src/hood.cpp.o" \
"CMakeFiles/testFra.dir/src/Label.cpp.o" \
"CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o" \
"CMakeFiles/testFra.dir/src/calc_bound.cpp.o" \
"CMakeFiles/testFra.dir/src/select_label.cpp.o" \
"CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o" \
"CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o" \
"CMakeFiles/testFra.dir/src/search_mcs.cpp.o" \
"CMakeFiles/testFra.dir/src/mol_mcs.cpp.o" \
"CMakeFiles/testFra.dir/src/mc_split.cpp.o" \
"CMakeFiles/testFra.dir/src/g2mol.cpp.o" \
"CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o" \
"CMakeFiles/testFra.dir/src/pair_vertex.cpp.o"

# External object files for target testFra
testFra_EXTERNAL_OBJECTS =

testFra: CMakeFiles/testFra.dir/src/testFra.cu.o
testFra: CMakeFiles/testFra.dir/src/gen_bond_labels.cpp.o
testFra: CMakeFiles/testFra.dir/src/gen_rotations.cpp.o
testFra: CMakeFiles/testFra.dir/src/select_vertex.cpp.o
testFra: CMakeFiles/testFra.dir/src/hood.cpp.o
testFra: CMakeFiles/testFra.dir/src/Label.cpp.o
testFra: CMakeFiles/testFra.dir/src/smiles_mcs.cpp.o
testFra: CMakeFiles/testFra.dir/src/calc_bound.cpp.o
testFra: CMakeFiles/testFra.dir/src/select_label.cpp.o
testFra: CMakeFiles/testFra.dir/src/gen_initial_labels.cpp.o
testFra: CMakeFiles/testFra.dir/src/gen_ring_classes.cpp.o
testFra: CMakeFiles/testFra.dir/src/search_mcs.cpp.o
testFra: CMakeFiles/testFra.dir/src/mol_mcs.cpp.o
testFra: CMakeFiles/testFra.dir/src/mc_split.cpp.o
testFra: CMakeFiles/testFra.dir/src/g2mol.cpp.o
testFra: CMakeFiles/testFra.dir/src/gpu_mc_split.cu.o
testFra: CMakeFiles/testFra.dir/src/pair_vertex.cpp.o
testFra: CMakeFiles/testFra.dir/build.make
testFra: CMakeFiles/testFra.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Linking CXX executable testFra"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testFra.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testFra.dir/build: testFra
.PHONY : CMakeFiles/testFra.dir/build

CMakeFiles/testFra.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testFra.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testFra.dir/clean

CMakeFiles/testFra.dir/depend:
	cd /home/fvirgulti/MaxSubGraph_CUDA_C/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fvirgulti/MaxSubGraph_CUDA_C /home/fvirgulti/MaxSubGraph_CUDA_C /home/fvirgulti/MaxSubGraph_CUDA_C/build /home/fvirgulti/MaxSubGraph_CUDA_C/build /home/fvirgulti/MaxSubGraph_CUDA_C/build/CMakeFiles/testFra.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/testFra.dir/depend

