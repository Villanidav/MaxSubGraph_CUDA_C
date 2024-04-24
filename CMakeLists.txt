cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(subgraph_search VERSION 5.0)
enable_language(CXX)

# ########################################################################
# ### Tweak the default configuration
# ########################################################################

# define an option for using native cpp code
set(TARGET_CPU_ARCHITECTURE "generic" CACHE STRING "The target cpu for architecture-specific optimizations")
set_property(CACHE TARGET_CPU_ARCHITECTURE PROPERTY STRINGS generic native)

# the idea is that we want to build the optimized version, unless the user says otherwise
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The type of build" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# enable additional warning flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(-Wall -Wextra -Wpedantic -Wl,-z,defs -Wl,-z,now -Wl,-z,relro)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion)
endif()

# enable colored compiler output (with Ninja)
if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-fdiagnostics-color=always)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-fcolor-diagnostics)
  endif()
endif()



# ########################################################################
# ### Add external dependency
# ########################################################################

# ########################################################################
# ### Define the application sources
# ########################################################################

# application headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND header_files
  "${header_path}/main.hpp"
  "${header_path}/test.hpp"
  "${header_path}/testFra.hpp"
)


# application sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_filess
  "${source_path}/Label.cpp"
  "${source_path}/gen_rotations.cpp"
)

# ########################################################################
# ### Define the compilation step
# ########################################################################

# define the target and its properties
add_executable(subgraph_search ${header_files} "${source_path}/main.cpp" ${source_files}
        src/mc_split.cpp
        src/mol_mcs.cpp
        src/g2mol.cpp
)
set_target_properties(subgraph_search
  PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)


add_executable(test ${header_files} "${source_path}/test.cpp"
        "${source_path}/gen_bond_labels.cpp"
        "${source_path}/gen_rotations.cpp"
        "${source_path}/select_vertex.cpp"
        "${source_path}/hood.cpp"
        "${source_path}/Label.cpp"
        "${source_path}/smiles_mcs.cpp"
        "${source_path}/calc_bound.cpp"
        "${source_path}/select_label.cpp"
        "${source_path}/gen_initial_labels.cpp"
        "${source_path}/gen_ring_classes.cpp"
        "${source_path}/gen_rotations.cpp"
        "${source_path}/search_mcs.cpp"
        "${source_path}/mc_split.cpp.cpp"
        "${source_path}/g2mol.cpp.cpp"

)
set_target_properties(test PROPERTIES
    COMPILE_DEFINITIONS "TESTING"
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

add_executable(testFra ${header_files} "${source_path}/testFra.cpp"
        "${source_path}/gen_bond_labels.cpp"
        "${source_path}/gen_rotations.cpp"
        "${source_path}/select_vertex.cpp"
        "${source_path}/hood.cpp"
        "${source_path}/Label.cpp"
        "${source_path}/smiles_mcs.cpp"
        "${source_path}/calc_bound.cpp"
        "${source_path}/select_label.cpp"
        "${source_path}/gen_initial_labels.cpp"
        "${source_path}/gen_ring_classes.cpp"
        "${source_path}/gen_rotations.cpp"
        "${source_path}/search_mcs.cpp"
        src/mc_split.cpp
        src/g2mol.cpp
)

set_target_properties(test PROPERTIES
    COMPILE_DEFINITIONS "TESTING"
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)




# enable link-time optimizations
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)

if(ipo_supported)
  set_property(TARGET subgraph_search PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# since RDKit doesn't handle the dependencies in a correct way, we need to improvise
# NOTE: by using this way of importing RDKit, we need to manually include its dependencies

# the rdkit to perform the heavy lifting
find_package(RDKit REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})




get_target_property(RDKIT_LIB_FULLPATH RDKit::RDGeneral LOCATION)
get_filename_component(RDKIT_LIB_DIRPATH "${RDKIT_LIB_FULLPATH}" DIRECTORY)
cmake_path(GET RDKIT_LIB_DIRPATH PARENT_PATH RDKIT_INSTALL_PREFIX)
include_directories(${RDKit_INCLUDE_DIRS}) #modified
include_directories(subgraph_search PUBLIC "${RDKIT_INSTALL_PREFIX}/include/rdkit")
target_link_directories(subgraph_search PUBLIC "${RDKIT_INSTALL_PREFIX}/lib")
target_link_directories(subgraph_search PUBLIC "${RDKIT_INSTALL_PREFIX}/lib64")
target_link_libraries(subgraph_search PUBLIC
        -lRDKit
        -lRDKitFileParsers
        -lRDKitGraphMol
        -lRDKitRDGeneral
        -lRDKitSmilesParse
        -lRDKitSubstructMatch
        -lRDKitChem
)
include_directories(test PUBLIC "${RDKIT_INSTALL_PREFIX}/include/rdkit")
target_link_directories(test PUBLIC "${RDKIT_INSTALL_PREFIX}/lib")
target_link_directories(test PUBLIC "${RDKIT_INSTALL_PREFIX}/lib64")
target_link_libraries(test PUBLIC
        -lRDKit
        -lRDKitFileParsers
        -lRDKitGraphMol
        -lRDKitRDGeneral
        -lRDKitSmilesParse
        -lRDKitSubstructMatch
        -lRDKitChem
)

#modificato togliendo target_ prima di link in ogni punto

include_directories(testFra PUBLIC "${RDKIT_INSTALL_PREFIX}/include/rdkit")
target_link_directories(testFra PUBLIC "${RDKIT_INSTALL_PREFIX}/lib")
target_link_directories(testFra PUBLIC "${RDKIT_INSTALL_PREFIX}/lib64")
target_link_libraries(testFra PUBLIC
        -lRDKit
        -lRDKitFileParsers
        -lRDKitGraphMol
        -lRDKitRDGeneral
        -lRDKitSmilesParse
        -lRDKitSubstructMatch
        -lRDKitChem
)

message(STATUS "RDKit Include Directories: ${RDKit_INCLUDE_DIRS}")
message(STATUS "RDKit Library Directories: ${RDKit_LIBRARY_DIRS}")
