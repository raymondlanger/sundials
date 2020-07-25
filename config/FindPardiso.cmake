# ---------------------------------------------------------------
# Programmer:  Raymond Langer @ RWTH Aachen Univeristy
# ---------------------------------------------------------------

# - Try to find pardiso 
# You can use the following configuration command to provide the
# the pardiso library location (you would specify the absolute path):
# cmake .. -DPARDISO_LIBRARY_DIR=$HOME/lib 
#
# Once done this will define
#  PARDISO_FOUND - True means system has pardiso
#  PARDISO_LIBRARIES - The libraries needed for using pardiso
#  PARDISO::PARDISO_CXX - The imported target (it's strongly recommended 
#                     to use this imported target and nothing else)
#  PARDISO::PARDISO_C - The imported target (it's strongly recommended 
#                     to use this imported target and nothing else)
# This is only provided for the "normal" Pardiso (not MKL_PARDISO)
#  PARDISO_NAME - The library file name pardiso
#  PARDISO_DIR - The directory of the pardiso library pardiso

if(NOT EXISTS ${PARDISO_LIBRARY_DIR_CACHED})
  message(STATUS "Pardiso directory does not exist \"${PARDISO_LIBRARY_DIR_CACHED}\"")
endif()

if(NOT MKL_PARDISO)
  if(PARDISO_GCC720)
  set(pardiso_names 
    pardiso600-GNU720-X86-64
    )
  else()
  set(pardiso_names 
    pardiso500-INTEL1301-X86-64
    pardiso500-GNU461-X86-64
    pardiso500-GNU472-X86-64
    pardiso500-GNU481-X86-64
    pardiso600-INTEL1301-X86-64
    pardiso600-GNU461-X86-64
    pardiso600-GNU472-X86-64
    pardiso600-GNU481-X86-64
    pardiso500-WIN-X86-64
    pardiso500-MACOS-X86-64
    pardiso600-WIN-X86-64
    pardiso600-MACOS-X86-64
    pardiso600-GNU720-X86-64
    pardiso600-GNU800-X86-64
    )
  endif(PARDISO_GCC720)
  
  find_library(PARDISO_LIBRARY NAMES ${pardiso_names} HINTS "${PARDISO_LIBRARY_DIR_CACHED}")
  find_package(GFortranLibs REQUIRED)
  find_package(OpenMP REQUIRED)
  find_package(Threads REQUIRED)
  set(pardiso_link_libraries "${LIBGFORTRAN_LIBRARIES}")
else(NOT MKL_PARDISO)
  set(pardiso_link_libraries "")
endif(NOT MKL_PARDISO)

if(ORK_INTEGRATION)
  include(utilities)
  find_blas_lapack_mkl_preferred("${BLA_VENDOR}" "${INTEL_MKL_DIR}" "${SEQ_LAPACK}" "sundials_pardiso")
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif(ORK_INTEGRATION)

list(APPEND pardiso_link_libraries "${BLAS_LIBRARIES}")
list(APPEND pardiso_link_libraries "${LAPACK_LIBRARIES}")

if(NOT MKL_PARDISO)
  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set PARDISO_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(Pardiso  DEFAULT_MSG
                                    PARDISO_LIBRARY LIBGFORTRAN_LIBRARIES)
  
  get_filename_component(PARDISO_DIR ${PARDISO_LIBRARY} DIRECTORY)
  message(STATUS "pardiso directory: ${PARDISO_DIR}")
  get_filename_component(PARDISO_NAME ${PARDISO_LIBRARY} NAME)
  message(STATUS "pardiso name: ${PARDISO_NAME}")
  
  mark_as_advanced( PARDISO_DIR PARDISO_NAME LIBGFORTRAN_LIBRARIES )
else()
  message(STATUS "Skip search for the pardiso library")
  find_package_handle_standard_args(Pardiso  DEFAULT_MSG)
endif(NOT MKL_PARDISO)

mark_as_advanced( PARDISO_LIBRARY )
set(PARDISO_LIBRARIES ${pardiso_link_libraries})

if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

if(NOT MKL_PARDISO)
add_library(PARDISO::PARDISO_CXX SHARED IMPORTED)
add_library(PARDISO::PARDISO_C SHARED IMPORTED)
else()
add_library(PARDISO::PARDISO_CXX INTERFACE IMPORTED)
add_library(PARDISO::PARDISO_C INTERFACE IMPORTED)
endif()

if(INTEL_MKL_DIR OR MKL_INCLUDE_DIR)
  if(EXISTS "${MKL_INCLUDE_DIR}")
    set_property(TARGET 
      PARDISO::PARDISO_CXX PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
    set_property(TARGET 
      PARDISO::PARDISO_C PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
endif()

if(NOT MKL_PARDISO)
set_property(TARGET PARDISO::PARDISO_CXX PROPERTY IMPORTED_NO_SONAME 1)
set_property(TARGET PARDISO::PARDISO_C PROPERTY IMPORTED_NO_SONAME 1)
set_target_properties(PARDISO::PARDISO_CXX PROPERTIES
  IMPORTED_LOCATION 
    "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}")

if(APPLE)
  set_property(TARGET PARDISO::PARDISO_CXX PROPERTY
    INTERFACE_LINK_LIBRARIES
       ${pardiso_link_libraries} ${OpenMP_CXX_FLAGS} ${OpenMP_libomp_LIBRARY} Threads::Threads)
else()
  set_property(TARGET PARDISO::PARDISO_CXX PROPERTY
    INTERFACE_LINK_LIBRARIES
       ${pardiso_link_libraries} ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()
set(OpenMP_CXX_FLAGS_CPY ${OpenMP_CXX_FLAGS})
string (REPLACE " " ";" OpenMP_CXX_FLAGS_CPY "${OpenMP_CXX_FLAGS_CPY}")
set_property(TARGET PARDISO::PARDISO_CXX
             PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS_CPY})

set_target_properties(PARDISO::PARDISO_C PROPERTIES
  IMPORTED_LOCATION 
    "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}")
if(APPLE)
  set_property(TARGET PARDISO::PARDISO_C PROPERTY
    INTERFACE_LINK_LIBRARIES
       ${pardiso_link_libraries} ${OpenMP_C_FLAGS} ${OpenMP_libomp_LIBRARY} Threads::Threads)
else()
  set_property(TARGET PARDISO::PARDISO_C PROPERTY
    INTERFACE_LINK_LIBRARIES
       ${pardiso_link_libraries} ${OpenMP_C_FLAGS} Threads::Threads)
endif()
set(OpenMP_C_FLAGS_CPY ${OpenMP_C_FLAGS})
string (REPLACE " " ";" OpenMP_C_FLAGS_CPY "${OpenMP_C_FLAGS_CPY}")
set_property(TARGET PARDISO::PARDISO_C
             PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_C_FLAGS_CPY})
else()
  set_property(TARGET PARDISO::PARDISO_CXX PROPERTY
    INTERFACE_LINK_LIBRARIES
       ${pardiso_link_libraries})
  set_property(TARGET PARDISO::PARDISO_C PROPERTY
    INTERFACE_LINK_LIBRARIES
       ${pardiso_link_libraries})
endif(NOT MKL_PARDISO)
