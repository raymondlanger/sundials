# Mostly copied from
# https://cmake-basis.github.io/apidoc/latest/FindUMFPACK_8cmake_source.html
 
if(UMFPACK_INCLUDE_DIR AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif()

if(DEFINED UMFPACKDIR)
  if(NOT IS_DIRECTORY "${UMFPACKDIR}")
    message(FATAL_ERROR "UMFPACK directory does not exist: \"${UMFPACKDIR}\"")
  else()
    message(STATUS "Preferred umfpack location: \"${UMFPACKDIR}\"")
  endif()
  if(NOT IS_DIRECTORY "${UMFPACKDIR}/include")
    message(FATAL_ERROR "UMFPACK directory does not exist: \"${UMFPACKDIR}/include\"")
  endif()
  if(NOT IS_DIRECTORY "${UMFPACKDIR}/lib")
    message(FATAL_ERROR "UMFPACK directory does not exist: \"${UMFPACKDIR}/lib\"")
  endif()
else()
  message(STATUS "Use -DUMFPACKDIR=/Some/Dir to specify an umfpack location")
endif()

 find_path(UMFPACK_INCLUDE_DIR
   NAMES
     umfpack.h
   HINTS
     ${UMFPACKDIR}/include
     $ENV{UMFPACKDIR}/include
     ${INCLUDE_INSTALL_DIR}
   REQUIRED
   PATH_SUFFIXES
     suitesparse
     ufsparse
 )
 
 find_library(UMFPACK_LIBRARIES
   NAMES umfpack
   REQUIRED
   HINTS
     ${UMFPACKDIR}/lib
     $ENV{UMFPACKDIR}/lib
     ${LIB_INSTALL_DIR}
 )
 
 if (UMFPACK_LIBRARIES)
 
   if (NOT UMFPACK_LIBDIR)
     get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARIES} PATH)
   endif(NOT UMFPACK_LIBDIR)
   message(STATUS "umfpack directory: ${UMFPACK_LIBDIR}")
   get_filename_component(UMFPACK_NAME ${UMFPACK_LIBRARIES} NAME)
   message(STATUS "umfpack name: ${UMFPACK_NAME}")
 
   find_library(COLAMD_LIBRARY colamd HINTS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR}/lib ${LIB_INSTALL_DIR})
   if (COLAMD_LIBRARY)
     set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${COLAMD_LIBRARY})
   endif (COLAMD_LIBRARY)
   
   find_library(AMD_LIBRARY amd HINTS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR}/lib ${LIB_INSTALL_DIR})
   if (AMD_LIBRARY)
     set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${AMD_LIBRARY})
   endif (AMD_LIBRARY)
 
   find_library(SUITESPARSE_LIBRARY SuiteSparse HINTS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR}/lib ${LIB_INSTALL_DIR})
   if (SUITESPARSE_LIBRARY)
     set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${SUITESPARSE_LIBRARY})
   endif (SUITESPARSE_LIBRARY)
 
 endif ()
 
 include(FindPackageHandleStandardArgs)
 find_package_handle_standard_args(UMFPACK DEFAULT_MSG
                                   UMFPACK_INCLUDE_DIR UMFPACK_LIBRARIES)
 
 if (UMFPACK_FOUND)
   string (REGEX REPLACE "/(suitesparse|ufsparse)/?" "" UMFPACK_DIR "${UMFPACK_INCLUDE_DIR}")
   get_filename_component(UMFPACK_DIR "${UMFPACK_DIR}" PATH)
 endif ()
 
 mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARIES AMD_LIBRARY COLAMD_LIBRARY SUITESPARSE_LIBRARY)


if(ORK_INTEGRATION)
  include(utilities)
  find_blas_lapack_mkl_preferred("${BLA_VENDOR}" "${INTEL_MKL_DIR}" "${SEQ_LAPACK}" "sundials_umfpack")
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif(ORK_INTEGRATION)

set(umfpack_link_libraries "")
list(APPEND umfpack_link_libraries "${BLAS_LIBRARIES}")
list(APPEND umfpack_link_libraries "${LAPACK_LIBRARIES}")
list(APPEND umfpack_link_libraries "${UMFPACK_LIBRARIES}")

add_library(SUITE_SPARSE::UMFPACK UNKNOWN IMPORTED)

set_target_properties(SUITE_SPARSE::UMFPACK PROPERTIES
  IMPORTED_LOCATION 
    "${UMFPACK_LIBDIR}/${UMFPACK_NAME}")

set_property(TARGET 
  SUITE_SPARSE::UMFPACK PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIR}")

set_property(TARGET SUITE_SPARSE::UMFPACK PROPERTY
  INTERFACE_LINK_LIBRARIES
     ${umfpack_link_libraries})


