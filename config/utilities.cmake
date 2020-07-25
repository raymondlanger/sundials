function(find_blas_lapack_mkl_preferred MY_VENDOR MY_INTEL_MKL_DIR MY_SEQ MY_FROM_WHERE)
  message(STATUS "${MY_FROM_WHERE} blas/lapack search: BLA_VENDOR=${MY_VENDOR}, sequential=${MY_SEQ}")

  if(MY_INTEL_MKL_DIR)
    if(MY_INTEL_MKL_DIR)
      message(STATUS "${MY_FROM_WHERE} blas/lapack search: mkl search directory=\"${MY_INTEL_MKL_DIR}\"")
    endif()
    if(NOT EXISTS "${MY_INTEL_MKL_DIR}")
      message(
        FATAL_ERROR
          "Could not find the specified intel mkl directory \"${MY_INTEL_MKL_DIR}\"")
    else()
      set(INTEL_MKL_DIR_CACHED "${MY_INTEL_MKL_DIR}"
          CACHE PATH "directory for the mkl"
          FORCE)
      set(ENV{DYLD_LIBRARY_PATH}
          "${INTEL_MKL_DIR_CACHED}:${INTEL_MKL_DIR_CACHED}/lib:ENV{DYLD_LIBRARY_PATH}")
      set(MKL_HEADER_SEARCH_HINTS "${INTEL_MKL_DIR_CACHED}")
      if(EXISTS "${INTEL_MKL_DIR_CACHED}/include")
        set(MKL_HEADER_SEARCH_HINTS "${MKL_HEADER_SEARCH_HINTS};${INTEL_MKL_DIR_CACHED}/include")
      elseif()
        get_filename_component(AUX "${MKL_HEADER_SEARCH_HINTS}/.." ABSOLUTE)
        if(EXISTS "${AUX}/include")
          set(MKL_HEADER_SEARCH_HINTS "${MKL_HEADER_SEARCH_HINTS};${AUX}/include")
        endif()
      endif()
    endif()
  else()
    message(
      STATUS
        "HINT: You can specify an additional directory which is used to find the MKL"
      )
    message(
      STATUS
        "  cmake ../ -DINTEL_MKL_DIR=\"/opt/intel/compilers_and_libraries_2019/mac/mkl/\""
      )
  endif(MY_INTEL_MKL_DIR)

  if(EXISTS "$ENV{MKLROOT}")
    message(STATUS "environment variable MKLROOT is to find the mkl")
    set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:$ENV{MKLROOT}/lib")
    set(MKL_HEADER_SEARCH_HINTS "${MKL_HEADER_SEARCH_HINTS};$ENV{MKLROOT}/include")
  elseif(EXISTS "$ENV{MKL_HOME}")
    message(STATUS "environment variable MKL_HOME is to find the mkl")
    set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:$ENV{MKL_HOME}/lib")
    set(MKL_HEADER_SEARCH_HINTS "${MKL_HEADER_SEARCH_HINTS};$ENV{MKL_HOME}/include")
  endif()

  # check some default locations where the mkl might be found
  if(NOT MY_INTEL_MKL_DIR)
    # checking every year...
    set(MKL_AUX "/opt/intel/compilers_and_libraries_2019/mac/mkl/")
    if(EXISTS "${MKL_AUX}")
      set(MY_INTEL_MKL_DIR "${MKL_AUX}")
    endif()
    # more recent versions are preferred...
    set(MKL_AUX "/opt/intel/compilers_and_libraries_2020/mac/mkl/")
    if(EXISTS "${MKL_AUX}")
      set(MY_INTEL_MKL_DIR "${MKL_AUX}")
    endif()
    if(EXISTS "${MY_INTEL_MKL_DIR}")
      message(STATUS "search for the mkl in \"${MY_INTEL_MKL_DIR}\"")
      set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:${MY_INTEL_MKL_DIR}/lib")
      set(MKL_HEADER_SEARCH_HINTS "${MY_INTEL_MKL_DIR}")
    endif()
  endif()


  if(NOT MY_VENDOR)
    message(STATUS "no blas vendor specified -- attempting to find the mkl")
    if(MY_SEQ)
      set(BLA_VENDOR "Intel10_64lp_seq")
    else()
      set(BLA_VENDOR "Intel10_64lp")
    endif()
    message(STATUS "BLA_VENDOR: \"${BLA_VENDOR}\" (use, e.g., -DBLA_VENDOR=Apple to override the preferred vendor explicitly)")
    option(BLA_STATIC ON)
    find_package(BLAS)
    find_package(LAPACK)
    if(LAPACK_FOUND AND BLAS_FOUND)
      find_path(MKL_INCLUDE_DIR mkl.h HINTS "${MY_INTEL_MKL_DIR}" PATH_SUFFIXES include/)
      set(EIGEN_CAN_USE_MKL ON PARENT_SCOPE)
      set(MKL_INCLUDE_DIR "${MKL_INCLUDE_DIR}" PARENT_SCOPE)
      set(BLAS_LIBRARIES "${BLAS_LIBRARIES}" PARENT_SCOPE)
      set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES}" PARENT_SCOPE)
      return()
    endif()
  endif(NOT MY_VENDOR)
  
  set(BLA_VENDOR "${MY_VENDOR}")
  message(STATUS "Using user supplied BLA_VENDOR=\"${BLA_VENDOR}\"")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  set(BLAS_LIBRARIES "${BLAS_LIBRARIES}" PARENT_SCOPE)
  set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES}" PARENT_SCOPE)
  if(MY_VENDOR MATCHES Intel)
    find_path(MKL_INCLUDE_DIR mkl.h HINTS "${MY_INTEL_MKL_DIR}" PATH_SUFFIXES include/)
    set(MKL_INCLUDE_DIR "${MKL_INCLUDE_DIR}" PARENT_SCOPE)
    set(EIGEN_CAN_USE_MKL ON PARENT_SCOPE)
  else()
    set(EIGEN_CAN_USE_MKL OFF PARENT_SCOPE)
  endif()
endfunction(find_blas_lapack_mkl_preferred)
