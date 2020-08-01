# -----------------------------------------------------------------
# Programmer:  Raymond Langer @ RWTH Aachen University
# -----------------------------------------------------------------
# ---------------------------------------------------------------
# CMakeLists.txt file for sunlinsol UMFPACK examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is
# 'develop' for examples excluded from 'make test' in releases

# Examples using the UMFPACK linear solver
set(sunlinsol_umfpack_examples
  "test_sunlinsol_umfpack\;300 0 0\;"
  "test_sunlinsol_umfpack\;300 1 0\;"
  "test_sunlinsol_umfpack\;1000 0 0\;"
  "test_sunlinsol_umfpack\;1000 1 0\;"
  )

# Dependencies for nvector examples
set(sunlinsol_umfpack_dependencies
  test_sunlinsol
  sundials_nvector
  sundials_matrix
  sundials_linearsolver
  )

if(POLICY CMP0054)
   cmake_policy(SET CMP0054 NEW)
endif()

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against (through the target that was used to 
# generate them) based on the value of the variable LINK_LIBRARY_TYPE
if(LINK_LIBRARY_TYPE MATCHES "static")
  set(NVEC_LIB sundials_nvecserial_static)
  set(SUNMAT_LIB sundials_sunmatrixdense_static)
  list(APPEND SUNMAT_LIB sundials_sunmatrixsparse_static)
  set(SUNLINSOL_LIB sundials_sunlinsolumfpack_static)
ELSE()
  set(NVEC_LIB sundials_nvecserial_shared)
  set(SUNMAT_LIB sundials_sunmatrixdense_shared)
  list(APPEND SUNMAT_LIB sundials_sunmatrixsparse_shared)
  set(SUNLINSOL_LIB sundials_sunlinsolumfpack_shared)
endif()

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVEC_LIB} ${SUNMAT_LIB} ${SUNLINSOL_LIB} ${EXTRA_LINK_LIBS})


# Add the build and install targets for each example
FOREACH(example_tuple ${sunlinsol_umfpack_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # This is used to get around DLL linkage issue since we are
  # manually including sundials_nvector.c here, which is normally in
  # a library that is included.  If this is not set build system
  # thinks nvector is externally linked.
  if(WIN32)
    add_definitions(-DBUILD_SUNDIALS_LIBRARY)
  endif(WIN32)

  # check if this example has already been added, only need to add
  # example source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files

    add_executable(${example} ${example}.c
      ../test_sunlinsol.c
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against    
    target_link_libraries(${example} ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  ELSE()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source files
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.c
      ../test_sunlinsol.h
      ../test_sunlinsol.c 
      ../../../src/sundials/sundials_matrix.c
      ../../../src/sundials/sundials_linearsolver.c
      ../../../src/sundials/sundials_nvector.c
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/umfpack)
  endif()

ENDFOREACH(example_tuple ${sunlinsol_umfpack_examples})

if(EXAMPLES_INSTALL)

  # Install the README file
  install(FILES DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/umfpack)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER_LIB "sundials_sunlinsolumfpack")
  set(LIBS "${LIBS} -lsundials_sunmatrixsparse -lsundials_sunmatrixdense")

  # Set the link directory for the sparse and dense sunmatrix libraries
  # The generated CMakeLists.txt does not use find_library() locate them
  set(EXTRA_LIBS_DIR "${libdir}")

  EXAMPLES2STRING(sunlinsol_umfpack_examples EXAMPLES_UMFPACK)
  EXAMPLES2STRING(sunlinsol_umfpack_dependencies EXAMPLES_DEPENDENCIES_UMFPACK)

  # Regardless of the platform we're on, we will generate and install 
  # CMakeLists.txt file for building the examples. This file  can then 
  # be used as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/sunlinsol/umfpack/CMakeLists.txt
    @ONLY
    )

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/umfpack/CMakeLists.txt
    DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/umfpack
    )

  # On UNIX-type platforms, we also  generate and install a makefile for 
  # building the examples. This makefile can then be used as a template 
  # for the user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/sunlinsol/umfpack/Makefile_ex
      @ONLY
      )
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/sunlinsol/umfpack/Makefile_ex 
      DESTINATION ${EXAMPLES_INSTALL_PATH}/sunlinsol/umfpack
      RENAME Makefile
      )
  endif(UNIX)

endif(EXAMPLES_INSTALL)
