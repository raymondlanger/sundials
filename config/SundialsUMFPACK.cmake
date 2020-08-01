# -----------------------------------------------------------------
# Programmer:  Raymond Langer @ RWTH Aachen University,
#                               Institute for Combustion Technology
# -----------------------------------------------------------------
# ---------------------------------------------------------------

### This is only set if running GUI - simply return first time enabled
if(UMFPACK_DISABLED)
  set(UMFPACK_DISABLED FALSE CACHE INTERNAL "GUI - UMFPACK now enabled" FORCE)
  return()
endif()

set(UMFPACK_FOUND FALSE)
find_package(UMFPACK REQUIRED)

# Create the UMFPACKTest directory
set(UMFPACKTest_DIR ${PROJECT_BINARY_DIR}/UMFPACKTest)
file(MAKE_DIRECTORY ${UMFPACKTest_DIR})
# Create a CMakeLists.txt file
file(WRITE ${UMFPACKTest_DIR}/CMakeLists.txt
  "CMAKE_MINIMUM_REQUIRED(VERSION 3.9)\n"
  "PROJECT(ltest C)\n"
  "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
  "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
  "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
  "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
  "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
  "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
  "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
  "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
  "SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)\n"
  "if (POLICY CMP0012)\n"
  "  cmake_policy(SET CMP0012 NEW)\n"
  "endif ()\n"
  "if(POLICY CMP0054)\n"
  "  cmake_policy(SET CMP0054 NEW)\n"
  "endif()\n"
  "\n"
  "if(COMMAND cmake_policy)\n"
  "  cmake_policy(SET CMP0003 NEW)\n"
  "endif(COMMAND cmake_policy)\n"
  "\n"
  "find_package(UMFPACK REQUIRED)\n"
  "\n"
  "ADD_EXECUTABLE(ltest ltest.c)\n"
  "TARGET_LINK_LIBRARIES(ltest SUITE_SPARSE::UMFPACK)\n")

file(WRITE ${UMFPACKTest_DIR}/ltest.c
  "\#include <stdio.h>\n"
  "\#include \"umfpack.h\"\n"
    "\n"
  "int n = 5 ;\n"
  "int Ap [ ] = {0, 2, 5, 9, 10, 12} ;\n"
  "int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;\n"
  "double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;\n"
  "double b [ ] = {8., 45., -3., 3., 19.} ;\n"
  "double x [5] ;\n"
  "int main (void)\n"
  "{\n"
  "double *null = (double *) NULL ;\n"
  "int i ;\n"
  "void *Symbolic, *Numeric ;\n"
  "(void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;\n"
  "(void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;\n"
  "umfpack_di_free_symbolic (&Symbolic) ;\n"
  "(void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;\n"
  "umfpack_di_free_numeric (&Numeric) ;\n"
  "for (i = 0 ; i < n ; i++) printf (\"x [%d] = %g\\n\", i, x [i]) ;\n"
  "return (0) ;\n"
  "}\n")


# Attempt to link the "ltest" executable
try_compile(LTEST_OK ${UMFPACKTest_DIR} ${UMFPACKTest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)

# To ensure we do not use stuff from the previous attempts,
# we must remove the CMakeFiles directory.
file(REMOVE_RECURSE ${UMFPACKTest_DIR}/CMakeFiles)
# Process test result
if(LTEST_OK)
  message(STATUS "Checking if UMFPACK works... OK")
  set(UMFPACK_FOUND TRUE)
else(LTEST_OK)
  message(STATUS "Checking if UMFPACK works... FAILED")
endif(LTEST_OK)
