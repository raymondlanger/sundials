# ---------------------------------------------------------------
# Programmer:  Raymond Langer @ RWTH Aachen University
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# CMakeLists.txt file for the pardiso SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_PARDISO\n\")")

# Add variable sunlinsolpardiso_SOURCES with the sources for the SUNLINSOLPARDISO lib
SET(sunlinsolpardiso_SOURCES sunlinsol_pardiso.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLPARDISO library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  )

# Add variable sunlinsolpardiso_HEADERS with the exported SUNLINSOLPARDISO header files
SET(sunlinsolpardiso_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_pardiso.h
  )

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLPARDISO library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLPARDISO library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolpardiso_static STATIC ${sunlinsolpardiso_HEADERS} ${sunlinsolpardiso_SOURCES} ${shared_SOURCES})
  target_link_libraries(sundials_sunlinsolpardiso_static PUBLIC PARDISO::PARDISO_C)
  target_include_directories(sundials_sunlinsolpardiso_static 
    PUBLIC
    $<BUILD_INTERFACE:${sundials_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
  )
  SET_TARGET_PROPERTIES(sundials_sunlinsolpardiso_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
  add_library(sundials::sundials_sunlinsolpardiso_static ALIAS sundials_sunlinsolpardiso_static)
  INSTALL(TARGETS sundials_sunlinsolpardiso_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLPARDISO library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLPARDISO library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolpardiso_shared SHARED ${sunlinsolpardiso_HEADERS} ${sunlinsolpardiso_SOURCES} ${shared_SOURCES})
  target_link_libraries(sundials_sunlinsolpardiso_shared PUBLIC PARDISO::PARDISO_C sundials_sunmatrixsparse_shared)
  target_include_directories(sundials_sunlinsolpardiso_shared 
    PUBLIC
    $<BUILD_INTERFACE:${sundials_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
  )

  message(STATUS "Building pardiso interface as shared library")

  SET_TARGET_PROPERTIES(sundials_sunlinsolpardiso_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolpardiso_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})

  add_library(sundials::sundials_sunlinsolpardiso_shared ALIAS sundials_sunlinsolpardiso_shared)
  INSTALL(TARGETS sundials_sunlinsolpardiso_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

if(BUILD_SHARED_LIBS)
  add_library(sundials::sundials_sunlinsolpardiso ALIAS sundials_sunlinsolpardiso_shared)
else()
  add_library(sundials::sundials_sunlinsolpardiso ALIAS sundials_sunlinsolpardiso_static)
endif()

# Install the SUNLINSOLPARDISO header files
INSTALL(FILES ${sunlinsolpardiso_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled and building with 64-bit indices, then build and install 
# the FSUNLINSOLPARDISO library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsolpardiso_SOURCES fsunlinsol_pardiso.c)

  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolpardiso_static STATIC ${fsunlinsolpardiso_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolpardiso_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolpardiso_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)

  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolpardiso_shared ${fsunlinsolpardiso_SOURCES})

    # fsunlinsolpardiso depends on fnvecserial, fsunmatrixsparse, sunlinsolpardiso
    TARGET_LINK_LIBRARIES(sundials_fsunlinsolpardiso_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixsparse_shared
      sundials_sunlinsolpardiso_shared)

    SET_TARGET_PROPERTIES(sundials_fsunlinsolpardiso_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolpardiso_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolpardiso_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)

ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_PARDISO module")
