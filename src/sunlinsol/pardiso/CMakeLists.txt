# ---------------------------------------------------------------
# Programmer:  Raymond Langer @ RWTH Aachen University
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# CMakeLists.txt file for the pardiso SUNLinearSolver library

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_PARDISO\n\")")

# Add variable sunlinsolpardiso_SOURCES with the sources for the SUNLINSOLPARDISO lib
set(sunlinsolpardiso_SOURCES sunlinsol_pardiso.c)

# Add variable sunlinsolpardiso_HEADERS with the exported SUNLINSOLPARDISO header files
set(sunlinsolpardiso_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_pardiso.h
  )

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLPARDISO library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLPARDISO library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunlinsolpardiso_static 
    STATIC 
    ${sunlinsolpardiso_HEADERS} 
    ${sunlinsolpardiso_SOURCES} 
    $<TARGET_OBJECTS:sundials_generic_static_obj>)
  target_link_libraries(sundials_sunlinsolpardiso_static 
    PUBLIC PARDISO::PARDISO_C)
  target_include_directories(sundials_sunlinsolpardiso_static 
    PUBLIC
    $<BUILD_INTERFACE:${sundials_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
  )
  set_target_properties(sundials_sunlinsolpardiso_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
  add_library(sundials::sundials_sunlinsolpardiso_static 
    ALIAS sundials_sunlinsolpardiso_static)
  install(TARGETS sundials_sunlinsolpardiso_static DESTINATION lib)
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLPARDISO library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLPARDISO library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunlinsolpardiso_shared 
    SHARED 
    ${sunlinsolpardiso_HEADERS} 
    ${sunlinsolpardiso_SOURCES} 
    $<TARGET_OBJECTS:sundials_generic_shared_obj>)
  target_link_libraries(sundials_sunlinsolpardiso_shared 
    PUBLIC PARDISO::PARDISO_C sundials_sunmatrixsparse_shared)
  target_include_directories(sundials_sunlinsolpardiso_shared 
    PUBLIC
    $<BUILD_INTERFACE:${sundials_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
  )

  set_target_properties(sundials_sunlinsolpardiso_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_sunlinsolpardiso_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})

  add_library(sundials::sundials_sunlinsolpardiso_shared 
    ALIAS sundials_sunlinsolpardiso_shared)
  install(TARGETS sundials_sunlinsolpardiso_shared DESTINATION lib)
endif(BUILD_SHARED_LIBS)

if(BUILD_SHARED_LIBS)
  add_library(sundials::sundials_sunlinsolpardiso 
    ALIAS sundials_sunlinsolpardiso_shared)
else()
  add_library(sundials::sundials_sunlinsolpardiso 
    ALIAS sundials_sunlinsolpardiso_static)
endif()

# Install the SUNLINSOLPARDISO header files
install(FILES ${sunlinsolpardiso_HEADERS} DESTINATION include/sunlinsol)


##########
#       WARNING: The FORTRAN interface is INCOMPLETE AND DOES NOT WORK.
##########
# If FCMIX is enabled and building with 64-bit indices, then build and install 
# the FSUNLINSOLPARDISO library
if(FCMIX_ENABLE AND F77_FOUND)
  set(fsunlinsolpardiso_SOURCES fsunlinsol_pardiso.c)

  if(BUILD_STATIC_LIBS)
    add_library(sundials_fsunlinsolpardiso_static STATIC ${fsunlinsolpardiso_SOURCES})
    set_target_properties(sundials_fsunlinsolpardiso_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
    install(TARGETS sundials_fsunlinsolpardiso_static DESTINATION lib)
  endif(BUILD_STATIC_LIBS)

  if(BUILD_SHARED_LIBS)
    add_library(sundials_fsunlinsolpardiso_shared ${fsunlinsolpardiso_SOURCES})

    # fsunlinsolpardiso depends on fnvecserial, fsunmatrixsparse, sunlinsolpardiso
    target_link_libraries(sundials_fsunlinsolpardiso_shared
      sundials_fnvecserial_shared
      sundials_fsunmatrixsparse_shared
      sundials_sunlinsolpardiso_shared)

    set_target_properties(sundials_fsunlinsolpardiso_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolpardiso CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(sundials_fsunlinsolpardiso_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    install(TARGETS sundials_fsunlinsolpardiso_shared DESTINATION lib)
  endif(BUILD_SHARED_LIBS)

endif(FCMIX_ENABLE AND F77_FOUND)

#
message(STATUS "Added SUNLINSOL_PARDISO module")
