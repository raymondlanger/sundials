# ---------------------------------------------------------------
# Programmer:  Raymond Langer @ RWTH Aachen University
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# CMakeLists.txt file for the umfpack SUNLinearSolver library

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_UMFPACK\n\")")

# Add variable sunlinsolumfpack_SOURCES with the sources for the SUNLINSOLUMFPACK lib
set(sunlinsolumfpack_SOURCES sunlinsol_umfpack.c)

# Add variable sunlinsolumfpack_HEADERS with the exported SUNLINSOLUMFPACK header files
set(sunlinsolumfpack_HEADERS
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_umfpack.h
  )

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLUMFPACK library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLUMFPACK library
if(BUILD_STATIC_LIBS)
  add_library(sundials_sunlinsolumfpack_static 
    STATIC 
    ${sunlinsolumfpack_HEADERS} 
    ${sunlinsolumfpack_SOURCES} 
    $<TARGET_OBJECTS:sundials_generic_static_obj>)
  target_link_libraries(sundials_sunlinsolumfpack_static 
    PUBLIC SUITE_SPARSE::UMFPACK)
  target_include_directories(sundials_sunlinsolumfpack_static 
    PUBLIC
    $<BUILD_INTERFACE:${sundials_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
  )
  set_target_properties(sundials_sunlinsolumfpack_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolumfpack CLEAN_DIRECT_OUTPUT 1)
  add_library(sundials::sundials_sunlinsolumfpack_static 
    ALIAS sundials_sunlinsolumfpack_static)
  install(TARGETS sundials_sunlinsolumfpack_static DESTINATION lib)
endif(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLUMFPACK library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLUMFPACK library
if(BUILD_SHARED_LIBS)
  add_library(sundials_sunlinsolumfpack_shared 
    SHARED 
    ${sunlinsolumfpack_HEADERS} 
    ${sunlinsolumfpack_SOURCES} 
    $<TARGET_OBJECTS:sundials_generic_shared_obj>)
  target_link_libraries(sundials_sunlinsolumfpack_shared 
    PUBLIC SUITE_SPARSE::UMFPACK sundials_sunmatrixsparse_shared)
  target_include_directories(sundials_sunlinsolumfpack_shared 
    PUBLIC
    $<BUILD_INTERFACE:${sundials_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
  )
  set_target_properties(sundials_sunlinsolumfpack_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolumfpack CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(sundials_sunlinsolumfpack_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  add_library(sundials::sundials_sunlinsolumfpack_shared 
    ALIAS sundials_sunlinsolumfpack_shared)
  install(TARGETS sundials_sunlinsolumfpack_shared DESTINATION lib)
endif(BUILD_SHARED_LIBS)

if(BUILD_SHARED_LIBS)
  add_library(sundials::sundials_sunlinsolumfpack 
    ALIAS sundials_sunlinsolumfpack_shared)
else()
  add_library(sundials::sundials_sunlinsolumfpack 
    ALIAS sundials_sunlinsolumfpack_static)
endif()

# Install the SUNLINSOLUMFPACK header files
install(FILES ${sunlinsolumfpack_HEADERS} DESTINATION include/sunlinsol)

message(STATUS "Added SUNLINSOL_UMFPACK module")
