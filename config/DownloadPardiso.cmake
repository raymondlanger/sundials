
# update this __first part__ of the complete link, if there is a new version of pardiso
set(DEFAULT_DOWNLOAD_URL "https://pardiso-project.org/download/ajan2018600bb6ex6hze7cv15")
# update this if the pardiso website changes
set(EXPECTED_URL_PART "https://pardiso-project.org/download/")

option(PARDISO_GCC720 "Download a pardiso version that links GLIBC_2.27" OFF)

if(PARDISO_URL)
  string(REGEX REPLACE "\n" "" PARDISO_URL "${PARDISO_URL}")
  string(REGEX REPLACE "\t" "" PARDISO_URL "${PARDISO_URL}")
  string(REGEX REPLACE "\r" "" PARDISO_URL "${PARDISO_URL}")
  string(REGEX REPLACE " "  "" PARDISO_URL "${PARDISO_URL}")
  message(STATUS "PARDISO URL: '${PARDISO_URL}'")
  string(REGEX MATCH "^https://pardiso-project.org/download/.+" EXPECTED_URL "${PARDISO_URL}")
  if(NOT EXPECTED_URL)
  	message(WARNING "PARDISO URL is expected to start with ${EXPECTED_URL_PART}...")
  endif()
  set(PARDISO_URL_CACHED "${PARDISO_URL}"
	  CACHE STRING "URL for pardiso download" FORCE)
  set(PARDISO_LIB_NAME "${PARDISO_URL}")
  string(REGEX REPLACE  "^${EXPECTED_URL_PART}.+/(libpardiso.+)" "\\1" PARDISO_LIB_NAME "${PARDISO_LIB_NAME}")
else()
	if(APPLE)
	  set(PARDISO_LIB_NAME libpardiso600-MACOS-X86-64.dylib)
	elseif (UNIX)
          if(PARDISO_GCC720)
	    set(PARDISO_LIB_NAME libpardiso600-GNU720-X86-64.so)
          else()
            set(PARDISO_LIB_NAME libpardiso600-GNU800-X86-64.so)
          endif()
	elseif(WIN32)
	  set(PARDISO_LIB_NAME libpardiso600-WIN-X86-64.dll)
	else()
		message(FATAL_ERROR "Unknown operating system. This information is required to download the correct version of pardiso")
	endif()
  set(PARDISO_URL_CACHED "${DEFAULT_DOWNLOAD_URL}/${PARDISO_LIB_NAME}"
	  CACHE STRING "URL for pardiso download" FORCE)
endif()






set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")

if(PARDISO_LIC_KEY)
  if(NOT PARDISO_URL)
  	message(STATUS "HINT: You can specify the pardiso library URL:")
  	message(STATUS "  cmake ../ -DPARDISO_URL=${PARDISO_URL_CACHED}")
  	message(STATUS "Using the default URL '${PARDISO_URL_CACHED}'")
  endif()

  message(STATUS "Downloading pardiso...")
  if(PARDISO_LIBRARY_DIR)
     set(PARDISO_LIBRARY_DIR_CACHED "${PARDISO_LIBRARY_DIR}"
	     CACHE STRING "pardiso directory" FORCE)
  else()
    set(PARDISO_LIBRARY_DIR_CACHED "${CMAKE_SOURCE_DIR}/generated_pardiso"
	  CACHE STRING "pardiso directory" FORCE)
  	message(STATUS "HINT: You can specify the directory used to store the pardiso library:")
  	message(STATUS "  cmake ../ -DPARDISO_LIBRARY_DIR=$HOME/Software/lib/pardiso6")
  	message(STATUS "Using the default directory '${PARDISO_LIBRARY_DIR}'")
  endif()

  file(MAKE_DIRECTORY ${PARDISO_LIBRARY_DIR_CACHED})
  if(NOT EXISTS "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}")
    file(DOWNLOAD 
     ${PARDISO_URL_CACHED}
     "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}")
    if(APPLE)
      execute_process(
        COMMAND 
        "${CMAKE_INSTALL_NAME_TOOL}" 
        -id 
        "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}" 
        "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}"
        RESULT_VARIABLE
        RES_RUN_INSTALL_NAME
        )
      if(NOT RES_RUN_INSTALL_NAME EQUAL 0)
        message(FATAL_ERROR 
          "install_name_tool could not modify the id of the the downloaded pardiso library\n"
          "FAILED COMMAND:\n"
          "\"${CMAKE_INSTALL_NAME_TOOL}\"" 
          "-id "
          "\"${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}\"" 
          "\"${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}\"\n")
      endif()
    endif()
  endif()

  if(NOT EXISTS "${PARDISO_LIBRARY_DIR_CACHED}/${PARDISO_LIB_NAME}")
    message(FATAL_ERROR "Pardiso library '${PARDISO_LIBRARY_DIR_CACHED}' does not exist.")
  endif()

  message(STATUS "Creating pardiso.lic...")
  string(REGEX REPLACE "\n" "" PARDISO_LIC_KEY "${PARDISO_LIC_KEY}")
  string(REGEX REPLACE "\t" "" PARDISO_LIC_KEY "${PARDISO_LIC_KEY}")
  string(REGEX REPLACE "\r" "" PARDISO_LIC_KEY "${PARDISO_LIC_KEY}")
  string(REGEX REPLACE " "  "" PARDISO_LIC_KEY "${PARDISO_LIC_KEY}")
  set(PARDISO_LIC_KEY_CACHE "${PARDISO_LIC_KEY}" 
	  CACHE STRING "Pardiso license key (register at https://pardiso-project.org)" FORCE)

  string(LENGTH "${PARDISO_LIC_KEY}" KEY_LEN)

  if(NOT KEY_LEN EQUAL 56)
    message(WARNING "Pardiso license key '${PARDISO_LIC_KEY}' consists of ${KEY_LEN} characters, but the expected length is 56")
  endif()

  set(ENV{PARDISO_LIC_PATH} "${PARDISO_LIBRARY_DIR_CACHED}")
  message(STATUS "IMPORTANT: SET THE ENVIRONMENT VARIABLE PARDISO_LIC_PATH (probably \"${PARDISO_LIBRARY_DIR_CACHED}\")")
  message(STATUS "\tExample command for bash/zsh: 'export PARDISO_LIC_PATH=\"${PARDISO_LIBRARY_DIR_CACHED}\"'")
  file(WRITE ${PARDISO_LIBRARY_DIR_CACHED}/pardiso.lic "${PARDISO_LIC_KEY}\n")

elseif(PARDISO_LIBRARY_DIR)
  if((EXISTS "${PARDISO_LIBRARY_DIR}"))
    message(STATUS "pardiso directory '${PARDISO_LIBRARY_DIR}' used to find pardiso library")
    set(PARDISO_LIBRARY_DIR_CACHED "${PARDISO_LIBRARY_DIR}"
	  CACHE STRING "pardiso directory" FORCE)
  else()
      message(FATAL_ERROR "pardiso directory '${PARDISO_LIBRARY_DIR}' does not exist")
  endif()
else()
  message(STATUS "pardiso library is installed to a system directory")
endif()

find_package(Pardiso REQUIRED)
