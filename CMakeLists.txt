cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# set project name and version (MAJOR.minor.patch)
project(hmm_decoder VERSION 0.2.0 LANGUAGES CXX C)

# WINDOWS WARNING
if (WIN32)
message(WARNING "Build/Use of this library on Windows is untested and
  will likely not work properly. If this is the case,
  please file an issue and/or make a pull request @${PROJECT_URL} .")
endif()

# Specify C++ standard
# n.b. this has to be before adding targets
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
message(STATUS "CMAKE_CXX_STANDARD = ${CMAKE_CXX_STANDARD}")

# Add project cmake modules to module path
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
# provides dump_cmake_variables() fn
include(util)

# CMake-provided
# provides: CHECK_INCLUDE_FILE_CXX(<include> <variable> [<flags>])
include(CheckIncludeFileCXX)

# provides: CHECK_LIBRARY_EXISTS(LIBRARY FUNCTION LOCATION VARIABLE)
include(CheckLibraryExists)

# Let user know where things getting installed since they did not specify
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
  message(STATUS "*** CMAKE_INSTALL_PREFIX was initalised
    to the default value of ${CMAKE_INSTALL_PREFIX} ...")
endif()

# for RTXI modules, install "plugin" libraries to /usr/local/lib/rtxi
# (i.e., not sticking to GNU conventions exactly ...)

# # Adhere to GNU filesystem layout conventions
# # (see install config below)
# include(GNUInstallDirs)

################################################################################
############ OPTIONS
################################################################################
# Hard-coded for usage later with clangd:
set(CMAKE_EXPORT_COMPILE_COMMANDS true)

message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "")
message(STATUS "*** Looking for external libraries")

################################################################################
############ Find Libraries
################################################################################
# For each dependecy, check whether exists and get
# include directories (list PROJECT_REQUIRED_INCLUDE_DIRS),
# libraries (list PROJECT_REQUIRED_LIBRARIES), etc.
set(PROJECT_REQUIRED_CXX_FLAGS) #CXXFLAGS required by dependencies
set(PROJECT_REQUIRED_INCLUDE_DIRS ";") #include dirs ""
set(PROJECT_REQUIRED_LIBRARY_DIRS ";") #library dirs ""
set(PROJECT_REQUIRED_LIBRARIES_ABSOLUTE_NAME ";")#absolute path to libraries
set(PROJECT_REQUIRED_LIBRARIES_SHORT_NAME ";")#short name of library
# (i.e. short_name in `/path/to/lib*short_name*.so`)

# There may be differences in how you handle this as a function of OS
if(APPLE)
  message(STATUS "Detected macOS")
  # set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
elseif(UNIX)
  message(STATUS "Detected Unix")
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
  # # from armadillo re: linux (ommitting as not explicitly using mkl here):
  # # "necessary when linking with Intel MKL on Linux systems"
  # set(CMAKE_SHARED_LINKER_FLAGS
  # "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")
endif()

# pkg-config used for libraries that were not installed with CMake config files.
# provides `pkg_check_modules`
# find_package(PkgConfig REQUIRED)
find_package(Doxygen COMPONENTS dot
             OPTIONAL_COMPONENTS mscgen dia)

# RTXI
include(RTXI)

# HMM
include(hmm)

# save the CXX flags configured for later use by dependency.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROJECT_REQUIRED_CXX_FLAGS}")

# after getting these lists, take out duplicates
set(PROJECT_REQUIRED_LIBRARIES ${PROJECT_REQUIRED_LIBRARIES_ABSOLUTE_NAME})
list(APPEND PROJECT_REQUIRED_LIBRARIES ${PROJECT_REQUIRED_LIBRARIES_SHORT_NAME})
list(REMOVE_DUPLICATES PROJECT_REQUIRED_LIBRARIES)
list(REMOVE_DUPLICATES PROJECT_REQUIRED_INCLUDE_DIRS)
list(REMOVE_DUPLICATES PROJECT_REQUIRED_LIBRARY_DIRS)

message(STATUS "")
message(STATUS "*** Project will use the following libraries:")
message(STATUS "*** PROJECT_REQUIRED_LIBRARIES = ${PROJECT_REQUIRED_LIBRARIES}")
message(STATUS "")

# cf target_include_directories()
include_directories(${PROJECT_SOURCE_DIR}/include/
  ${PROJECT_REQUIRED_INCLUDE_DIRS})

# Set optimization level
# (have to check if a compiler that supports)
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
else()
  if(DEFINED CMAKE_CXX_COMPILER_ID)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    endif()
  endif()
endif()

# Finished configuring: print status
message(STATUS "CMAKE_CXX_FLAGS               = ${CMAKE_CXX_FLAGS}")
message(STATUS "PROJECT_REQUIRED_INCLUDE_DIRS = \
  ${PROJECT_REQUIRED_INCLUDE_DIRS}")

################################################################################
############ Add Targets
################################################################################
# library itself
add_library(${CMAKE_PROJECT_NAME} SHARED)
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC hmm ${PROJECT_REQUIRED_LIBRARIES})
target_link_directories(${CMAKE_PROJECT_NAME}
  PUBLIC ${PROJECT_REQUIRED_LIBRARY_DIRS})
target_include_directories(${CMAKE_PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES PREFIX "" VERSION
  ${CMAKE_PROJECT_VERSION_MAJOR}.${CMAKE_PROJECT_VERSION_MINOR}.${CMAKE_PROJECT_VERSION_PATCH}
  SOVERSION ${CMAKE_PROJECT_VERSION_MAJOR})

add_subdirectory(src) #add source files...

get_target_property(LIB_INCLUDES ${CMAKE_PROJECT_NAME} INCLUDE_DIRECTORIES)
message("LIB_INCLUDES: ${LIB_INCLUDES}")

# build doxygen documention
# build documentation
if(DOXYGEN_FOUND)
  set(DOXYGEN_FULL_PATH_NAMES "NO")
  set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/doxygen")
  set(DOXYGEN_HTML_OUTPUT "${PROJECT_BINARY_DIR}/doxygen/html")
  set(DOXYGEN_GENERATE_XML "YES")
  set(DOXYGEN_XML_OUTPUT "${PROJECT_BINARY_DIR}/doxygen/xml")
  # set(DOXYGEN_EXAMPLE_PATH "examples")
  set(DOXYGEN_EXAMPLE_PATTERNS "*.cpp")
  set(DOXYGEN_EXCLUDE_PATTERNS "*.m")

  doxygen_add_docs(doxygen
      ${PROJECT_SOURCE_DIR}/include/ ${PROJECT_SOURCE_DIR}/src/ ALL)
      # ${PROJECT_SOURCE_DIR}/examples/
elseif()
  message("Doxygen not found. Will not generate docs ... ")
endif()

################################################################################
############ Install Configuration
################################################################################
set(CMAKE_INSTALL_LIBDIR "lib/rtxi")
set(CMAKE_INSTALL_INCLUDEDIR "include/rtxi/plugins")

message(STATUS ""                                                      )
message(STATUS "CMAKE_INSTALL_PREFIX     = ${CMAKE_INSTALL_PREFIX}"    )
message(STATUS "CMAKE_INSTALL_LIBDIR     = ${CMAKE_INSTALL_LIBDIR}"    )
message(STATUS "CMAKE_INSTALL_INCLUDEDIR = ${CMAKE_INSTALL_INCLUDEDIR}")

# # add installation includedir
# list(APPEND PROJECT_REQUIRED_INCLUDE_DIRS
#   ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
# list(REMOVE_DUPLICATES PROJECT_REQUIRED_INCLUDE_DIRS)

# # Note that the trailing / character in "include/" is critical
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#   PATTERN ".svn" EXCLUDE
#   PATTERN "*.cmake" EXCLUDE
#   PATTERN "*~" EXCLUDE
#   PATTERN "*orig" EXCLUDE
#   )

install(
  TARGETS ${CMAKE_PROJECT_NAME}
  EXPORT ${CMAKE_PROJECT_NAME}#-> file <export-name>.cmake by default
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
