## packages are either found or built, depending on if
## they are statically or dynamically linked

set(DC_DEPENDENCY_LIBS "")
set(DC_EXTERNAL_DEPS "")

## ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT
    COMPONENTS MathCore
    RIO
    Hist
    Tree
    Net)

if(NOT ROOT_FOUND)
    # we will look using root-config
    find_program(ROOT_CONFIG root-config PATHS
        ${ROOTSYS}/bin
        )

    if(ROOT_CONFIG)
        set(ROOT_FOUND TRUE)

        execute_process(
            COMMAND ${ROOT_CONFIG} --prefix
            OUTPUT_VARIABLE ROOT_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )

        execute_process(
            COMMAND ${ROOT_CONFIG} --version 
            OUTPUT_VARIABLE ROOT_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        execute_process(
            COMMAND ${ROOT_CONFIG} --incdir
            OUTPUT_VARIABLE ROOT_INCLUDE_DIRS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )

        execute_process(
            COMMAND ${ROOT_CONFIG} --libs
            OUTPUT_VARIABLE ROOT_LIBRARIES
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )

        execute_process(
            COMMAND ${ROOT_CONFIG} --libdir
            OUTPUT_VARIABLE ROOT_LIBRARY_DIRS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )

        string(SUBSTRING ${ROOT_VERSION} 0 1 ROOT_MAJOR_VERSION)
        string(SUBSTRING ${ROOT_VERSION} 2 2 ROOT_MINOR_VERSION)
        string(SUBSTRING ${ROOT_VERSION} 5 2 ROOT_PATCH_VERSION)

        if(${ROOT_MAJOR_VERSION} STREQUAL 5)
          set(ROOT_V5 TRUE)
        else(${ROOT_MAJOR_VERSION} STREQUAL 5)
          set(ROOT_V5 FALSE)
        endif(${ROOT_MAJOR_VERSION} STREQUAL 5)

        if(${ROOT_MAJOR_VERSION} STREQUAL 6)
          set(ROOT_V6 TRUE)
        else(${ROOT_MAJOR_VERSION} STREQUAL 6)
          set(ROOT_V6 FALSE)
        endif(${ROOT_MAJOR_VERSION} STREQUAL 6)
    endif(ROOT_CONFIG)
endif(NOT ROOT_FOUND)

if(NOT ROOT_FOUND)
  MESSAGE(FATAL_ERROR "could not find root")
endif(NOT ROOT_FOUND)

link_directories(${ROOT_LIBRARY_DIRS})
dc_include_directories(${ROOT_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${ROOT_LIBRARIES})
if(${ROOT_USE_FILE})
  include(${ROOT_USE_FILE})
endif(${ROOT_USE_FILE})
message(STATUS "Found ROOT")

## boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS filesystem)
dc_include_directories(${Boost_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${Boost_LIBRARIES})

## fastjet
include("cmake/external/fastjet.cmake")
dc_include_directories(${FASTJET_INCLUDE_DIRS})
link_directories(${FASTJET_LIBRARY_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${FASTJET_LIBRARIES})
list(APPEND DC_EXTERNAL_DEPS ${FASTJET_LIBRARIES})
message(STATUS "fastjet libs: ")
message(STATUS ${DC_DEPENDENCY_LIBS})
message(STATUS ${DC_EXTERNAL_DEPS})

## eventstructure
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/eventStructuredAu)
list(APPEND DC_DEPENDENCY_LIBS ${PICO_LIBS})
dc_include_directories(${PICO_INCLUDE_DIRS})

## gflags
include("cmake/external/gflags.cmake")
dc_include_directories(${GFLAGS_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${GFLAGS_LIBRARIES})
list(APPEND DC_EXTERNAL_DEPS ${GFLAGS_LIBRARIES})

## glog
include("cmake/external/glog.cmake")
dc_include_directories(${GLOG_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${GLOG_LIBRARIES})
list(APPEND DC_EXTERNAL_DEPS ${GLOG_LIBRARIES})

## testing is done via gtest, gmock (currently not used)
## and google benchmark. They are compiled as static libraries
## and embedded in the test binaries
if(BUILD_TEST)
  set(TEMP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(BUILD_SHARED_LIBS OFF)
  set(BUILD_GTEST ON CACHE BOOL "build core gtest")
  set(INSTALL_GTEST OFF CACHE BOOL "do not install gtest to install directory")
  ## gmock currently not used
  set(BUILD_GMOCK OFF CACHE BOOL "do not build gmock")
  add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/googletest)
  dc_include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/googletest/include)

  # We will not need to test benchmark lib itself.
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing.")
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Disable benchmark install to avoid overwriting.")
  add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/benchmark)
  dc_include_directories(${PROJECT_SOURCE_DIR}/third_party/benchmark/include)
  list(APPEND DC_DEPENDENCY_LIBS benchmark)

  # restore the build shared libs option.
  set(BUILD_SHARED_LIBS ${TEMP_BUILD_SHARED_LIBS})
endif()

## set more verbose variable names
set(DC_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(DC_BUILD_TEST ${BUILD_TEST})
