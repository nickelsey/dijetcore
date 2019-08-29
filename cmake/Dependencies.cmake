## packages are either found or built, depending on if
## they are statically or dynamically linked

set(DC_DEPENDENCY_LIBS "")
set(DC_EXTERNAL_DEPS "")
message(STATUS "fastjet")
## fastjet
include("cmake/external/fastjet.cmake")
dc_include_directories(${FASTJET_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${FASTJET_LIBRARIES})
list(APPEND DC_EXTERNAL_DEPS ${FASTJET_LIBRARIES})
message(STATUS "eventstructure")
## eventstructure
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/eventStructuredAu)
list(APPEND DC_DEPENDENCY_LIBS ${PICO_LIBS})
dc_include_directories(${PICO_INCLUDE_DIRS})
message(STATUS "boost")
## boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS filesystem)
dc_include_directories(${Boost_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${Boost_LIBRARIES})
message(STATUS "root")
## ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
list(APPEND DC_DEPENDENCY_LIBS ${ROOT_LIBRARIES})
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT")
message(STATUS "gflags")
## gflags
include("cmake/external/gflags.cmake")
dc_include_directories(${GFLAGS_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${GFLAGS_LIBRARIES})
list(APPEND DC_EXTERNAL_DEPS ${GFLAGS_LIBRARIES})
message(STATUS "glog")
## glog
include("cmake/external/glog.cmake")
dc_include_directories(${GLOG_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${GLOG_LIBRARIES})
list(APPEND DC_EXTERNAL_DEPS ${GLOG_LIBRARIES})

message(STATUS "testing")
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
