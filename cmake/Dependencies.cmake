## packages are either found or built, depending on if
## they are statically or dynamically linked

set(DC_DEPENDENCY_LIBS "")

## fastjet
find_package(FastJet REQUIRED)
list(APPEND DC_DEPENDENCY_LIBS ${FASTJET_LIBRARIES})

## pythia
find_package(Pythia8 REQUIRED)
list(APPEND DC_DEPENDENCY_LIBS ${PYTHIA8_LIBRARIES})

## eventstructure
find_package(TStarJetPico REQUIRED)
list(APPEND DC_DEPENDENCY_LIBS ${TSTARJETPICO_LIBRARIES})

## boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS filesystem)
dc_include_directories(${Boost_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${Boost_LIBRARIES})

## ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
list(APPEND DC_DEPENDENCY_LIBS ${ROOT_LIBRARIES})
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT")

## gflags
find_package(gflags)
if(GFLAGS_FOUND)
  dc_include_directories(${GFLAGS_INCLUDE_DIRS})
  list(APPEND DC_DEPENDENCY_LIBS ${GFLAGS_LIBRARIES})
else(GFLAGS_FOUND)
  message(FATAL_ERROR "gflags library not found")
endif(GFLAGS_FOUND)

## glog
find_package(glog)
if(GLOG_FOUND)
  set(DC_USE_GLOG 1)
  dc_include_directories(${GLOG_INCLUDE_DIRS})
  list(APPEND DC_DEPENDENCY_LIBS ${GLOG_LIBRARIES})
else(GLOG_FOUND)
  message(FATAL_ERROR "glog library not found")
endif(GLOG_FOUND)

## protobuf
find_package(Protobuf)
if (Protobuf_FOUND)
  dc_include_directories(${Protobuf_INCLUDE_DIRS})
  list(APPEND SFD_DEPENDENCY_LIBS ${Protobuf_LIBRARIES})
else(Protobuf_FOUND)
  message(FATAL_ERROR "protobuf library not found")
endif(Protobuf_FOUND)

## testing is done via gtest, gmock (currently not used)
## and google benchmark. They are compiled as static libraries
## and embedded in the test binaries
if(BUILD_TEST)
  set(TEMP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(BUILD_SHARED_LIBS OFF)
  set(BUILD_GTEST ON)
  set(INSTALL_GTEST OFF)
  ## gmock currently not used
  set(BUILD_GMOCK OFF)
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
