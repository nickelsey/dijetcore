## packages are either found or built, depending on if
## they are statically or dynamically linked

set(DC_DEPENDENCY_LIBS "")
set(DC_TEST_LIBS "")
set(DC_LINKER_LIBS "")
set(DC_EXTERNAL_DEPS "")

## ROOT
include("cmake/external/root.cmake")
link_directories(${ROOT_LIBRARY_DIRS})
dc_include_directories(${ROOT_INCLUDE_DIRS})
list(APPEND DC_LINKER_LIBS ${ROOT_LIBRARIES})

## eventstructure
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/eventStructuredAu)
list(APPEND DC_LINKER_LIBS ${PICO_LIBS})
dc_include_directories(${PICO_INCLUDE_DIRS})

## boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem)
dc_include_directories(${Boost_INCLUDE_DIRS})
list(APPEND DC_LINKER_LIBS ${Boost_LIBRARIES})

## glog
include("cmake/external/glog.cmake")
dc_include_directories(${GLOG_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS PUBLIC ${GLOG_LIBRARIES})

## gflags
include("cmake/external/gflags.cmake")
dc_include_directories(${GFLAGS_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS PUBLIC ${GFLAGS_LIBRARIES})

## fastjet
include("cmake/external/fastjet.cmake")
dc_include_directories(${FASTJET_INCLUDE_DIRS})
link_directories(${FASTJET_LIBRARY_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${FASTJET_LIBRARIES})

if(BUILD_TEST)
  set(TEMP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(TEMP_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  set(BUILD_SHARED_LIBS OFF)
  set(BUILD_GTEST ON CACHE BOOL "build core gtest")
  set(INSTALL_GTEST OFF CACHE BOOL "do not install gtest to install directory")
  if(BUILD_32BIT)
    add_compile_options(-m32)
    set(CMAKE_CXX_FLAGS "-m32" CACHE STRING "set 32 bit for googletest and benchmark")
  endif(BUILD_32BIT)
  ## gmock currently not used
  set(BUILD_GMOCK OFF CACHE BOOL "do not build gmock")
  add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/googletest)
  dc_include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/googletest/include)
  list(APPEND DC_TEST_LIBS gtest)

  # We will not need to test benchmark lib itself.
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing.")
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Disable benchmark install to avoid overwriting.")
  add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/benchmark)
  dc_include_directories(${PROJECT_SOURCE_DIR}/third_party/benchmark/include)
  list(APPEND DC_TEST_LIBS benchmark)

  # restore the build shared libs option.
  set(BUILD_SHARED_LIBS ${TEMP_BUILD_SHARED_LIBS})
  set(CMAKE_CXX_FLAGS ${TEMP_CXX_FLAGS})
endif(BUILD_TEST)