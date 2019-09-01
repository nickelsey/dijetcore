## packages are either found or built, depending on if
## they are statically or dynamically linked

set(DC_DEPENDENCY_LIBS "")
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
find_package(Boost REQUIRED COMPONENTS filesystem)
dc_include_directories(${Boost_INCLUDE_DIRS})
list(APPEND DC_DEPENDENCY_LIBS ${Boost_LIBRARIES})

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