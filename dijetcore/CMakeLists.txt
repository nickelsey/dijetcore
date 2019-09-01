## declare our source files
## DC_SRCS is used to build libdijetcore
## DC_TEST_SRCS contains all test sources
## DC_TEST_MAIN contains gtest main()
## DC_BENCH_SRCS contains benchmark routine sources
## DC_BIN_SRCS contains binary source files

set(DC_SRCS)
set(DC_HDRS)
set(DC_TEST_SRCS)
set(DC_TEST_MAIN)
set(DC_BENCH_SRCS)
set(DC_BINARY_SRCS)

## add all subdirectories
add_subdirectory(lib)
add_subdirectory(util)
add_subdirectory(worker)
add_subdirectory(dijet_imbalance)
add_subdirectory(resources)

## compile libraries

add_library(dijetcore ${DC_SRCS})
target_link_libraries(dijetcore PRIVATE ${DC_DEPENDENCY_LIBS} ${DC_LINKER_LIBS})
target_include_directories(dijetcore INTERFACE $<INSTALL_INTERFACE:include>)
target_compile_options(dijetcore INTERFACE "--std=c++11")
install(TARGETS dijetcore EXPORT DCTargets DESTINATION lib)

## dijetcore depends on our external projects
if(DC_EXTERNAL_DEPS)
  add_dependencies(dijetcore ${DC_EXTERNAL_DEPS})
endif(DC_EXTERNAL_DEPS)

### declare libraries
set(DC_LIBS dijetcore)

## build binaries
if (BUILD_BINARIES)
  foreach(bin_src ${DC_BINARY_SRCS})
    get_filename_component(bin_name ${bin_src} NAME_WE)

    ### extract the directory
    get_filename_component(dir_name ${bin_src} DIRECTORY)
    ## extract directory relative to dijetcore source structure
    string(REGEX MATCH "${CMAKE_SOURCE_DIR}/dijetcore/([^ ]*)" TMP_PATH ${dir_name})
    set(bin_path ${CMAKE_MATCH_1})

    add_executable(${bin_name} "${bin_src}")
    add_dependencies(${bin_name} ${DC_LIBS})
    target_link_libraries(${bin_name}
                          ${DC_LIBS}
                          ${DC_LINKER_LIBS})
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${bin_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${bin_name} COMMAND $<TARGET_FILE:${bin_name}>)
    install(TARGETS ${bin_name} DESTINATION bin/${bin_path})
  endforeach()
endif()

## export libraries to parent scope
set(DC_LIBS ${DC_LIBS} PARENT_SCOPE)
