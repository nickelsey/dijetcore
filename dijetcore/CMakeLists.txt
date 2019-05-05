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
add_subdirectory(proto)
add_subdirectory(worker)
add_subdirectory(test)
add_subdirectory(dijet_imbalance)
add_subdirectory(resources)
if (BUILD_DATA_QA_BINARIES)
  add_subdirectory(data)
endif(BUILD_DATA_QA_BINARIES)

## compile libraries

add_library(dijetcore ${DC_SRCS} ${DC_HDRS})
target_link_libraries(dijetcore PRIVATE ${DC_DEPENDENCY_LIBS})
if (TARGET dijetcore_proto)
  target_link_libraries(dijetcore PRIVATE dijetcore_proto)
endif(TARGET dijetcore_proto)
target_include_directories(dijetcore INTERFACE $<INSTALL_INTERFACE:include>)
target_compile_options(dijetcore INTERFACE "--std=c++17")
install(TARGETS dijetcore EXPORT DCTargets DESTINATION lib)

### declare libraries
set(DC_LIBS dijetcore)
if (TARGET dijetcore_proto)
  set(DC_LIBS ${DC_LIBS} dijetcore_proto)
endif(TARGET dijetcore_proto)

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
    target_link_libraries(${bin_name} ${DC_LIBS}
                          ${DC_DEPENDENCY_LIBS})
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${bin_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${bin_name} COMMAND $<TARGET_FILE:${bin_name}>)
    install(TARGETS ${bin_name} DESTINATION bin/${bin_path})
  endforeach()
endif()

## test binaries
if (BUILD_TEST)
  foreach(test_src ${DC_TEST_SRCS})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} "${test_src}")
    add_dependencies(${test_name} ${DC_LIBS})
    target_link_libraries(${test_name} ${DC_LIBS}
                    ${DC_DEPENDENCY_LIBS} gtest_main)
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${test_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    install(TARGETS ${test_name} DESTINATION test)
  endforeach()

  ## build gtest main
  add_executable(dijetcore_gtest_main ${DC_TEST_MAIN} ${DC_TEST_SRCS})
  add_dependencies(dijetcore_gtest_main ${DC_LIBS})
  target_link_libraries(dijetcore_gtest_main ${DC_LIBS} 
                        ${DC_DEPENDENCY_LIBS} gtest_main)
  #add_test(NAME dijetcore_gtest_main COMMAND $<TARGET_FILE:$dijetcore_gtest_main>)
  install(TARGETS dijetcore_gtest_main DESTINATION test)
endif()

## benchmark binaries
if (BUILD_TEST)
  foreach(bm_src ${DC_BENCH_SRCS})
    get_filename_component(bm_name ${bm_src} NAME_WE)
    add_executable(${bm_name} "${bm_src}")
    add_dependencies(${bm_name} ${DC_LIBS})
    target_link_libraries(${bm_name} ${DC_LIBS}
                  ${DC_DEPENDENCY_LIBS})
    if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.0)
      target_compile_features(${bm_name} PRIVATE cxx_range_for)
    endif()
    add_test(NAME ${bm_name} COMMAND $<TARGET_FILE:${bm_name}>)
    install(TARGETS ${bm_name} DESTINATION test)
  endforeach()
endif()

## export libraries to parent scope
set(DC_LIBS ${DC_LIBS} PARENT_SCOPE)
