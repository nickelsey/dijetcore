#############################################
## tests, including simple tests to make sure
## everything has linked properly

## tests if ROOT include directory and libraries are set properly
set(ROOT_BASIC_TESTING_SRCS basic_root_test.cc)
add_executable(root_test ${ROOT_BASIC_TESTING_SRCS})
target_link_libraries(root_test ${ROOT_LIBRARIES})
set_target_properties(root_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

## tests if FastJet include directory and libraries are set properly
set(FASTJET_TESTING_SRCS basic_fastjet_test.cc)
add_executable(fastjet_test ${FASTJET_TESTING_SRCS})
target_link_libraries(fastjet_test ${FASTJET_LIBRARIES} ${ROOT_LIBRARIES})
set_target_properties(fastjet_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

## tests if TStarJetPico include directory and libraries are set properly
set(TS_TESTING_SRCS basic_ts_test.cc)
add_executable(ts_test ${TS_TESTING_SRCS})
target_link_libraries(ts_test ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES})
set_target_properties(ts_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(JET_DEF_TESTING_SRCS jet_def_test.cc)
add_executable(jet_def_test ${JET_DEF_TESTING_SRCS})
target_link_libraries(jet_def_test DijetWorker)
set_target_properties(jet_def_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(MATCH_DEF_TESTING_SRCS match_def_test.cc)
add_executable(match_def_test ${MATCH_DEF_TESTING_SRCS})
target_link_libraries(match_def_test DijetWorker)
set_target_properties(match_def_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(DIJET_DEFINITION_TESTING_SRCS dijet_definition_test.cc)
add_executable(dijet_def_test ${DIJET_DEFINITION_TESTING_SRCS})
target_link_libraries(dijet_def_test DijetWorker)
set_target_properties(dijet_def_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(DIJET_MATRIX_TESTING_SRCS dijet_matrix_test.cc)
add_executable(dijet_matrix_test ${DIJET_MATRIX_TESTING_SRCS})
target_link_libraries(dijet_matrix_test DijetWorker)
set_target_properties(dijet_matrix_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(DIJET_WORKER_TESTING_SRCS dijet_worker_test.cc)
add_executable(dijet_worker_test ${DIJET_WORKER_TESTING_SRCS})
target_link_libraries(dijet_worker_test DijetWorker)
set_target_properties(dijet_worker_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(ARG_HELPER_TESTING_SRCS arg_helper_test.cc)
add_executable(arg_helper_test ${ARG_HELPER_TESTING_SRCS})
set_target_properties(arg_helper_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(TRIGGER_LOOKUP_TESTING_SRCS trigger_lookup_test.cc)
add_executable(trig_lookup_test ${TRIGGER_LOOKUP_TESTING_SRCS})
set_target_properties(trig_lookup_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(SELECTOR_COMPARE_TESTING_SRCS selector_compare_test.cc)
add_executable(selector_compare_test ${SELECTOR_COMPARE_TESTING_SRCS})
target_link_libraries(selector_compare_test ${FASTJET_LIBRARIES})
set_target_properties(selector_compare_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(DIJET_KEY_TESTING_SRCS dijet_key_test.cc)
add_executable(dijet_key_test ${DIJET_KEY_TESTING_SRCS})
target_link_libraries(dijet_key_test DijetWorker)
set_target_properties(dijet_key_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

## reader util testing has a few text files it depends on,
## which are copied into the binary directory
configure_file(csv_test.csv ${CMAKE_BINARY_DIR}/jet_analysis/test/csv_test.csv)
configure_file(reader_util_test.in.txt ${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.txt)
configure_file(reader_util_test.in.txt ${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.txt)
configure_file(example_bad_run_list.txt ${CMAKE_BINARY_DIR}/jet_analysis/test/example_bad_run_list.txt)
configure_file(example_bad_tower_list.txt ${CMAKE_BINARY_DIR}/jet_analysis/test/example_bad_tower_list.txt)

## and also must be configured itself
configure_file(reader_util_test.in.cc ${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.cc)

set(READER_UTIL_TESTING_SRCS ${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.cc)
add_executable(reader_util_test ${READER_UTIL_TESTING_SRCS})
target_link_libraries(reader_util_test ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(reader_util_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(RUN4_EFF_TESTING_SRCS run4_eff_test.cc)
add_executable(run4_eff_test ${RUN4_EFF_TESTING_SRCS})
target_link_libraries(run4_eff_test ${ROOT_LIBRARIES} Efficiency)
set_target_properties(run4_eff_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)

set(PT_DIST_TESTING_SRCS pt_distribution_test.cc)
add_executable(pt_distribution_test ${PT_DIST_TESTING_SRCS})
target_link_libraries(pt_distribution_test ${ROOT_LIBRARIES})
set_target_properties(pt_distribution_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test)


#########################################
## Add tests
add_test(NAME star_basic COMMAND ts_test )
add_test(NAME fastjet_basic COMMAND fastjet_test)
add_test(NAME root_basic COMMAND root_test)
add_test(NAME arg_helper COMMAND arg_helper_test)
add_test(NAME trig_lookup COMMAND trig_lookup_test)
add_test(NAME jet_def COMMAND jet_def_test)
add_test(NAME match_def COMMAND match_def_test)
add_test(NAME dijet_definition COMMAND dijet_def_test)
add_test(NAME selector_compare COMMAND selector_compare_test)
add_test(NAME dijet_key COMMAND dijet_key_test)
add_test(NAME dijet_matrix COMMAND dijet_matrix_test)
add_test(NAME dijet_worker COMMAND dijet_worker_test)
add_test(NAME reader_util COMMAND reader_util_test)
add_test(NAME run4_eff COMMAND run4_eff_test)
add_test(NAME pt_dist_test COMMAND pt_distribution_test)

#########################################
## install targets
install(TARGETS ts_test fastjet_test root_test
         DESTINATION bin/test)
install(TARGETS jet_def_test arg_helper_test
         DESTINATION bin/test)
install(TARGETS match_def_test dijet_def_test
         DESTINATION bin/test)
install(TARGETS selector_compare_test dijet_key_test
         DESTINATION bin/test)
install(TARGETS dijet_matrix_test dijet_worker_test
         DESTINATION bin/test)
install(TARGETS reader_util_test run4_eff_test
         DESTINATION bin/test)
install(TARGETS pt_distribution_test
        DESTINATION bin/test)

## install the configured targets
install(FILES ${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.cc
        DESTINATION jet_analysis/test)
install(FILES ${CMAKE_BINARY_DIR}/jet_analysis/test/reader_util_test.txt
        DESTINATION jet_analysis/test)
