## first run over all proto files
file(GLOB dijetcore_PROTOBUF_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.proto")
if(dijetcore_PROTOBUF_FILES)
protobuf_generate_cpp(DC_PROTO_SRCS DC_PROTO_HDRS ${dijetcore_PROTOBUF_FILES})
endif(dijetcore_PROTOBUF_FILES)

add_library(dijetcore_proto STATIC ${DC_PROTO_HDRS} ${DC_PROTO_SRCS})

## get source files for libdijetcore
file(GLOB tmp *.cc)
set(DC_SRCS ${DC_SRCS} ${tmp})

## get header files for libdijetcore
file(GLOB tmp *.hh)
set(DC_HDRS ${DC_HDRS} ${tmp})

## get test source files
file(GLOB tmp *_test.cc)
set(DC_TEST_SRCS ${DC_TEST_SRCS} ${tmp})

## and exclude the test files from libdijetcore
if (DC_TEST_SRCS)
list(REMOVE_ITEM DC_SRCS ${DC_TEST_SRCS})
endif(DC_TEST_SRCS)

## get test benchmark files
file(GLOB tmp *_benchmark.cc)
set(DC_BENCH_SRCS ${DC_BENCH_SRCS} ${tmp})

## and exclude the benchmark files from libdijetcore
if (DC_BENCH_SRCS)
list(REMOVE_ITEM DC_SRCS ${DC_BENCH_SRCS})
endif(DC_BENCH_SRCS)

## install headers
install(FILES ${DC_PROTO_HDRS} DESTINATION include/dijetcore/lib)
install(FILES ${DC_HDRS} DESTINATION include/dijetcore/lib)

## add subdirectories
## none

## export file lists to parent scope
set(DC_SRCS ${DC_SRCS} PARENT_SCOPE)
set(DC_HDRS ${DC_HDRS} PARENT_SCOPE)
set(DC_TEST_SRCS ${DC_TEST_SRCS} PARENT_SCOPE)
set(DC_BENCH_SRCS ${DC_BENCH_SRCS} PARENT_SCOPE)
