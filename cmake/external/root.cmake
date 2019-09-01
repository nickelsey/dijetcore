## ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT
    COMPONENTS MathCore
    RIO
    Hist
    Tree
    Net)

if(${ROOT_USE_FILE})
  include(${ROOT_USE_FILE})
endif(${ROOT_USE_FILE})

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
  MESSAGE(FATAL_ERROR "could not find root, and no external project is defined")
endif(NOT ROOT_FOUND)