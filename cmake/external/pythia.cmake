if (NOT _PYTHIA8_INCLUDE)
    set(_PYTHIA8_INCLUDE TRUE)

    # if pythia is installed on the system, we will use that
    find_package(Pythia8)

    IF(NOT PYTHIA8_FOUND)
        # otherwise we build from an external download

        # we exclude pythia8 targets from cmake clean
        set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM TRUE)
    
        # set the version of pythia we want
        set(PYTHIA8_VERSION 8303)
        
        set(PY8_FILENAME pythia${PYTHIA8_VERSION}.tgz)
        set(PY8_URL http://home.thep.lu.se/~torbjorn/pythia8/${PY8_FILENAME})
        set(PY8_DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/external/pythia8-source)
        set(PY8_SOURCE ${PY8_DOWNLOAD_DIR}/${PY8_FILENAME})
        set(PY8_SOURCE_DIR ${CMAKE_BINARY_DIR}/external/pythia${PYTHIA8_VERSION})
        set(PY8_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/pythia8)

        # check if cmake build type is debug
        if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
            set (configure_command ${configure_command} --enable-debug)
        endif ()

        if(UNIX)
            set(PY8_EXTRA_COMPILER_FLAGS "${PY8_EXTRA_COMPILER_FLAGS} -fPIC")
        endif(UNIX)

        if(BUILD_32BIT)
            set(PY8_EXTRA_COMPILER_FLAGS "${PY8_EXTRA_COMPILER_FLAGS} -m32")
        endif(BUILD_32BIT)

        set(PY8_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PY8_EXTRA_COMPILER_FLAGS}")
        set(PY8_C_FLAGS ${CMAKE_C_FLAGS} ${PY8_EXTRA_COMPILER_FLAGS})

        # build the configure command 
        set(PY8_CONFIGURE ./configure
            --prefix=${PY8_INSTALL_DIR}
            --cxx-common="${PY8_CXX_FLAGS} -O2 -std=c++11"
            )

        # define the external project
        ExternalProject_Add(Pythia8
            URL ${PY8_URL}
            DOWNLOAD_DIR ${PY8_DOWNLOAD_DIR}
            SOURCE_DIR ${PY8_SOURCE_DIR}
            CONFIGURE_COMMAND env "CFLAGS=${PY8_C_FLAGS}" ${PY8_CONFIGURE}
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make
            INSTALL_COMMAND make install
            )

        set(PYTHIA8_FOUND TRUE)
        set(PYTHIA8_INCLUDE_DIRS ${PY8_INSTALL_DIR}/include)
        set(PYTHIA8_LIBRARY_DIRS ${PY8_INSTALL_DIR}/lib)
        set(PYTHIA8_EXTERNAL TRUE)

        # now we have to identify the libraries
        set(PYTHIA8_LIBRARIES)
        set(REQ_LIB_NAMES pythia8)
        
        foreach(lib ${REQ_LIB_NAMES})
            list(APPEND PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY_DIRS}/lib${lib}.a)
        endforeach(lib ${REQ_LIB_NAMES})

        list(APPEND DC_EXTERNAL_DEPS Pythia8)

    endif(NOT PYTHIA8_FOUND)

endif(NOT _PYTHIA8_INCLUDE)
