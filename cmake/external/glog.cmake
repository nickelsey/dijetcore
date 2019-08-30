if (NOT _GLOG_INCLUDE)
    set(_GLOG_INCLUDE TRUE)

    # if glog is installed on the system, we will use that
    find_package(glog)
    if (GLOG_FOUND)
        set(GLOG_EXTERNAL FALSE)
    else(GLOG_FOUND)
        # otherwise we build from github

        set(glog_PREFIX ${CMAKE_BINARY_DIR}/external/glog-build)
        #set(glog_INSTALL ${CMAKE_BINARY_DIR}/external/glog-install)
        set(glog_INSTALL ${CMAKE_INSTALL_PREFIX})

        set(glog_SOURCE_DIR ${CMAKE_BINARY_DIR}/external/glog-build/GLog)

        # we build statically and link into a shared object - requires position independent code
        if (UNIX)
            set(GLOG_EXTRA_COMPILER_FLAGS "-fPIC")
        endif()
        set(GLOG_EXTRA_COMPILER_FLAGS ${GLOG_EXTRA_COMPILER_FLAGS} "-m32")

        set(GLOG_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${GLOG_EXTRA_COMPILER_FLAGS})
        set(GLOG_C_FLAGS ${CMAKE_C_FLAGS} ${GLOG_EXTRA_COMPILER_FLAGS})

        # if gflags is being built internally, it must be built before glog
        if (GFLAGS_EXTERNAL)
            set(GLOG_DEPENDS GFlags)
        endif()

        # build the configure and autogen commands
        set(GLOG_AUTOGEN
            ./autogen.sh
            )
        set(GLOG_CONFIGURE 
            ./configure
            --prefix=${glog_INSTALL}
            --enable-shared=yes 
            --enable-static=no 
            --with-gflags=${GFLAGS_LIBRARY_DIRS}/..
            )
        
        ExternalProject_Add(GLog
            DEPENDS ${GLOG_DEPENDS}
            PREFIX ${glog_PREFIX}
            GIT_REPOSITORY "https://github.com/google/glog"
            GIT_TAG "v0.4.0"
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ${GLOG_AUTOGEN}
            COMMAND env "${GLOG_CONFIGURE}" "CFLAGS=\"${GLOG_C_FLAGS}\"" "CXXFLAGS=\"${GLOG_CXX_FLAGS}\""
            BUILD_COMMAND make
            INSTALL_COMMAND make install
            BUILD_IN_SOURCE 1
            LOG_DOWNLOAD 1
            LOG_CONFIGURE 1
            LOG_BUILD 1
            LOG_INSTALL 1
            )
            #INSTALL_DIR ${glog_INSTALL}

        #ExternalProject_Add(GLog
        #    DEPENDS ${GLOG_DEPENDS}
        #    PREFIX ${glog_PREFIX}
        #    GIT_REPOSITORY "https://github.com/google/glog"
        #    GIT_TAG "v0.4.0"
        #    UPDATE_COMMAND ""
        #    INSTALL_DIR ${glog_INSTALL}
        #    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        #    -DCMAKE_INSTALL_PREFIX=${glog_INSTALL}
        #    -DBUILD_SHARED_LIBS=OFF
        #    -DCMAKE_C_FLAGS=${GLOG_C_FLAGS}
        #    -DCMAKE_CXX_FLAGS=${GLOG_CXX_FLAGS}
        #    LOG_DOWNLOAD 1
        #    LOG_INSTALL 1
        #    )

        set(GLOG_FOUND TRUE)
        set(GLOG_INCLUDE_DIRS ${glog_INSTALL}/include)
        set(GLOG_LIBRARY_DIRS ${glog_INSTALL}/lib)
        set(GLOG_EXTERNAL TRUE)

        # identify libraries
        set(GLOG_LIBRARIES ${GFLAGS_LIBRARIES})
        find_library(lib glog PATHS ${GLOG_LIBRARY_DIRS} NO_DEFAULT_PATH)
        if(NOT lib)
            message(FATAL_ERROR "could not find glog library")
        else(NOT lib)
            set(GLOG_LIBRARIES ${GLOG_LIBRARIES} ${lib})
        endif(NOT lib)

    endif(GLOG_FOUND)

endif(NOT _GLOG_INCLUDE)
