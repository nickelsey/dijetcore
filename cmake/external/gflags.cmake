if (NOT _GFLAGS_INCLUDE) 
    set(_GFLAGS_INCLUDE TRUE)

    # look if gflags has been installed on the system
    find_package(GFlags)
    if (GFLAGS_FOUND)
        set(GFLAGS_EXTERNAL FALSE)
    else(GFLAGS_FOUND)
        # otherwise we build from scratch

        # gflags will use pthreads if its available
        find_package(Threads)

        set(gflags_PREFIX ${CMAKE_BINARY_DIR}/external/gflags-build)
        set(gflags_INSTALL ${CMAKE_BINARY_DIR}/external/gflags-install)

        # we build statically and link into a shared object - requires position independent code
        if (UNIX)
            set(GFLAGS_EXTRA_COMPILER_FLAGS "-fPIC")
        endif()

        set(GFLAGS_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${GFLAGS_EXTRA_COMPILER_FLAGS})
        set(GFLAGS_C_FLAGS ${CMAKE_C_FLAGS} ${GFLAGS_EXTRA_COMPILER_FLAGS})

        ExternalProject_Add(gflags
            PREFIX ${gflags_PREFIX}
            GIT_REPOSITORY "https://github.com/gflags/gflags.git"
            GIT_TAG "v2.2.2"
            UPDATE_COMMAND ""
            INSTALL_DIR ${gflags_INSTALL}
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${gflags_INSTALL}
            -DBUILD_SHARED_LIBS=OFF
            -DBUILD_STATIC_LIBS=ON
            -DBUILD_PACKAGING=OFF
            -DBUILD_TESTING=OFF
            -DBUILD_NC_TESTS=OFF
            -BUILD_CONFIG_TESTS=OFF
            -DINSTALL_HEADERS=ON
            -DCMAKE_C_FLAGS=${GFLAGS_C_FLAGS}
            -DCMAKE_CXX_FLAGS=${GFLAGS_CXX_FLAGS}
            LOG_DOWNLOAD 1
            LOG_INSTALL 1
            )

        set(GFLAGS_FOUND TRUE)
        set(GFLAGS_INCLUDE_DIRS ${gflags_INSTALL}/include)
        set(GFLAGS_LIBRARIES ${gflags_INSTALL}/lib/libgflags.a ${CMAKE_THREAD_LIBS_INIT})    
        set(GFLAGS_LIBRARY_DIRS ${gflags_INSTALL}/lib)
        set(GFLAGS_EXTERNAL TRUE)

    endif(GFLAGS_FOUND)

endif(NOT _GFLAGS_INCLUDE)
