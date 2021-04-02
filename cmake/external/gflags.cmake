if (NOT _GFLAGS_INCLUDE) 
  set(_GFLAGS_INCLUDE TRUE)

  # if GFlags is present on the system, use that
  find_package(GFlags)
  if (GFLAGS_FOUND)
    set(GFLAGS_EXTERNAL FALSE)
  # otherwise we will build it internally
  else(GFLAGS_FOUND)
    # gflags will use pthreads if the system supports it
    find_package(Threads)

    set(gflags_PREFIX ${CMAKE_BINARY_DIR}/external/gflags-build)
    set(gflags_INSTALL ${CMAKE_BINARY_DIR}/external/gflags-install)

    # if on a unix system, need to compile with -fPIC so that we can link into a shared library
    if(UNIX)
        set(GFLAGS_EXTRA_COMPILER_FLAGS "-fPIC")
    endif(UNIX)

    if(BUILD_32BIT)
        set(GFLAGS_EXTRA_COMPILER_FLAGS "${GFLAGS_EXTRA_COMPILER_FLAGS} -m32")
    endif(BUILD_32BIT)

    set(GFLAGS_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GFLAGS_EXTRA_COMPILER_FLAGS}")
    set(GFLAGS_C_FLAGS "${CMAKE_C_FLAGS} ${GFLAGS_EXTRA_COMPILER_FLAGS}")

    ExternalProject_Add(gflags
      PREFIX ${gflags_PREFIX}
      GIT_REPOSITORY "https://github.com/gflags/gflags.git"
      GIT_TAG "v2.2.2"
      UPDATE_COMMAND ""
      INSTALL_DIR ${gflags_INSTALL}
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_INSTALL_PREFIX=${gflags_INSTALL}
                 -DCMAKE_C_FLAGS=${GFLAGS_C_FLAGS}
                 -DCMAKE_CXX_FLAGS=${GFLAGS_CXX_FLAGS}
                 -DBUILD_SHARED_LIBS=OFF
                 -DBUILD_STATIC_LIBS=ON
                 -DBUILD_PACKAGING=OFF
                 -DBUILD_TESTING=OFF
                 -DBUILD_NC_TESTS=OFF
                 -BUILD_CONFIG_TESTS=OFF
                 -DINSTALL_HEADERS=ON
      LOG_DOWNLOAD 1
      LOG_INSTALL 1
      )

    set(GFLAGS_FOUND TRUE)
    set(GFLAGS_INCLUDE_DIRS ${gflags_INSTALL}/include)
    set(GFLAGS_LIBRARIES ${gflags_INSTALL}/lib/libgflags.a ${CMAKE_THREAD_LIBS_INIT})
    set(GFLAGS_LIBRARY_DIRS ${gflags_INSTALL}/lib)
    set(GFLAGS_EXTERNAL TRUE)

    list(APPEND DC_EXTERNAL_DEPS gflags)
  endif(GFLAGS_FOUND)

endif(NOT _GFLAGS_INCLUDE)
