# make sure we include gflags, because glog depends on it
include("cmake/external/gflags.cmake")

if (NOT _GLOG_INCLUDE)
  set(_GLOG_INCLUDE TRUE)

  # if GLog is present on the system, use that
  find_package(Glog)
  if (GLOG_FOUND)
      set(GLOG_EXTERNAL FALSE)
  # otherwise, we will build from source
  else(GLOG_FOUND)
    
    set(glog_PREFIX ${CMAKE_BINARY_DIR}/external/glog-build)
    set(glog_INSTALL ${CMAKE_BINARY_DIR}/external/glog-install)

    # if on a unix system, need to compile with -fPIC so that we can link into a shared library
    if (UNIX)
      set(GLOG_EXTRA_COMPILER_FLAGS "-fPIC")
    endif()

    if(BUILD_32BIT)
        set(GLOG_EXTRA_COMPILER_FLAGS "${GLOG_EXTRA_COMPILER_FLAGS} -m32")
    endif(BUILD_32BIT)

    set(GLOG_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLOG_EXTRA_COMPILER_FLAGS}")
    set(GLOG_C_FLAGS "${CMAKE_C_FLAGS} ${GLOG_EXTRA_COMPILER_FLAGS}")

    # depend on gflags if we're also building it
    if (GFLAGS_EXTERNAL)
      set(GLOG_DEPENDS gflags)
    endif()

    ExternalProject_Add(glog
      DEPENDS ${GLOG_DEPENDS}
      PREFIX ${glog_PREFIX}
      GIT_REPOSITORY "https://github.com/google/glog"
      GIT_TAG "v0.4.0"
      UPDATE_COMMAND ""
      INSTALL_DIR ${gflags_INSTALL}
      PATCH_COMMAND autoreconf -i ${glog_PREFIX}/src/glog
      CONFIGURE_COMMAND env "CFLAGS=${GLOG_C_FLAGS}" "CXXFLAGS=${GLOG_CXX_FLAGS}" ${glog_PREFIX}/src/glog/configure --prefix=${glog_INSTALL} --enable-shared=no --enable-static=yes --with-gflags=${GFLAGS_LIBRARY_DIRS}/..
      LOG_DOWNLOAD 1
      LOG_CONFIGURE 1
      LOG_INSTALL 1
      )

    set(GLOG_FOUND TRUE)
    set(GLOG_INCLUDE_DIRS ${glog_INSTALL}/include)
    set(GLOG_LIBRARIES ${GFLAGS_LIBRARIES} ${glog_INSTALL}/lib/libglog.a)
    set(GLOG_LIBRARY_DIRS ${glog_INSTALL}/lib)
    set(GLOG_EXTERNAL TRUE)

    list(APPEND DC_EXTERNAL_DEPS glog)
  endif(GLOG_FOUND)

endif(NOT _GLOG_INCLUDE)

