if (NOT _FASTJET_INCLUDE)
    set(_FASTJET_INCLUDE TRUE)

    # if fastjet is installed on the system, we will use that
    find_package(fastjet)
    if (FASTJET_FOUND)
        set(FASTJET_EXTERNAL FALSE)
        find_program(FJ_CONFIG fastjet-config)
        
        execute_process(
            COMMAND ${FJ_CONFIG} --version 
            OUTPUT_VARIABLE FJ_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        
        string(SUBSTRING ${FJ_VERSION} 0 1 FJ_MAJOR_VERSION)
        string(SUBSTRING ${FJ_VERSION} 2 1 FJ_MINOR_VERSION)
        string(SUBSTRING ${FJ_VERSION} 4 1 FJ_PATCH_VERSION)

        if (FJ_MAJOR_VERSION STRLESS 2 OR FJ_MINOR_VERSION STRLESS 3)
            set(FASTJET_FOUND FALSE)
        endif(FJ_MAJOR_VERSION STRLESS 2 OR FJ_MINOR_VERSION STRLESS 3)
    endif(FASTJET_FOUND)

    IF(NOT FASTJET_FOUND)
        # otherwise we build from github

        # we exclude fastjet targets from cmake clean
        set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM TRUE)

        # set the version of fastjet we want
        set(FASTJET_VERSION 3.3.2)


        set(FJ_MD5 ca3708785c9194513717a54c1087bfb0)
        set(FJ_FILENAME fastjet-${FASTJET_VERSION}.tar.gz)
        set(FJ_URL http://fastjet.fr/repo/${FJ_FILENAME})
        set(FJ_DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/external/fastjet_source)
        set(FJ_SOURCE ${FJ_DOWNLOAD_DIR}/${FJ_FILENAME})
        set(FJ_SOURCE_DIR ${CMAKE_BINARY_DIR}/external/fastjet-${FASTJET_VERSION})
        set(FJ_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/fastjet)

        # build the configure command 
        set(FJ_CONFIGURE ./configure
            --prefix=${FJ_INSTALL_DIR}
            --with-pic
            --disable-auto-ptr
            )

        # check if cmake build type is debug
        if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
            set (configure_command ${configure_command} --enable-debug)
        endif ()

        # define the external project
        ExternalProject_Add(FastJet
            URL ${FJ_URL}
            URL_MD5 ${FJ_MD5}
            DOWNLOAD_DIR ${FJ_DOWNLOAD_DIR}
            SOURCE_DIR ${FJ_SOURCE_DIR}
            CONFIGURE_COMMAND ${FJ_CONFIGURE}
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make
            INSTALL_COMMAND make install
            )

        ExternalProject_Add_Finish(FastJet ${FJ_INSTALL_DIR})
        ExternalProject_Get_Property(FastJet stamp_dir)

        # custom uninstall script for fastjet
        file(WRITE ${stamp_dir}/FastJet-uninstall.cmake
            "
            set(bins ${FJ_INSTALL_DIR}/bin/fastjet-config)
            foreach(bin ${bins})
                message (STATUS \"Removing: \${bin}\")
                file (REMOVE ${bin})
            endforeach(bin)
            file(GLOB libs ${FJ_INSTALL_DIR}/lib/libfastjet* ${FJ_INSTALL_DIR}/lib/libsiscone*)
            foreach (lib ${libs})
                message (STATUS \"Removing: \${lib}\")
                file (REMOVE ${lib})
            endforeach (lib)
            message (STATUS \"Removing: \${FJ_INSTALL_DIR}/include/fastjet\")
            file (REMOVE_RECURSE ${FJ_INSTALL_DIR}/include/fastjet)
            "
            )
        # add default uninstall and clean targets
        ExternalProject_Add_Uninstall(FastJet
            ${CMAKE_COMMAND} -P ${stamp_dir}/FastJet-uninstall.cmake
            )
        ExternalProject_Add_Clean(FastJet
            make clean
            )
        ExternalProject_Add_Distclean(FastJet)
        ExternalProject_Add_Fullclean(FastJet
            ${CMAKE_COMMAND} -E remove -f ${source_file}
            )

        set(FASTJET_FOUND TRUE)
        set(FASTJET_INCLUDE_DIRS ${FJ_INSTALL_DIR}/include)
        set(FASTJET_LIBRARY_DIRS ${FJ_INSTALL_DIR}/lib)
        set(FASTJET_EXTERNAL TRUE)

        # now we have to identify the libraries
        set(FASTJET_LIBRARIES)
        set(REQ_LIB_NAMES fastjet fastjettools fastjetplugins)
        set(OPT_LIB_NAMES siscone siscone_spherical)

        foreach(lib_name ${REQ_LIB_NAMES})
            unset(lib CACHE)
            find_library(lib ${lib_name} PATHS ${FASTJET_LIBRARY_DIRS} NO_DEFAULT_PATH)
            find_library(lib ${lib_name})
            if(lib)
                set(FASTJET_LIBRARIES ${FASTJET_LIBRARIES} ${lib_name})
            ELSE(lib)
                message(FATAL_ERROR "Can't find FastJet library ${lib_name}")
            endif(lib)
        endforeach(lib_name)
        foreach(lib_name ${OPT_LIB_NAMES})
            unset(lib CACHE)
            find_library(lib ${lib_name} PATHS ${FASTJET_LIBRARY_DIRS} NO_DEFAULT_PATH)
            find_library(lib ${lib_name})
            if(lib)
                set(FASTJET_LIBRARIES ${FASTJET_LIBRARIES} ${lib_name})
            ELSE(lib)
                message(STATUS "Can't find optional FastJet library ${lib_name}")
            endif(lib)
        endforeach(lib_name)

    endif(NOT FASTJET_FOUND)

endif(NOT _FASTJET_INCLUDE)
