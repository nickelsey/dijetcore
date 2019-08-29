########################################
# Set up Fastjet environment
# requires fastjet-config to be in $PATH
# this find script is a modified version of the SpartyJet fastjet cmake package

EXECUTE_PROCESS ( COMMAND which fastjet-config OUTPUT_VARIABLE FASTJET_CONFIG )

if ( FASTJET_CONFIG )
    EXECUTE_PROCESS ( COMMAND fastjet-config --prefix OUTPUT_VARIABLE FASTJET_PREFIX )
    if ( FASTJET_PREFIX )
        MESSAGE ( STATUS "Found Fastjet" )
        STRING ( STRIP ${FASTJET_PREFIX} FASTJET_PREFIX )
        SET ( FASTJET_LIB_PATH ${fastjet_LIB_PATH} ${FASTJET_PREFIX}/lib )
        SET ( FASTJET_INCLUDE_DIRS ${FASTJET_INCLUDE_DIRS} ${FASTJET_PREFIX}/include )
        SET ( FASTJET_FOUND ON )

        # We'll set up all the libraries, including optional libraries
        SET ( REQ_LIB_NAMES fastjet fastjettools fastjetplugins )
        SET ( OPT_LIB_NAMES siscone siscone_spherical )

        FOREACH ( lib_name ${REQ_LIB_NAMES} )
            UNSET (lib CACHE)
            FIND_LIBRARY ( lib ${lib_name} PATHS ${FASTJET_LIB_PATH} NO_DEFAULT_PATH)
            FIND_LIBRARY ( lib ${lib_name} )
            IF ( lib )
                SET ( FASTJET_LIBRARIES ${FASTJET_LIBRARIES} ${lib_name} )
            ELSE ( lib )
                MESSAGE ( FATAL_ERROR "Can't find FastJet library ${lib_name}" )
            ENDIF ( lib )
        ENDFOREACH ( lib_name )

        FOREACH ( lib_name ${OPT_LIB_NAMES} )
            UNSET (lib CACHE)
            FIND_LIBRARY ( lib ${lib_name} PATHS ${FASTJET_LIB_PATH} NO_DEFAULT_PATH)
            FIND_LIBRARY ( lib ${lib_name} )
            IF ( lib )
                SET ( FASTJET_LIBRARIES ${FASTJET_LIBRARIES} ${lib_name} )
            ELSE ( lib )
                MESSAGE ( FATAL_ERROR "Can't find FastJet library ${lib_name}" )
            ENDIF ( lib )
        ENDFOREACH ( lib_name )

        LINK_DIRECTORIES ( ${FASTJET_LIB_PATH} )

    ENDIF ( FASTJET_PREFIX )
ENDIF ( FASTJET_CONFIG )
