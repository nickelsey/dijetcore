########################################
## Set up TStarJetPico library environment
## STARPICODIR or STARPICOPATH have to be set as an environment variable
#  STARPICO_FOUND - System has TStarJetPico
#  STARPICO_INCLUDE_DIRS - The TStarJetPico include directories
#  STARPICO_LIBRARIES - The libraries of the TStarJetPico framework

# - Try to find TStarJetPico
# Once done this will define

find_path(TSTARJETPICO_INCLUDE_DIR TStarJetPicoReader.h
          HINTS $ENV{STARPICOPATH} $ENV{STARPICODIR})

find_library(TSTARJETPICO_LIBRARY
             NAMES TStarJetPico
             HINTS $ENV{STARPICODIR} $ENV{STARPICOPATH})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set STARPICO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TStarJetPico
                                  TSTARJETPICO_INCLUDE_DIR TSTARJETPICO_LIBRARY)

mark_as_advanced(TSTARJETPICO_LIBRARY TSTARJETPICO_INCLUDE_DIR)

set(TSTARJETPICO_LIBRARIES ${TSTARJETPICO_LIBRARY})
set(TSTARJETPICO_INCLUDE_DIRS ${TSTARJETPICO_INCLUDE_DIR})

## Add TStarJetPico include path to cmake
include_directories(${TSTARJETPICO_INCLUDE_DIRS})
link_directories(${TSTARJETPICO_INCLUDE_DIRS})
