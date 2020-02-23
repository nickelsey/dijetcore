########################################
## Set up ROOUNFOLD library environment
## ROOUNFOLD_PATH or ROOUNFOLDPATH have to be set as an environment variable
#  ROOUNFOLD_FOUND - System has TStarJetPico
#  ROOUNFOLD_INCLUDE_DIRS - The TStarJetPico include directories
#  ROOUNFOLD_LIBRARIES - The libraries of the TStarJetPico framework

# - Try to find TStarJetPico
# Once done this will define

find_path(ROOUNFOLD_INCLUDE_DIRS RooUnfold/RooUnfold.h
          HINTS $ENV{ROOUNFOLDPATH} $ENV{ROOUNFOLD_PATH}
                $ENV{ROOUNFOLDPATH}/include $ENV{ROOUNFOLD_PATH}/include)

find_library(ROOUNFOLD_LIBRARY
             NAMES RooUnfold
             HINTS $ENV{ROOUNFOLDPATH} $ENV{ROOUNFOLD_PATH}
             PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set STARPICO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(RooUnfold
           ROOUNFOLD_INCLUDE_DIRS ROOUNFOLD_LIBRARY)

mark_as_advanced(ROOUNFOLD_LIBRARY ROOUNFOLD_INCLUDE_DIR)

set(ROOUNFOLD_LIBRARIES ${ROOUNFOLD_LIBRARY})
set(ROOUNFOLD_INCLUDE_DIRS ${ROOUNFOLD_INCLUDE_DIRS})

## Add TStarJetPico include path to cmake
include_directories(${ROOUNFOLD_INCLUDE_DIRS})
