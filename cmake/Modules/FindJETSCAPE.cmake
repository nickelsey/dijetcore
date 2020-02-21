# Find the JetScape includes and library.
#
# This module defines
# JETSCAPE_INCLUDE_DIR   where to locate Pythia.h file
# JETSCAPE_LIBRARY       where to find the libJetScape library
# JETSCAPE_<lib>_LIBRARY Additional libraries
# JETSCAPE_LIBRARIES     (not cached) the libraries to link against to use JETSCAPE
# JETSCAPE_FOUND         if false, you cannot build anything that requires JETSCAPE
# JETSCAPE_VERSION       version of JETSCAPE if found

set(_jetscapedirs
    ${JETSCAPE}
    $ENV{JETSCAPE_DIR}
    $ENV{JETSCAPE}
    ${JETSCAPEDIR}
    $ENV{JETSCAPEDIR}
    /usr/local)

find_path(JETSCAPE_INCLUDE_DIR
          NAMES JETSCAPE/JetScape.h
          HINTS ${_jetscapedirs}
          PATH_SUFFIXES include
          DOC "Specify the directory containing Pythia.h.")

find_library(JETSCAPE_LIBRARY
             NAMES JetScape JETSCAPE jetscape
             HINTS ${_jetscapedirs}
             PATH_SUFFIXES lib
             DOC "Specify the JETSCAPE library here.")

find_library(JETSCAPE_reader_LIBRARY
             NAMES JetScapeReader jetscapereader
             HINTS ${_jetscapedirs}
             PATH_SUFFIXES lib)

find_library(JETSCAPE_third_LIBRARY
             NAMES JetScapeThird jetscapethird
             HINTS ${_jetscapedirs}
             PATH_SUFFIXES lib)

find_library(JETSCAPE_gtl_LIBRARY
             NAMES GTL gtl
             HINTS ${_jetscapedirs}
             PATH_SUFFIXES lib)

find_library(JETSCAPE_trento_LIBRARY
             NAMES trento Trento
             HINTS ${_jetscapedirs}
             PATH_SUFFIXES lib)

find_libraryh(JETSCAPE_hydrofromfile_LIBRARY
              NAMES hydroFromFile hydrofromfile
              HINTS ${_jetscapedirs}
              PATH_SUFFIXES lib)

foreach(_lib JETSCAPE_LIBRARY PYTHIA8_hepmcinterface_LIBRARY PYTHIA8_lhapdfdummy_LIBRARY)
  if(${_lib})
    set(JETSCAPE_LIBRARIES ${PYTHIA8_LIBRARIES} ${${_lib}})
    set(JETSCAPE_FOUND TRUE)
  endif()
endforeach()
set(JETSCAPE_INCLUDE_DIRS ${PYTHIA8_INCLUDE_DIR} ${PYTHIA8_INCLUDE_DIR}/JETSCAPE )

if(JETSCAPE_FOUND)
  message (STATUS "Found JETSCAPE")
else(JETSCAPE_FOUND)
  message(FATAL_ERROR "Didn't find JETSCAPE")
endif(JETSCAPE_FOUND)

# Add JetScape include path to cmake
include_directories(${JETSCAPE_INCLUDE_DIRS})

