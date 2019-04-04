## get all files, assume that if they are in scripts
## directory, they should be installed into install_path/submit
file(GLOB SUBMITSRCS *)

## get test source files
file(GLOB CMAKEFILES *.txt)

## and exclude the CMakeLists file from install list
list(REMOVE_ITEM SUBMITSRCS ${CMAKEFILES})

install(PROGRAMS ${SUBMITSRCS} PERMISSIONS  DESTINATION submit)

