## anaconda has quite a few packages that can be
## in conflict with the OS - this script attempts
## to avoid the problem, we always put anaconda
## include paths after system include paths

function(dc_include_directories)
  foreach(dir IN LISTS ARGN)
    if (${dir} MATCHES "/anaconda")
      include_directories(AFTER SYSTEM ${dir})
    else()
      include_directories(BEFORE SYSTEM ${dir})
    endif()
  endforeach()
endfunction()

function(build_dijetcore_binary bin_name)
  if (${ARGN})
    set(__target ${bin_name})
    prepend(__srcs "${CMAKE_CURRENT_SOURCE_DIR}/" "${ARGN}")
  else()
    get_filename_component(__target ${bin_name} NAME_WE)
    prepend(__srcs "${CMAKE_CURRENT_SOURCE_DIR}/" "${bin_name}")
  endif()
  add_executable(__target __srcs)
  target_link_libraries(__target ${SFD_LIBS})
  install(TARGETS ${__target} DESTINATION bin)
endfunction()
