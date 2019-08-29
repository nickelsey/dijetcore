#-------------------------------------------------------------------------------
# Some useful addons to ExternalProject module.
#
# Author(s):
#    Frank Tackmann
#
# Copyright:
#    Copyright (C) 2017 DESY
#
# Defines the following functions which define targets for the ExternalProject
# 'name' with the corresponding actions:
#
#  'ExternalProject_Add_Finish (name prefix_dir)'
#
#  Adds a step after installation to automatically rerun cmake with
#  name_ROOT set to prefix_dir
#
#  'ExternalProject_Add_Uninstall (name uninstall_command)'
#
#  Defines a target 'name-uninstall' which uses the given uninstall_command.
#  Note that the command should not depend on the SOURCE being present so it can
#  be called even after a call to 'name-distclean'.
#
#  'ExternalProject_Add_Clean (name clean_command)'
#
#  Defines a target 'name-clean' which uses the given clean_command with the
#  working directory set to the project's SOURCE_DIR.
#
#  'ExternalProject_Add_Distclean (name)'
#
#  Defines a target 'name-distclean' which removes the project's SOURCE_DIR but
#  keeps the downloaded tarball if it was stored somewhere else.
#
#  'ExternalProject_Add_Fullclean (name fullclean_command)'
#
#  Defines a target 'name-fullclean' which runs the given fullclean_command. The
#  target depends on the targets 'name-uninstall' and 'name-distclean', so these
#  will be run first.
#-------------------------------------------------------------------------------

function (ExternalProject_Add_Finish name prefix_dir)
   ExternalProject_Add_Step (${name} finish
      COMMAND ${CMAKE_COMMAND} -D${name}_ROOT=${prefix_dir} ${CMAKE_BINARY_DIR}
      COMMENT "Performing finish step for '${name}'"
      DEPENDEES install
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
   )
endfunction ()

function (ExternalProject_Add_Uninstall name)
   ExternalProject_Get_Property(${name} stamp_dir)
   _ep_get_configuration_subdir_suffix(cfgdir)

   add_custom_target (${name}-uninstall
      COMMAND ${ARGN}
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-install
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-custom-install
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-done
      COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${name}-complete
      COMMENT "Performing uninstall step for '${name}'"
      VERBATIM
   )
endfunction ()

function (ExternalProject_Add_Clean name)
   ExternalProject_Get_Property(${name} stamp_dir binary_dir)
   _ep_get_configuration_subdir_suffix(cfgdir)

   add_custom_target (${name}-clean
      COMMAND ${ARGN}
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-build
      WORKING_DIRECTORY ${binary_dir}
      COMMENT "Performing clean step for '${name}'"
      VERBATIM
   )
endfunction ()

function (ExternalProject_Add_Distclean name)
   ExternalProject_Get_Property(${name} stamp_dir source_dir binary_dir)
   _ep_get_configuration_subdir_suffix(cfgdir)

   add_custom_target (${name}-distclean
      COMMAND ${CMAKE_COMMAND} -E remove_directory ${source_dir}
      COMMAND ${CMAKE_COMMAND} -E remove_directory ${binary_dir}
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-mkdir
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-download
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-update
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-patch
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-configure
      COMMAND ${CMAKE_COMMAND} -E remove -f ${stamp_dir}${cfgdir}/${name}-build
      COMMENT "Performing distclean step for '${name}'"
      VERBATIM
   )
endfunction ()

function (ExternalProject_Add_Fullclean name)
   ExternalProject_Get_Property(${name} stamp_dir source_dir)
   _ep_get_configuration_subdir_suffix(cfgdir)

   add_custom_target (${name}-fullclean
      COMMAND ${ARGN}
      COMMENT "Performing fullclean step for '${name}'"
      VERBATIM
   )

   add_dependencies (${name}-fullclean ${name}-distclean ${name}-uninstall)
endfunction ()
