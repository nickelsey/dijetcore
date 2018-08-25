## install runtime resources into install directory

install(DIRECTORY bad_run_lists DESTINATION resources)
install(DIRECTORY bad_tower_lists DESTINATION resources)
install(DIRECTORY data_lists DESTINATION resources)
install(DIRECTORY efficiencies DESTINATION resources)
install(DIRECTORY runid_lists DESTINATION resources)

add_subdirectory(submit)
