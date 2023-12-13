if(TARGET slot_map::slot_map)
    return()
endif()

message(STATUS "Third-party (external): creating target 'slot_map::slot_map'")

include(CPM)
CPMAddPackage(
  NAME slot_map
  GITHUB_REPOSITORY qnzhou/slot_map
  GIT_TAG 98717cc2b63c3866a00272894fcf4ea8ce1fe484
  DOWNLOAD_ONLY YES
)

add_subdirectory(${slot_map_SOURCE_DIR}/slot_map ${slot_map_BINARY_DIR})
set_target_properties(slot_map PROPERTIES SYSTEM ON)
set_target_properties(slot_map PROPERTIES FOLDER third_party)
add_library(slot_map::slot_map ALIAS slot_map)
