if(TARGET slot_map::slot_map)
    return()
endif()

message(STATUS "Third-party (external): creating target 'slot_map::slot_map'")

include(CPM)
CPMAddPackage(
  NAME slot_map
  GITHUB_REPOSITORY SergeyMakeev/slot_map
  GIT_TAG 4d6c2ba92d55d4726154a290304d10aefd09a89b
  DOWNLOAD_ONLY YES
)

add_subdirectory(${slot_map_SOURCE_DIR}/slot_map)
set_target_properties(slot_map PROPERTIES SYSTEM ON)
set_target_properties(slot_map PROPERTIES FOLDER third_party)
add_library(slot_map::slot_map ALIAS slot_map)
