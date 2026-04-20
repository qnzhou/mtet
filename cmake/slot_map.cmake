if(TARGET slot_map::slot_map)
    return()
endif()

message(STATUS "Third-party (external): creating target 'slot_map::slot_map'")

include(CPM)
CPMAddPackage(
  NAME slot_map
  GITHUB_REPOSITORY SergeyMakeev/SlotMap
  GIT_TAG 122fb6960bc417f82f65a0469d38bfa05b506a25
  DOWNLOAD_ONLY YES
)

add_subdirectory(${slot_map_SOURCE_DIR}/slot_map ${slot_map_BINARY_DIR})
set_target_properties(slot_map PROPERTIES SYSTEM ON)
set_target_properties(slot_map PROPERTIES FOLDER third_party)
add_library(slot_map::slot_map ALIAS slot_map)
