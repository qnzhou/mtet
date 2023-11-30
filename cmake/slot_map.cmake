if(TARGET slot_map::slot_map)
    return()
endif()

message(STATUS "Third-party (external): creating target 'slot_map::slot_map'")

include(CPM)
CPMAddPackage(
  NAME slot_map
  GITHUB_REPOSITORY qnzhou/slot_map
  GIT_TAG bc21042de3b20dd96a7d42e454669350b2b74388
  DOWNLOAD_ONLY YES
)

add_subdirectory(${slot_map_SOURCE_DIR}/slot_map)
set_target_properties(slot_map PROPERTIES SYSTEM ON)
set_target_properties(slot_map PROPERTIES FOLDER third_party)
add_library(slot_map::slot_map ALIAS slot_map)
