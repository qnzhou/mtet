if(TARGET nanothread::nanothread)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nanothread::nanothread'")

include(CPM)
CPMAddPackage(
  NAME nanothread
  GITHUB_REPOSITORY mitsuba-renderer/nanothread
  GIT_TAG 9073b959f02da3395cdae8ed7e0e0f86b1c2ddb8
)

set_target_properties(nanothread PROPERTIES FOLDER third_party)
add_library(nanothread::nanothread ALIAS nanothread)
