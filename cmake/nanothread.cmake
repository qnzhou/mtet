if(TARGET nanothread::nanothread)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nanothread::nanothread'")

include(CPM)
CPMAddPackage(
  NAME nanothread
  GITHUB_REPOSITORY mitsuba-renderer/nanothread
  GIT_TAG 94b3237e4c777f044a5118b10c321af5be57e88e
)

set_target_properties(nanothread PROPERTIES FOLDER third_party)
add_library(nanothread::nanothread ALIAS nanothread)
