if(TARGET nanothread::nanothread)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nanothread::nanothread'")

include(CPM)
CPMAddPackage(
  NAME nanothread
  GITHUB_REPOSITORY mitsuba-renderer/nanothread
  GIT_TAG e76778dca2bad7e399cc9e408e371dc73569e125
)

set_target_properties(nanothread PROPERTIES FOLDER third_party)
add_library(nanothread::nanothread ALIAS nanothread)
