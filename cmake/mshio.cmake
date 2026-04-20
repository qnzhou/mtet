if (TARGET mshio::mshio)
    return()
endif()

message(STATUS "Third-party (external): creating target 'mshio::mshio'")

include(CPM)
CPMAddPackage(
  NAME mshio
  GITHUB_REPOSITORY qnzhou/MshIO
  GIT_TAG v0.1.1
)

set_target_properties(mshio PROPERTIES
    FOLDER third_party
    POSITION_INDEPENDENT_CODE On
)
