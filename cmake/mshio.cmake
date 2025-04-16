if (TARGET mshio::mshio)
    return()
endif()

message(STATUS "Third-party (external): creating target 'mshio::mshio'")

include(CPM)
CPMAddPackage(
  NAME mshio
  GITHUB_REPOSITORY qnzhou/MshIO
  GIT_TAG 003ee572f9c4b3ac36f546501a2e157e6f47a4fd
)

set_target_properties(mshio PROPERTIES
    FOLDER third_party
    POSITION_INDEPENDENT_CODE On
)
