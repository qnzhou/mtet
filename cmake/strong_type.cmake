if(TARGET strong_type::strong_type)
    return()
endif()

message(STATUS "Third-party (external): creating target 'strong_type::strong_type'")

include(CPM)
CPMAddPackage(
    NAME strong_type
    GITHUB_REPOSITORY rollbear/strong_type
    GIT_TAG v13
)

set_target_properties(strong_type PROPERTIES SYSTEM ON)
set_target_properties(strong_type PROPERTIES FOLDER third_party)
