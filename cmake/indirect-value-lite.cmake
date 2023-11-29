if(TARGET nonstd::indirect-value-lite)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nonstd::indirect-value-lite'")

include(CPM)
CPMAddPackage(
    NAME indirect-value-lite
    GITHUB_REPOSITORY martinmoene/indirect-value-lite
    GIT_TAG v0.1.0
)

set_target_properties(indirect-value-lite PROPERTIES FOLDER third_party)
