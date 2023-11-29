if(TARGET nonstd::value-ptr-lite)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nonstd::value-ptr-lite'")

include(CPM)
CPMAddPackage(
    NAME value-ptr-lite
    GITHUB_REPOSITORY martinmoene/value-ptr-lite
    GIT_TAG v0.2.1
)

set_target_properties(value_ptr-lite PROPERTIES FOLDER third_party)
