if(TARGET nonstd::span-lite)
    return()
endif()

message(STATUS "Third-party (external): creating target 'nonstd::span-lite'")

include(CPM)
CPMAddPackage(
    NAME span-lite
    GITHUB_REPOSITORY martinmoene/span-lite
    GIT_TAG bc08bf87258d881aaa83b50c54dea67ea33d0e8e
)

set_target_properties(span-lite PROPERTIES FOLDER third_party)
