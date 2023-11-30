cmake_minimum_required(VERSION 3.14)

project(mtet)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(indirect-value-lite)
include(slot_map)
include(unordered_dense)
include(nanothread)
include(sanitizers)

option(MTET_BUILD_TESTS "Build tests" OFF)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(mtet STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(mtet PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_compile_features(mtet PRIVATE cxx_std_17)
target_link_libraries(mtet PRIVATE
    slot_map::slot_map
    unordered_dense::unordered_dense
    nanothread::nanothread
    PUBLIC
    nonstd::indirect-value-lite
)
add_library(mtet::mtet ALIAS mtet)
target_compile_features(mtet PRIVATE cxx_std_20)

if(MTET_BUILD_TESTS)
    include(CTest)
    enable_testing()
    include(catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
    add_executable(mtet_tests ${TEST_FILES})
    target_include_directories(mtet_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR}/src)
    target_link_libraries(mtet_tests PRIVATE
        mtet::mtet
        slot_map::slot_map
        unordered_dense::unordered_dense
        nanothread::nanothread
        Catch2::Catch2WithMain)
    target_compile_definitions(mtet_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
    target_compile_features(mtet_tests PRIVATE cxx_std_20)

    catch_discover_tests(mtet_tests)
endif()
