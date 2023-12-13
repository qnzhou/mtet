cmake_minimum_required(VERSION 3.14)

project(mtet)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(indirect-value-lite)
include(mshio)
include(slot_map)
include(unordered_dense)
include(nanothread)
include(sanitizers)
include(smallvector)
include(strong_type)

option(MTET_BUILD_TESTS "Build tests" OFF)
option(MTET_BUILD_EXAMPLES "Build examples" OFF)

file(GLOB INC_FILES ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(mtet STATIC ${SRC_FILES} ${INC_FILES})
target_include_directories(mtet PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_compile_features(mtet PRIVATE cxx_std_17)
target_link_libraries(mtet PRIVATE
    slot_map::slot_map
    unordered_dense::unordered_dense
    nanothread::nanothread
    mshio::mshio
    smallvector::smallvector
    PUBLIC
    nonstd::indirect-value-lite
    strong_type::strong_type
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
        smallvector::smallvector
        nanothread::nanothread
        Catch2::Catch2WithMain)
    target_compile_definitions(mtet_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
    target_compile_features(mtet_tests PRIVATE cxx_std_20)

    if (WIN32)
        add_custom_command(
            TARGET mtet_tests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE_DIR:nanothread>/nanothread.dll
            $<TARGET_FILE_DIR:mtet_tests>)
    endif()

    catch_discover_tests(mtet_tests)
endif()

if (MTET_BUILD_EXAMPLES)
    add_executable(longest_edge_split ${CMAKE_CURRENT_LIST_DIR}/examples/longest_edge_split.cpp)
    target_link_libraries(longest_edge_split PRIVATE mtet::mtet)
    target_compile_features(longest_edge_split PRIVATE cxx_std_20)

    if (WIN32)
        add_custom_command(
            TARGET longest_edge_split POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE_DIR:nanothread>/nanothread.dll
            $<TARGET_FILE_DIR:longest_edge_split>)
    endif()
endif()
