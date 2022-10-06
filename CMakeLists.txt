cmake_minimum_required(VERSION 3.24)

project(delaunaypan3d_test CXX)
set (CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.12.0.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_executable(DelaunayPan3DTest test/mesher_test.cpp)
target_link_libraries(DelaunayPan3DTest gtest_main)

include(GoogleTest)
gtest_discover_tests(DelaunayPan3DTest)
