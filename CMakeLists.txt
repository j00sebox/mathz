cmake_minimum_required(VERSION 3.2)

project(mathz)

add_subdirectory(external/catch2)

add_executable(${PROJECT_NAME} tests/test.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
