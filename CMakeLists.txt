cmake_minimum_required(VERSION 3.0.0)
project(gradient_proect VERSION 0.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


add_subdirectory(matrix)
add_executable(gradient_proect gradient_proect.cpp)


target_link_libraries(gradient_proect PUBLIC matrix)


target_include_directories(gradient_proect PUBLIC
                          "${PROJECT_BINARY_DIR}/matrix"
                          "${PROJECT_SOURCE_DIR}/matrix"
                          )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
