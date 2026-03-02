cmake_minimum_required(VERSION 3.20.0)
project (CADTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17) 

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(samples)
