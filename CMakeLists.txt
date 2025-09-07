cmake_minimum_required(VERSION 3.8)
project(path_planning)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Find Dependencies
find_package(Eigen3 REQUIRED)

add_executable(cubic_spline
    src/cubic_spline.cpp 
    include/cubic_spline.h 
    include/solve_tridiagonal.h
)

add_executable(b_spline
    src/b-spline.cpp 
    include/b-spline.h 
    include/solve_tridiagonal.h
)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
)