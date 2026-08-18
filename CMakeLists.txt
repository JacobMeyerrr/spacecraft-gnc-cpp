cmake_minimum_required(VERSION 3.10)
project(SpacecraftGNC VERSION 1.0 LANGUAGES CXX)

# Use modern C++17 standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Tell CMake where to look for your .hpp header files
include_directories(include)

# Include the source code directory
add_subdirectory(src)

# ---------------------------------------------------
# HOMEWORK EXECUTABLES (Your runnable programs)
# ---------------------------------------------------

# Homework 1: Orbit and Euler Problems
add_executable(hw1_orbit homework/hw1/problem1_orbit.cpp)
add_executable(hw1_euler homework/hw1/problem2_euler.cpp)

# Homework 2: Quaternions, Attitude, and Pointing
add_executable(hw2_quaternions homework/hw2/problem1_quaternions.cpp)
add_executable(hw2_attitude homework/hw2/problem2_attitude.cpp)
add_executable(hw2_pointing homework/hw2/problem3_pointing.cpp)


