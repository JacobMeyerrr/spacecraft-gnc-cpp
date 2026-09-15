cmake_minimum_required(VERSION 3.20)

project(spacecraft_gnc_cpp
    VERSION 1.0.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

# ------------------------------------------------------------------
# Reusable spacecraft GNC library
#
# The validation programs below are clients of this library. The
# library itself is not organized around the homework assignments.
# ------------------------------------------------------------------

add_library(spacecraft_gnc
    src/attitude/AttitudeKinematics.cpp
    src/attitude/EulerAngles.cpp
    src/attitude/Quaternion.cpp
    src/attitude/RotationMatrix.cpp

    src/control/ControlSystem.cpp
    src/control/PIDController.cpp
    src/control/StateSpace.cpp

    src/dynamics/GravityGradient.cpp
    src/dynamics/RigidBodyDynamics.cpp
    src/dynamics/SpacecraftProperties.cpp

    src/math/Polynomial.cpp

    src/models/EarthPointingSpacecraft.cpp
    src/models/MomentumBiasedSpacecraft.cpp

    src/orbital/Orbit.cpp
    src/orbital/OrbitPropagator.cpp
    src/orbital/ReferenceFrames.cpp

    src/simulation/DataLogger.cpp
    src/simulation/Simulator.cpp
)

target_include_directories(spacecraft_gnc
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(spacecraft_gnc
    PUBLIC
        Eigen3::Eigen
)

target_compile_options(spacecraft_gnc PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

# ------------------------------------------------------------------
# Validation: HW1
# ------------------------------------------------------------------

add_executable(validate_hw1_orbit
    validation/hw1/problem1_orbit.cpp
)
target_link_libraries(validate_hw1_orbit PRIVATE spacecraft_gnc)

add_executable(validate_hw1_attitude
    validation/hw1/problem2_euler.cpp
)
target_link_libraries(validate_hw1_attitude PRIVATE spacecraft_gnc)

# ------------------------------------------------------------------
# Validation: HW2
# ------------------------------------------------------------------

add_executable(validate_hw2_quaternions
    validation/hw2/problem1_quaternions.cpp
)
target_link_libraries(validate_hw2_quaternions PRIVATE spacecraft_gnc)

add_executable(validate_hw2_attitude
    validation/hw2/problem2_attitude.cpp
)
target_link_libraries(validate_hw2_attitude PRIVATE spacecraft_gnc)

add_executable(validate_hw2_pointing
    validation/hw2/problem3_pointing.cpp
)
target_link_libraries(validate_hw2_pointing PRIVATE spacecraft_gnc)

# ------------------------------------------------------------------
# Validation: HW3
# ------------------------------------------------------------------

add_executable(validate_hw3
    validation/hw3/main.cpp
)
target_link_libraries(validate_hw3 PRIVATE spacecraft_gnc)

# ------------------------------------------------------------------
# Validation: HW4
# ------------------------------------------------------------------

add_executable(validate_hw4
    validation/hw4/main.cpp
)
target_link_libraries(validate_hw4 PRIVATE spacecraft_gnc)
