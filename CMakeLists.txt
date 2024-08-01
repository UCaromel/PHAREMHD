cmake_minimum_required (VERSION 3.20.1) # released April 8, 2021 - https://www.kitware.com/cmake-3-20-1-available-for-download/

project(PHAREMHD VERSION 0.1 LANGUAGES CXX C)

# Release mode default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()


SET(PYBIND_MIN_VERSION "2.5.0")
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
find_program(Git git)
include(pybind.cmake)

set(PHAREMHD_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(SOURCES
    src/ConservativeVariables.cpp 
    src/ConstrainedTransport.cpp
    src/EquationOfState.cpp
    src/GodunovFlux.cpp
    src/Interface.cpp
    src/ModularityUtils.cpp
    src/PhareMHD.cpp
    src/PhysicalConstants.cpp
    src/PrimitiveVariables.cpp
    src/RiemannSolver.cpp
    src/TimeIntegrator.cpp)



pybind11_add_module(pyMHD pyMHD/pyMHD.cpp
    src/ConservativeVariables.cpp 
    src/ConstrainedTransport.cpp
    src/EquationOfState.cpp
    src/GodunovFlux.cpp
    src/Interface.cpp
    src/ModularityUtils.cpp
    src/PhareMHD.cpp
    src/PhysicalConstants.cpp
    src/PrimitiveVariables.cpp
    src/RiemannSolver.cpp
    src/TimeIntegrator.cpp)





