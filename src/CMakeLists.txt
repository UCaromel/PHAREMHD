cmake_minimum_required(VERSION 3.10)
project(MHDProblem)

set(SOURCES
    ConservativeVariablesCC.cpp
    PrimitiveVariablesCC.cpp
    ConstainedTransport.hpp
    GodunovFlux.hpp
    Initialisation.hpp
    Interface.hpp
    ReconstructedValues.hpp
    RusanovRiemannSolver.hpp
    TimeIntegrator.hpp
    EquationOfState.cpp
    AddGhostCells.hpp
)

foreach(SOURCE_FILE ${SOURCES})
    get_filename_component(FILE_EXT ${SOURCE_FILE} EXT)
    if(${FILE_EXT} STREQUAL ".cpp")
        list(APPEND CPP_SOURCES ${SOURCE_FILE})
    endif()
endforeach()

# Create object libraries for each source file
foreach(SOURCE_FILE ${SOURCES})
    get_filename_component(FILE_NAME ${SOURCE_FILE} NAME_WE)
    add_library(${FILE_NAME} OBJECT ${SOURCE_FILE})

    # Check if the file is a .cpp file and include directories if needed
    list(FIND CPP_SOURCES ${SOURCE_FILE} IS_CPP)
    if(IS_CPP GREATER -1)
        target_include_directories(${FILE_NAME} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    endif()

    # Collect all object libraries
    list(APPEND OBJECT_LIBRARIES ${FILE_NAME})
endforeach()

# Add the executable
add_executable(MHDProblem MHDProblem.cpp)

# Link all object libraries to the executable
target_link_libraries(MHDProblem ${OBJECT_LIBRARIES})

# Include directories if required
target_include_directories(MHDProblem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
