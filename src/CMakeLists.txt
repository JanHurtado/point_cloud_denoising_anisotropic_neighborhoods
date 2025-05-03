project(PointCloudDenoising)
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0020 NEW)

set(CMAKE_CXX_STANDARD 14)

add_definitions(-D_USE_MATH_DEFINES)

find_package(OpenMP)

if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include "${PROJECT_SOURCE_DIR}/../dependencies/OpenMesh/include" "${PROJECT_SOURCE_DIR}/../dependencies/Eigen/include" "${PROJECT_SOURCE_DIR}/../dependencies/cplex/include" "${PROJECT_SOURCE_DIR}/../dependencies/geogram/include/geogram1")

link_directories("${PROJECT_SOURCE_DIR}/../dependencies/OpenMesh/lib" "${PROJECT_SOURCE_DIR}/../dependencies/cplex/lib" "${PROJECT_SOURCE_DIR}/../dependencies/Geogram/lib")

file(GLOB_RECURSE PCDENOISING_SOURCES *.cpp *.h)

add_executable(pcdenoising ${PCDENOISING_SOURCES})

target_include_directories(pcdenoising
    PUBLIC ${INCLUDE_DIRS}
)

target_link_libraries(
    pcdenoising
	debug OpenMeshCored
	optimized OpenMeshCore
	optimized cplex2010
	debug cplex2010d
	geogram
)

install( TARGETS pcdenoising DESTINATION bin )

file( GLOB CPLEX_DLLS ../dependencies/cplex/bin/*.dll )
file( GLOB GEOGRAM_DLLS ../dependencies/geogram/bin/*.dll )

install( FILES ${CPLEX_DLLS} DESTINATION bin )
install( FILES ${GEOGRAM_DLLS} DESTINATION bin )