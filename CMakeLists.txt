cmake_minimum_required(VERSION 3.1)
project(discrete_elastic_rods)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Find Eigen
find_package(Eigen3 REQUIRED)

# Find Polyscope
find_package(Polyscope REQUIRED)
  
file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} polyscope)