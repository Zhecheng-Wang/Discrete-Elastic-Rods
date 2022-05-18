cmake_minimum_required(VERSION 3.16)
project(deram)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "DEBUG")
endif()

# Find libigl
include(libigl)

# Find Polyscope
find_package(Polyscope REQUIRED)

# Find protobuf
include(FindProtobuf)
find_package(Protobuf REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp src/rod.pb.h src/rod.pb.cc .src/proto/*.cc .src/proto/*.h .src/proto/*.proto)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories( ${PROJECT_NAME} PUBLIC
						    ${EIGEN3_INCLUDE_DIR}
							${LIBIGL_INCLUDE_DIR}
							${CMAKE_CURRENT_BINARY_DIR}/proto
							./src
							./src/proto
							)

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} igl::core polyscope protobuf::libprotoc protobuf::libprotobuf protobuf::libprotobuf-lite)
