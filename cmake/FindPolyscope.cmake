# - Try to find the Polyscope library
# Once done this will define
#
#  POLYSCOPE_FOUND - system has POLYSCOPE 
#  POLYSCOPE_INCLUDE_DIR - **the** POLYSCOPE include directory
if(POLYSCOPE_FOUND)
    return()
endif()

find_path(POLYSCOPE_INCLUDE_DIR polyscope/polyscope.h
    HINTS
        ENV POLYSCOPE_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/polyscope
        ${CMAKE_SOURCE_DIR}/../polyscope
        ${CMAKE_SOURCE_DIR}/../tools/polyscope
        ${CMAKE_SOURCE_DIR}/../../polyscope
        ${CMAKE_SOURCE_DIR}/lib/polyscope
        /usr
        /usr/local
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Polyscope 
    "\npolyscope not found --- You can download it using:\n\tgit clone --recursive https://github.com/nmwsharp/polyscope.git ${CMAKE_SOURCE_DIR}/libs/polyscope"
    POLYSCOPE_INCLUDE_DIR)
mark_as_advanced(POLYSCOPE_INCLUDE_DIR)

add_subdirectory("${POLYSCOPE_INCLUDE_DIR}/../" "polyscope")