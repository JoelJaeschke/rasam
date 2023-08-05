cmake_minimum_required(VERSION 3.15...3.25)

project(parsam
  VERSION
    0.0.1
  DESCRIPTION
    "Parallel raster sampling using GDAL"
  LANGUAGES
    CXX
)

SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
find_package(GDAL REQUIRED)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)