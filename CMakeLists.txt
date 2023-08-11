cmake_minimum_required(VERSION 3.22)

project(rasam
  VERSION
    0.0.1
  DESCRIPTION
    "Parallel raster sampling using GDAL"
  LANGUAGES
    CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(GDAL REQUIRED)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)