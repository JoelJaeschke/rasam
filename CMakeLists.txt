cmake_minimum_required(VERSION 3.15...3.25)

project(rasterio-parsam
  VERSION
    0.0.1
  DESCRIPTION
    "Parallel raster sampling using GDAL"
  LANGUAGES
    CXX
)

#find_package(GDAL CONFIG REQUIRED)
find_package(GDAL REQUIRED)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)