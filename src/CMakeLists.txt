# add_executable(main main.cpp)
# target_link_libraries(main PRIVATE GDAL::GDAL)

add_executable(parsam parsam.cpp)
target_link_libraries(parsam PRIVATE GDAL::GDAL)