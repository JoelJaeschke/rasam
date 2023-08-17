add_executable(rasam rasam.cpp)
target_link_libraries(rasam PUBLIC GDAL::GDAL)

if(MSVC)
  target_compile_options(rasam PRIVATE /W4 /WX)
else()
  target_compile_options(rasam PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()