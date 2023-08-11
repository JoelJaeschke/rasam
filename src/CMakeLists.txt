add_executable(rasam rasam.cpp)
target_link_libraries(rasam PUBLIC GDAL::GDAL)

if(MSVC)
  target_compile_options(rarsam PRIVATE /W4 /WX)
else()
  target_compile_options(rarsam PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()