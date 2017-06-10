set(CMAKE_CXX_STANDARD 11)
find_package(Threads REQUIRED)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(Equestrotopia SHARED input.cpp geometry.cpp brdf.cpp kdtree.cpp light.cpp texture.cpp ppm.cpp)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/lodepng)
add_library(lodepng SHARED lodepng/lodepng.cpp)

ADD_EXECUTABLE(raytracing raytracing.cpp)
target_link_libraries(raytracing Equestrotopia lodepng)

ADD_EXECUTABLE(photontracing photontracing.cpp)
target_link_libraries(photontracing Equestrotopia)

ADD_EXECUTABLE(updation updation.cpp)
target_link_libraries(updation Equestrotopia lodepng ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(updation PUBLIC lodepng)

ADD_EXECUTABLE(place place.cpp)
target_link_libraries(place Equestrotopia)
