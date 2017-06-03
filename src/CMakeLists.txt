set(CMAKE_CXX_STANDARD 11)
find_package(Threads REQUIRED)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(Equestrotopia SHARED input.cpp geometry.cpp brdf.cpp kdtree.cpp light.cpp texture.cpp)

ADD_EXECUTABLE(raytracing raytracing.cpp)
target_link_libraries(raytracing Equestrotopia)

ADD_EXECUTABLE(photontracing photontracing.cpp)
target_link_libraries(photontracing Equestrotopia)

ADD_EXECUTABLE(updation updation.cpp)
target_link_libraries(updation Equestrotopia ${CMAKE_THREAD_LIBS_INIT})

ADD_EXECUTABLE(place place.cpp)
target_link_libraries(place Equestrotopia)
