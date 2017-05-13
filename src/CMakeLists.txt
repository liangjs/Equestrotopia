set(CMAKE_CXX_STANDARD 11)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
ADD_EXECUTABLE(Equestrotopia main.cpp input.cpp geometry.cpp brdf.cpp kdtree.cpp light.cpp)
target_link_libraries(Equestrotopia ${OPENGL_LIBRARY} ${GLUT_LIBRARIES})
