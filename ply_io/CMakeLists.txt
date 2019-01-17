cmake_minimum_required (VERSION 3.6)
project (ply_io)

set(EXTERNAL_DIR "../external")
include_directories(${EXTERNAL_DIR}/include)
link_directories(${EXTERNAL_DIR}/lib)

set(PLY_IO_SOURCES
ooc_ply_io.cpp 
ooc_triangle_mesh.cpp
ooc_ply_io.hpp 
ooc_triangle_mesh.hpp
std_ply_io.cpp 
std_triangle_mesh.cpp
std_ply_io.hpp 
std_triangle_mesh.hpp
 )


include_directories(${CMAKE_SOURCE_DIR})

add_library( ply_io ${PLY_IO_SOURCES})
target_link_libraries (ply_io
	sl
	ply)
	
SOURCE_GROUP( Source\ Files FILES ${PLY_SOURCES} )
