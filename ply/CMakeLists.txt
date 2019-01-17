cmake_minimum_required (VERSION 3.6)
project (ply)

set(PLY_SOURCES
 ply.c
 ply.h
 )
 
include_directories(${CMAKE_SOURCE_DIR})
add_library( ply ${PLY_SOURCES})

SOURCE_GROUP( Source\ Files FILES ${PLY_SOURCES} )
