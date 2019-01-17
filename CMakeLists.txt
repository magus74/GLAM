cmake_minimum_required (VERSION 3.6)
project (glam)

set(EXTERNAL_DIR "../../../external")

# BOOST
 set(BOOST_LIBRARYDIR "C:/Boost/lib")
 set(Boost_ADDITIONAL_VERSIONS "1.66" "1.66.0")
# set(Boost_DEBUG ON)
# set(Boost_USE_STATIC_LIBS       ON)
# set(Boost_USE_MULTITHREADED      ON)
# set(Boost_USE_STATIC_RUNTIME    ON)
	
 find_package( Boost 1.50 REQUIRED COMPONENTS system thread program_options )
 include_directories ( ${Boost_INCLUDE_DIRS} )
 link_directories ( ${Boost_LIBRARY_DIRS} )

#
# GLAM - Glycogen Lactate Absorption Modelling
#
set(GLAM_CORE_SOURCES
 main.cpp
 scalar_field.cpp
 parser.cpp
 influence_map_builder.cpp
 cluster_analysis.cpp
 scalar_field.hpp
 parser.hpp
 influence_map_builder.hpp
 cluster_analysis.hpp
 )

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${EXTERNAL_DIR}/include)
link_directories(${EXTERNAL_DIR}/lib)

add_executable (glam ${GLAM_CORE_SOURCES})

target_link_libraries (glam
	sl
	ply
	ply_io)

SOURCE_GROUP( Source\ Files FILES ${GLAM_CORE_SOURCES} )
