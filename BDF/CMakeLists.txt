project ( bdf2.1 )
cmake_minimum_required ( VERSION 2.8 )

find_package ( Boost REQUIRED )

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" ) 

include_directories ( ${Boost_INCLUDE_DIRS} )

add_subdirectory (examples/system1)
add_subdirectory (examples/system2)
add_subdirectory (examples/system3)

