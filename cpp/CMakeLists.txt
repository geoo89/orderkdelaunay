# CMake for order-k Delaunay triangulations using CGAL.
# Includes a commandline tool (orderk) and tests (test).

project (orderk)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_definitions(-std=c++11)

find_package(CGAL COMPONENTS)
include(${CGAL_USE_FILE}) 

if ( CGAL_FOUND )

set(DEPENDENCIES  ) # no dependencies at this point

add_executable(orderk src/main.cpp ${DEPENDENCIES})
add_executable(tests src/tests.cpp ${DEPENDENCIES})
add_to_cached_list( CGAL_EXECUTABLE_TARGETS orderk )
target_link_libraries( orderk ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL, and will not be compiled.")

endif()