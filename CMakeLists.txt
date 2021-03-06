cmake_minimum_required (VERSION 2.6)
project (SymUniverse)

# Basic project definitions
set (SymUniverse_PROJECT_NAME "SymUniverse")
set (SymUniverse_PROJECT_AUTHOR "J. Lowell Wofford")
set (SymUniverse_VERSION_MAJOR 1)
set (SymUniverse_VERSION_MINOR 0)

set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
if( APPLE )
	add_definitions( -DAPPLE )
elseif( CMAKE_SYSTEM_NAME STREQUAL "Linux" )
	add_definitions( -DLINUX )
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
else( APPLE )
	message( FATAL_ERROR "System type ${CMAKE_SYSTEM_NAME} is not
supported!" )
endif( )

# Setup a compiled .h file
configure_file (
	"${PROJECT_SOURCE_DIR}/include/SymUniverseConfig.h.in"
	"${PROJECT_BINARY_DIR}/include/SymUniverseConfig.h"
)

# Setup includes
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB_RECURSE INCLUDES "include/*.h" "include/*.h.in")

# Define LIBRARIES varaible, set here because other builds need it.
set(LIBRARIES universe boundaries)

# Add subdirectories
add_subdirectory(src)

add_custom_target(docs ALL SOURCES "LICENSE" "README.md")

# // TODO: Setup install directives
#install(TARGETS sym ubuild umerge ucat uextract udivide DESTINATION bin)
#install(FILES universe.h DESTINATION include)
#install(FILES dummy.mod DESTINATION modules)
#install(FILES "${PROJECT_BINARY_DIR}/SymUniverseConfig.h" DESTINATION include)

# // TODO: Setup installer
# build a CPack driven installer package
#include (InstallRequiredSystemLibraries)
#set (CPACK_RESOURCE_FILE_LICENSE
#	"${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
#set (CPACK_PACKAGE_VERSION_MAJOR "${SymUniverse_VERSION_MAJOR}")
#set (CPACK_PACKAGE_VERSION_MINOR "${SymUniverse_VERSION_MINOR}")
#include(CPack)
