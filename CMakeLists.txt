cmake_minimum_required (VERSION 2.6)
project (SymUniverse)

set (SymUniverse_PROJECT_NAME "SymUniverse")
set (SymUniverse_PROJECT_AUTHOR "J. Lowell Wofford")
set (SymUniverse_VERSION_MAJOR 1)
set (SymUniverse_VERSION_MINOR 0)

configure_file (
	"${PROJECT_SOURCE_DIR}/include/SymUniverseConfig.h.in"
	"${PROJECT_BINARY_DIR}/include/SymUniverseConfig.h"
)

include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(universe src/universe.c)

add_executable(sym src/sym.c)
add_executable(ubuild src/ubuild.c)
#add_executable(umerge src/umerge.c)
#add_executable(ucat src/ucat.c)
#add_executable(uextract src/uextract.c)
#add_executable(udivide src/udivide.c)
target_link_libraries(sym universe)
target_link_libraries(ubuild universe)
#target_link_libraries(umerge universe)
#target_link_libraries(ucat universe)
#target_link_libraries(uextract universe)
#target_link_libraries(udivide universe)

add_library(dummy MODULE src/mod_dummy.c)
target_link_libraries(dummy universe)

if(APPLE)
    set_property(TARGET dummy PROPERTY PREFIX "")
    set_property(TARGET dummy PROPERTY SUFFIX ".mod")
endif()

#install(TARGETS sym ubuild umerge ucat uextract udivide DESTINATION bin)
install(TARGETS sym ubuild DESTINATION bin)
install(FILES universe.h DESTINATION include)
install(FILES dummy.mod DESTINATION modules)
install(FILES "${PROJECT_BINARY_DIR}/SymUniverseConfig.h" DESTINATION include)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
	"${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${SymUniverse_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${SymUniverse_VERSION_MINOR}")
include(CPack)
