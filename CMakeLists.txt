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


file(GLOB_RECURSE INCS "include/*.h" "include/*.h.in")

add_library(universe src/universe.c ${INCS})

add_executable(sym src/sym.c ${INCS})
add_executable(ubuild src/ubuild.c ${INCS})
add_executable(umerge src/umerge.c ${INCS})
add_executable(ucat src/ucat.c ${INCS})
add_executable(uextract src/uextract.c ${INCS})
add_executable(udivide src/udivide.c ${INCS})
target_link_libraries(sym universe)
target_link_libraries(ubuild universe)
target_link_libraries(umerge universe)
target_link_libraries(ucat universe)
target_link_libraries(uextract universe)
target_link_libraries(udivide universe)

add_library(dummy MODULE src/mod_dummy.c ${INCS})
target_link_libraries(dummy universe)
add_library(integrate MODULE src/mod_integrate.c ${INCS})
target_link_libraries(integrate universe)

if(APPLE)
    set_property(TARGET dummy PROPERTY PREFIX "")
    set_property(TARGET dummy PROPERTY SUFFIX ".mod")
    set_property(TARGET integrate PROPERTY PREFIX "")
    set_property(TARGET integrate PROPERTY SUFFIX ".mod")
endif()

install(TARGETS sym ubuild umerge ucat uextract udivide DESTINATION bin)
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
