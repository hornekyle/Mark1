cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
project("Mark1")

set(source
	"src/integrate.f90"
	"src/kinds.f90"
	"src/main.f90"
	"src/output.f90"
	#"src/plplotlib.f90"
	"src/system.f90"
	"src/settings.f90"
	"src/units.f90"
	"src/utilities.f90"
	"src/lmpIntegrator.f90"
	)

set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan -fbackslash")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mtune=native -march=native")

#include(FindPkgConfig)
#pkg_search_module(PLPLOT REQUIRED plplot-f95 plplotd-f95)
#include_directories(${PLPLOT_INCLUDE_DIRS})
#link_directories(${PLPLOT_LIBRARY_DIRS})

add_custom_target(symlink
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/scripts" "${CMAKE_BINARY_DIR}/scripts"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/lammps" "${CMAKE_BINARY_DIR}/lammps"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

add_executable(main ${source})
#target_link_libraries(main ${PLPLOT_LIBRARIES})
#add_dependencies(main symlink)

set(FORD_EXECUTABLE "ford")
add_custom_target(ford
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/pages" "${CMAKE_BINARY_DIR}/pages"
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}/src"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/project.md" "${CMAKE_BINARY_DIR}"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/README.md" "${CMAKE_BINARY_DIR}"
	COMMAND ${FORD_EXECUTABLE} "project.md"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
