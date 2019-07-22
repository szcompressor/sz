cmake_minimum_required (VERSION 3.14)
project (SZ VERSION 2.1.5)
enable_testing()

include(GNUInstallDirs)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

find_package (FFTW)
if (FFTW_FOUND)
  include_directories (${FFTW_INCLUDES})
  add_definitions (-DHAVE_FFTW3)
endif ()

find_package (NETCDF)
if (NETCDF_FOUND)
  include_directories (${NETCDF_INCLUDE_DIRS})
endif ()

find_package(OpenCL)
if(OpenCL_FOUND)
	include_directories(${OpenCL_INCLUDE_DIRS})
	add_definitions(-DHAVE_OPENCL)
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_search_module(ZLIB IMPORTED_TARGET GLOBAL zlib)
  if(NOT ZLIB_FOUND)
    set(USE_PACKAGED_ZLIB 1)
  else()
    add_library(ZLIB ALIAS PkgConfig::ZLIB)
  endif()

  pkg_search_module(ZSTD IMPORTED_TARGET GLOBAL libzstd)
  if(NOT ZSTD_FOUND)
    set(USE_PACKAGED_ZSTD 1)
  else()
    add_library(ZSTD ALIAS PkgConfig::ZSTD)
  endif()

else()
  set(USE_PACKAGED_ZLIB 1)
  set(USE_PACKAGED_ZSTD 1)
endif()

if(USE_PACKAGED_ZSTD)
  add_subdirectory(zstd)
endif()


if(USE_PACKAGED_ZLIB)
  add_subdirectory(zlib)
endif()

find_program(TAGS ctags)
if(TAGS)
	add_custom_target(tags ALL
    COMMAND ${TAGS} --exclude=${CMAKE_BINARY_DIR} -f ${CMAKE_BINARY_DIR}/tags --c++-kinds=+p --fields=+iaS --extras=+q -R 
		COMMENT Generating Tag files
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()


add_subdirectory (sz)
add_subdirectory (example)
add_subdirectory (test)
add_subdirectory (swig)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/sz.pc.in
  ${CMAKE_BINARY_DIR}/sz.pc
  @ONLY
)
install(FILES ${CMAKE_BINARY_DIR}/sz.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pkgconfig)
