cmake_minimum_required (VERSION 3.1.3)
project (sz)
enable_testing()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=c99 -lm -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -lm -fopenmp")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#find_package (FFTW)
if (FFTW_FOUND)
  include_directories (${FFTW_INCLUDES})
  add_definitions (-DHAVE_FFTW3)
endif ()

#find_package (netCDF)
if (netCDF_FOUND)
  include_directories (${NETCDF_INCLUDE_DIRS})
endif ()

find_program(TAGS ctags)
if(TAGS)
	add_custom_target(tags ALL
		COMMAND ${TAGS} -f ${CMAKE_BINARY_DIR}/tags -R
		COMMENT Generating Tag files
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()

find_package (ZLIB)
if (ZLIB_FOUND)
  include_directories (${ZLIB_INCLUDE_DIRS})
else ()
  set (ZLIB_LIBRARIES zlib)
  add_subdirectory (zlib)
  include_directories (${CMAKE_CURRENT_SOURCE_DIR}/zlib)
endif ()

find_package (ZSTD)
if (ZSTD_FOUND)
  include_directories (${ZSTD_INCLUDE_DIRS})
else ()
  set (ZSTD_LIBRARIES zstd)
  add_subdirectory (zstd)
  include_directories (${CMAKE_CURRENT_SOURCE_DIR}/zstd)
endif ()

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/sz/include)

add_subdirectory (sz)
add_subdirectory (example)
#add_subdirectory (test)

