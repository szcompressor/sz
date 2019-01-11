aux_source_directory(src sz_sources)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/include/sz_opencl_config.h.in
	${CMAKE_CURRENT_BINARY_DIR}/include/sz_opencl_config.h @ONLY)

add_library (SZ SHARED ${sz_sources})
target_include_directories(SZ PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries (SZ zlib zstd m ${OpenCL_LIBRARIES})
target_compile_options(SZ
	PUBLIC -DHAVE_RANDOMACCESS
	PRIVATE $<$<CONFIG:Debug>:-Wall -Wextra -Wpedantic -Wno-unused-parameter>
	)

install (TARGETS SZ LIBRARY DESTINATION lib)

