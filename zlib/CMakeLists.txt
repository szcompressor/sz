aux_source_directory(. zlib_sources)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

add_library (zlib SHARED ${zlib_sources})
target_link_libraries (zlib ${zlib_dependencies})

install (TARGETS zlib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
