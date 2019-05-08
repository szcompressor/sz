aux_source_directory(. zlib_sources)

add_library (ZLIB SHARED ${zlib_sources})
target_link_libraries (ZLIB ${zlib_dependencies})
target_include_directories(ZLIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIRECTORY})

install (TARGETS ZLIB LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
