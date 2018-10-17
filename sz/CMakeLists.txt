aux_source_directory(src sz_sources)

list (APPEND sz_dependencies zlib zstd)

add_library (SZ SHARED ${sz_sources})
target_link_libraries (SZ ${sz_dependencies})

install (TARGETS SZ LIBRARY DESTINATION lib)
