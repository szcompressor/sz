aux_source_directory(. zlib_sources)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

add_library (zlib SHARED ${zlib_sources})
target_link_libraries (zlib ${zlib_dependencies})

install (TARGETS zlib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

set (zlib_headers
      crc32.h  deflate.h  gzguts.h  inffast.h  inffixed.h  inflate.h  inftrees.h  trees.h  zconf.h  zlib.h  zutil.h  )

install (FILES ${zlib_headers} DESTINATION include)
