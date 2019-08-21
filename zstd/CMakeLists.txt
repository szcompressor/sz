include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}/
  ${CMAKE_CURRENT_SOURCE_DIR}/common
  ${CMAKE_CURRENT_SOURCE_DIR}/compress
  ${CMAKE_CURRENT_SOURCE_DIR}/decompress
  ${CMAKE_CURRENT_SOURCE_DIR}/deprecated
  ${CMAKE_CURRENT_SOURCE_DIR}/dictBuilder
  ${CMAKE_CURRENT_SOURCE_DIR}/dll
  ${CMAKE_CURRENT_SOURCE_DIR}/legacy)

set(zstd_sources
./decompress/zstd_decompress.c
./decompress/huf_decompress.c
./compress/zstd_lazy.c
./compress/zstdmt_compress.c
./compress/zstd_double_fast.c
./compress/zstd_fast.c
./compress/hist.c
./compress/fse_compress.c
./compress/zstd_opt.c
./compress/zstd_compress.c
./compress/huf_compress.c
./compress/zstd_ldm.c
./common/xxhash.c
./common/fse_decompress.c
./common/pool.c
./common/zstd_common.c
./common/error_private.c
./common/debug.c
./common/threading.c
./common/entropy_common.c
./deprecated/zbuff_compress.c
./deprecated/zbuff_decompress.c
./deprecated/zbuff_common.c
./dictBuilder/zdict.c
./dictBuilder/divsufsort.c
./dictBuilder/cover.c
./legacy/zstd_v03.c
./legacy/zstd_v07.c
./legacy/zstd_v06.c
./legacy/zstd_v01.c
./legacy/zstd_v02.c
./legacy/zstd_v04.c
./legacy/zstd_v05.c
	)

add_library (ZSTD SHARED ${zstd_sources})
target_link_libraries (ZSTD ${zstd_dependencies})
target_include_directories(ZSTD 
  PUBLIC common compress decompress deprecated legacy dictBuilder
  )

install (TARGETS ZSTD LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
set (zstd_headers
     common/bitstream.h  common/error_private.h  common/pool.h         common/zstd_internal.h             compress/zstd_fast.h        compress/zstd_opt.h       legacy/zstd_legacy.h  legacy/zstd_v04.h
common/compiler.h   common/fse.h            common/threading.h    compress/hist.h                    compress/zstd_lazy.h        deprecated/zbuff.h        legacy/zstd_v01.h     legacy/zstd_v05.h
common/cpu.h        common/huf.h            common/xxhash.h       compress/zstd_compress_internal.h  compress/zstd_ldm.h         dictBuilder/divsufsort.h  legacy/zstd_v02.h     legacy/zstd_v06.h
common/debug.h      common/mem.h            common/zstd_errors.h  compress/zstd_double_fast.h        compress/zstdmt_compress.h  dictBuilder/zdict.h       legacy/zstd_v03.h     legacy/zstd_v07.h   
zstd.h )

install (FILES ${zstd_headers} DESTINATION include)
