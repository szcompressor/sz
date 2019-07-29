aux_source_directory(src sz_sources)

add_library (SZ SHARED ${sz_sources})
target_link_libraries (SZ ${ZLIB_LIBRARIES} ${ZSTD_LIBRARIES})
if (UNIX)
  target_link_libraries (SZ m)
endif ()

target_compile_options(SZ
	PRIVATE $<$<CONFIG:Debug>:-O3 -Wall -Wextra -Wpedantic -Wno-unused-parameter>
	)

install (TARGETS SZ LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin)

set (sz_headers
include/ArithmeticCoding.h    include/DynamicFloatArray.h                 include/rw.h              include/szd_int8.h       include/sz_float_pwr.h  include/sz_uint32.h
include/ByteToolkit.h         include/DynamicIntArray.h                   include/szd_double.h      include/sz_double.h      include/sz_float_ts.h   include/sz_uint64.h
include/CacheTable.h          include/Huffman.h                           include/szd_double_pwr.h  include/sz_double_pwr.h  include/sz.h            include/sz_uint8.h
include/callZlib.h            include/szd_double_ts.h                     include/sz_double_ts.h    include/sz_int16.h       include/TightDataPointStorageD.h
include/CompressElement.h     include/MultiLevelCacheTable.h              include/szd_float.h       include/szd_uint16.h     include/sz_int32.h      include/TightDataPointStorageF.h
include/conf.h                include/MultiLevelCacheTableWideInterval.h  include/szd_float_pwr.h   include/szd_uint32.h     include/sz_int64.h      include/TightDataPointStorageI.h
include/dataCompression.h     include/pastriD.h                           include/szd_float_ts.h    include/szd_uint64.h     include/sz_int8.h       include/TypeManager.h
include/pastriF.h             include/szd_int16.h                         include/szd_uint8.h       include/sz_omp.h         include/utility.h
include/DynamicByteArray.h    include/pastriGeneral.h                     include/szd_int32.h       include/szf.h            include/sz_opencl.h     include/VarSet.h
include/DynamicDoubleArray.h  include/pastri.h                            include/szd_int64.h       include/sz_float.h       include/sz_uint16.h )

install (FILES ${sz_headers} DESTINATION include)
