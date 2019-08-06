
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/include/sz_opencl_config.h.in
	${CMAKE_CURRENT_BINARY_DIR}/include/sz_opencl_config.h @ONLY)

add_library (SZ SHARED 
  src/ArithmeticCoding.c
  src/ByteToolkit.c
  src/CompressElement.c
  src/DynamicByteArray.c
  src/DynamicDoubleArray.c
  src/DynamicFloatArray.c
  src/DynamicIntArray.c
  src/Huffman.c
  src/TightDataPointStorageD.c
  src/TightDataPointStorageF.c
  src/TightDataPointStorageI.c
  src/TypeManager.c
  src/VarSet.c
  src/callZlib.c
  src/conf.c
  src/dataCompression.c
  src/dictionary.c
  src/iniparser.c
  src/pastri.c
  src/rw.c
  src/rw_interface.F90
  src/rwf.c
  src/sz.c
  src/sz_cuda.cu
  src/sz_double.c
  src/sz_double_pwr.c
  src/sz_double_ts.c
  src/sz_float.c
  src/sz_float_pwr.c
  src/sz_float_ts.c
  src/sz_int16.c
  src/sz_int32.c
  src/sz_int64.c
  src/sz_int8.c
  src/sz_interface.F90
  src/sz_omp.c
  src/sz_opencl.cc
  src/sz_opencl.cl
  src/sz_opencl_host_utils.cc
  src/sz_uint16.c
  src/sz_uint32.c
  src/sz_uint64.c
  src/sz_uint8.c
  src/szd_double.c
  src/szd_double_pwr.c
  src/szd_double_ts.c
  src/szd_float.c
  src/szd_float_pwr.c
  src/szd_float_ts.c
  src/szd_int16.c
  src/szd_int32.c
  src/szd_int64.c
  src/szd_int8.c
  src/szd_uint16.c
  src/szd_uint32.c
  src/szd_uint64.c
  src/szd_uint8.c
  src/szf.c
  src/utility.c
  src/CacheTable.c
  src/MultiLevelCacheTable.c
  src/MultiLevelCacheTableWideInterval.c
  )
target_include_directories(SZ 
  PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries (SZ PUBLIC ZLIB ZSTD m ${OpenCL_LIBRARIES})
target_compile_options(SZ
	PUBLIC -DHAVE_RANDOMACCESS
	PRIVATE $<$<CONFIG:Debug>:-Wall -Wextra -Wpedantic -Wno-unused-parameter>
	)

install (TARGETS SZ EXPORT SZConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install (EXPORT SZConfig DESTINATION share/SZ/cmake)
install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sz)
install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sz)
export(TARGETS SZ FILE SZ.cmake)
