#ifdef __cplusplus
#include <CL/cl_platform.h>
#else
#define cl_ulong ulong
#define cl_float float
#endif

struct sz_opencl_sizes
{

#ifdef __cplusplus
  sz_opencl_sizes(int block_size, int r1, int r2, int r3)
    : r1(r1)
    , r2(r2)
    , r3(r3)
    , block_size(block_size)
    , num_x((r1 - 1) / block_size + 1)
    , num_y((r2 - 1) / block_size + 1)
    , num_z((r3 - 1) / block_size + 1)
    , max_num_block_elements(block_size * block_size * block_size)
    , num_blocks(num_x * num_y * num_z)
    , num_elements(r1 * r2 * r3)
    , dim0_offset(r2 * r3)
    , dim1_offset(r3)
    , params_offset_b(num_blocks)
    , params_offset_c(2 * num_blocks)
    , params_offset_d(3 * num_blocks)
    , pred_buffer_block_size(block_size + 1)
    , strip_dim0_offset(pred_buffer_block_size * pred_buffer_block_size)
    , strip_dim1_offset(pred_buffer_block_size)
    , unpred_data_max_size(max_num_block_elements)
    , reg_params_buffer_size(num_blocks * 4)
    , pred_buffer_size(pred_buffer_block_size * pred_buffer_block_size * pred_buffer_block_size)
    , parallel_pred_buffer_size(pred_buffer_size * num_blocks)
    , block_dim0_offset(pred_buffer_block_size*pred_buffer_block_size)
    , block_dim1_offset(pred_buffer_block_size)
  {}
#endif

  const cl_ulong r1, r2, r3;
  const cl_ulong block_size;
  const cl_ulong num_x;
  const cl_ulong num_y;
  const cl_ulong num_z;
  const cl_ulong max_num_block_elements;
  const cl_ulong num_blocks;
  const cl_ulong num_elements;
  const cl_ulong dim0_offset;
  const cl_ulong dim1_offset;
  const cl_ulong params_offset_b;
  const cl_ulong params_offset_c;
  const cl_ulong params_offset_d;
  const cl_ulong pred_buffer_block_size;
  const cl_ulong strip_dim0_offset;
  const cl_ulong strip_dim1_offset;
  const cl_ulong unpred_data_max_size;
  const cl_ulong reg_params_buffer_size;
  const cl_ulong pred_buffer_size;
  const cl_ulong parallel_pred_buffer_size;
  const cl_ulong block_dim0_offset;
  const cl_ulong block_dim1_offset;
};

struct sz_opencl_decompress_positions {
#ifdef __cplusplus
  sz_opencl_decompress_positions(int block_size, int s1, int s2, int s3, int e1, int e2, int e3)
      : start_elm1(s1),
  start_elm2(s2),
  start_elm3(s3),
  end_elm1(e1),
  end_elm2(e2),
  end_elm3(e3),
  start_block1(s1 / block_size),
  start_block2(s2 / block_size),
  start_block3(s3 / block_size),
  end_block1((e1 - 1) / block_size + 1),
  end_block2((e2 - 1) / block_size + 1),
  end_block3((e3 - 1) / block_size + 1),
  num_data_blocks1(end_block1-start_block1),
  num_data_blocks2(end_block2-start_block2),
  num_data_blocks3(end_block3-start_block3),
  data_elms1(e1-s1),
  data_elms2(e2-s2),
  data_elms3(e3-s3),
  data_buffer_size(data_elms1 * data_elms2 * data_elms3),
  dec_block_dim1_offset(num_data_blocks3 * block_size),
  dec_block_dim0_offset(dec_block_dim1_offset * num_data_blocks2 * block_size)
  {}

#endif
  const cl_ulong start_elm1, start_elm2, start_elm3;
  const cl_ulong end_elm1, end_elm2, end_elm3;
  const cl_ulong start_block1, start_block2, start_block3;
  const cl_ulong end_block1, end_block2, end_block3;
  const cl_ulong num_data_blocks1, num_data_blocks2, num_data_blocks3;
  const cl_ulong data_elms1, data_elms2, data_elms3;
  const cl_ulong data_buffer_size;
  const cl_ulong dec_block_dim1_offset;
  const cl_ulong dec_block_dim0_offset;
};
