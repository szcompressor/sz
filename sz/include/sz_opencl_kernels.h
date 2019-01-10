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
    , pred_buffer_size((block_size + 1) * (block_size + 1) * (block_size + 1))
    , parallel_pred_buffer_size(pred_buffer_size * num_blocks)
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
};
