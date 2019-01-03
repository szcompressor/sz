#include "sz_opencl_kernels.h"

kernel void
calculate_regression_coefficents(
		__global const float* oriData,
		__global struct sz_opencl_sizes const* sizes,
        __global float* reg_params_pos,
		__global float* pred_buffer,
        __global const float* data_pos,
		__global float* pred_buffer_pos)
{
  ulong i = get_global_id(0);
  ulong j = get_global_id(1);
  ulong k = get_global_id(2);
  data_pos = oriData + i * sizes->block_size * sizes->dim0_offset +
             j * sizes->block_size * sizes->dim1_offset + k * sizes->block_size;
  pred_buffer_pos = pred_buffer;
  for (size_t ii = 0; ii < sizes->block_size; ii++) {
    for (size_t jj = 0; jj < sizes->block_size; jj++) {
      for (size_t kk = 0; kk < sizes->block_size; kk++) {
        int ii_ = (i * sizes->block_size + ii < sizes->r1)
                    ? ii
                    : sizes->r1 - 1 - i * sizes->block_size;
        int jj_ = (j * sizes->block_size + jj < sizes->r2)
                    ? jj
                    : sizes->r2 - 1 - j * sizes->block_size;
        int kk_ = (k * sizes->block_size + kk < sizes->r3)
                    ? kk
                    : sizes->r3 - 1 - k * sizes->block_size;
        *pred_buffer_pos = *(data_pos + ii_ * sizes->dim0_offset +
                             jj_ * sizes->dim1_offset + kk_);
        pred_buffer_pos++;
      }
    }
  }

  /*Calculate regression coefficients*/
  {
    __global float* cur_data_pos = pred_buffer;
    float fx = 0.0;
    float fy = 0.0;
    float fz = 0.0;
    float f = 0;
    float sum_x, sum_y;
    float curData;
    for (size_t i = 0; i < sizes->block_size; i++) {
      sum_x = 0;
      for (size_t j = 0; j < sizes->block_size; j++) {
        sum_y = 0;
        for (size_t k = 0; k < sizes->block_size; k++) {
          curData = *cur_data_pos;
          sum_y += curData;
          fz += curData * k;
          cur_data_pos++;
        }
        fy += sum_y * j;
        sum_x += sum_y;
      }
      fx += sum_x * i;
      f += sum_x;
    }
    float coeff =
      1.0 / (sizes->block_size * sizes->block_size * sizes->block_size);
    reg_params_pos[0] = (2 * fx / (sizes->block_size - 1) - f) * 6 * coeff /
                        (sizes->block_size + 1);
    reg_params_pos[sizes->params_offset_b] =
      (2 * fy / (sizes->block_size - 1) - f) * 6 * coeff /
      (sizes->block_size + 1);
    reg_params_pos[sizes->params_offset_c] =
      (2 * fz / (sizes->block_size - 1) - f) * 6 * coeff /
      (sizes->block_size + 1);
    reg_params_pos[sizes->params_offset_d] =
      f * coeff -
      ((sizes->block_size - 1) * reg_params_pos[0] / 2 +
       (sizes->block_size - 1) * reg_params_pos[sizes->params_offset_b] / 2 +
       (sizes->block_size - 1) * reg_params_pos[sizes->params_offset_c] / 2);
  }
  reg_params_pos++;
}
