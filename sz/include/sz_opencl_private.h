/**
 * this file is for any non-public methods or structures used by the sz_opencl implementation
 * that do not also need to be included on a device kernel
 */
#include <CL/cl.hpp>

struct sz_opencl_state
{
  struct error_t
  {
    cl_int code = 0;
    const char* str = nullptr;
  } error;
  cl::Platform platform;
  cl::Device device;
  cl::Context context;
  cl::CommandQueue queue;
  cl::Kernel calculate_regression_coefficents;
};
