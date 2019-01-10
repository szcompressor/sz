/**
 * Contains methods that load or launch kernels
 */
#include <cstddef>
#include <stdexcept>
#include <string>
#include <CL/cl.hpp>

//forward declare to avoid dependency on sz_opencl.h or sz_opencl_private.h
class sz_opencl_state;
class sz_opencl_sizes;

enum class copy_mode {
  TO,
  FROM,
  ZERO_FROM,
  ZERO,
  TO_FROM
};

struct buffer_copy_info
{
  void* ptr;
  std::size_t size;
  cl_mem_flags flags;
  copy_mode copy;
};

class sz_opencl_exception : public std::runtime_error
{
  using std::runtime_error::runtime_error;
};

std::string get_sz_kernel_sources();


std::vector<cl::Event> run_kernel(cl::Kernel kernel,
                                  cl::NDRange const &global,
                                  const sz_opencl_state *state,
                                  const sz_opencl_sizes *sizes,
                                  const std::vector<buffer_copy_info> &buffer_info);
