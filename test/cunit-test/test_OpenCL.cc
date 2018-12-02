#include <vector>
#include <random>
#include <algorithm>

#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"
#include "zlib.h"

namespace {
	struct sz_opencl_state* state = nullptr;
	const int size = 64*64*64;
	const int num_random_test_cases = 32;

	template <class InputIt, class InputIt2>
	unsigned int count_non_equal(InputIt it1, InputIt2 it2, size_t size)
	{
		unsigned int not_equal_elms = 0; 
		for (size_t i = 0; i < size; ++i) {
			if(it1[i] != it2[i]) not_equal_elms++;
		}
		return not_equal_elms;
	}
}

extern "C" {
int init_suite()
{
	int rc = sz_opencl_init(&state);
	rc |= sz_opencl_error_code(state);
	return rc;
}

int clean_suite()
{
	int rc = sz_opencl_release(&state);
	return rc;
}

inline sz_params
sz_default_config()
{
  sz_params params;
  params.dataType = SZ_FLOAT;
  params.max_quant_intervals = 65536;
  params.quantization_intervals = 0;
  params.predThreshold = 0.99;
  params.szMode = SZ_BEST_SPEED;
  params.gzipMode = Z_BEST_SPEED;
  params.errorBoundMode = ABS;
  params.absErrBound = 1e-4;
  params.relBoundRatio = 1e-4;
  params.psnr = 80;
  params.pw_relBoundRatio = 1e-2;
  params.segment_size = 25;
  params.pwr_type = SZ_PWR_MIN_TYPE;
  params.sampleDistance = 100;
  params.sol_ID = SZ;
  params.predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	params.randomAccess = true;
  return params;
}


void test_valid_opencl()
{
	int rc = sz_opencl_check(state);
	CU_ASSERT_EQUAL(rc, 0);
}

void test_identical_output()
{
	std::seed_seq seed{0};
	std::mt19937 gen{seed};
	std::uniform_real_distribution<float> dist;
	auto rand = [&gen,&dist](){return dist(gen);};

	sz_params conf_params{sz_default_config()};
	SZ_Init_Params(&conf_params);

	for (int test_case = 0; test_case < num_random_test_cases; ++test_case) {
		std::vector<float> test_array(size);
		std::generate(std::begin(test_array), std::end(test_array), rand);
		std::vector<float> test_arraycl(test_array);

		size_t size_existing = 0;
		sz_params conf_params{sz_default_config()};
		SZ_Init_Params(&conf_params);
		unsigned char* data = SZ_compress_float_3D_MDQ_decompression_random_access_with_blocked_regression(test_array.data(), 64, 64, 64, .99, &size_existing);
		SZ_Finalize();

		size_t size_cl = 0;
		sz_params conf_params_cl{sz_default_config()};
		SZ_Init_Params(&conf_params_cl);
		unsigned char* data_cl = sz_compress_float3d_opencl(test_arraycl.data(), 64, 64, 64, .99, &size_cl);
		SZ_Finalize();

		auto num_not_equal = count_non_equal(data, data_cl, size_existing);
		CU_ASSERT_EQUAL(num_not_equal, 0);

		free(data);
		free(data_cl);
	}

}

int main(int argc, char *argv[])
{
	unsigned int num_failures = 0;
	if (CUE_SUCCESS != CU_initialize_registry())
	{
		return CU_get_error();
	}

	CU_pSuite suite = CU_add_suite("test_opencl_suite", init_suite, clean_suite);
	if(suite == nullptr) {
		goto error;
	}

	if(CU_add_test(suite, "test_valid_opencl", test_valid_opencl) == nullptr ||
		 CU_add_test(suite, "test_identical_output", test_identical_output) == nullptr) {
		goto error;
	}

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_basic_show_failures(CU_get_failure_list());
	num_failures = CU_get_number_of_failures();

error:
	CU_cleanup_registry();
	return  num_failures ||  CU_get_error();
}

}
