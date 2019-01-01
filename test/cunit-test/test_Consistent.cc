#include <vector>
#include <random>
#include <algorithm>
#include <functional>

#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"
#include "RegressionTest.hpp"

extern "C" {
int init_suite()
{
	return 0;
}

int clean_suite()
{
	return 0;
}

void test_all_functions() {

	std::vector<Compressors> compressor_functions{
			{
					SZ_compress_float_3D_MDQ_decompression_random_access_with_blocked_regression,
					SZ_decompress_args_randomaccess_float
			}
	};


	auto num_random_test_cases = 4;

	for(auto functions: compressor_functions) {
		auto compressor = functions.compressor;
		auto decompressor = functions.decompressor;

		test_identical_output_compression_random(num_random_test_cases, compressor, compressor, decompressor,
												 decompressor);
		test_identical_output_compression_deterministic(compressor, compressor, decompressor, decompressor);
	}
}
}

int main(int argc, char *argv[])
{
	unsigned int num_failures = 0;
	CU_ErrorCode error_code = CUE_SUCCESS;
	if (CUE_SUCCESS != CU_initialize_registry())
	{
		return CU_get_error();
	}

	CU_pSuite suite = CU_add_suite("test_opencl_suite", init_suite, clean_suite);
	if(suite == nullptr) {
		goto error;
	}

	if(CU_add_test(suite, "test_all_functions", test_all_functions) == nullptr) {
		goto error;
	}

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	error_code = (CU_ErrorCode)(error_code | CU_get_error());
	CU_basic_show_failures(CU_get_failure_list());
	num_failures = CU_get_number_of_failures();

error:
	CU_cleanup_registry();
	error_code = (CU_ErrorCode)(error_code | CU_get_error());
	return  num_failures ||  error_code;
}
