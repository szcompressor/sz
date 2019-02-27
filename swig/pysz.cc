#include "pysz.h"
#include <sstream>
#include <stdexcept>
#include <zlib.h>

ConfigBuilder::ConfigBuilder() {
  building.params.absErrBound = 1e-4;
  building.params.dataType = SZ_FLOAT;
  building.params.errorBoundMode = ABS;
  building.params.gzipMode = Z_BEST_SPEED;
  building.params.losslessCompressor = ZSTD_COMPRESSOR;
  building.params.maxRangeRadius = 32768;
  building.params.max_quant_intervals = 65536;
  building.params.predThreshold = 0.99;
  building.params.predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
  building.params.psnr = 80;
  building.params.pw_relBoundRatio = 1e-2;
  building.params.pwr_type = SZ_PWR_MIN_TYPE;
  building.params.quantization_intervals = 0;
  building.params.randomAccess = 0;
  building.params.relBoundRatio = 1e-4;
  building.params.sampleDistance = 100;
  building.params.segment_size = 25;
  building.params.snapshotCmprStep = 5;
  building.params.sol_ID = SZ;
  building.params.szMode = SZ_BEST_SPEED;
}


ConfigBuilder& ConfigBuilder::absErrBound(double value) noexcept { building.params.absErrBound = value; return *this; }
ConfigBuilder& ConfigBuilder::dataType(int value) noexcept { building.params.dataType = value; return *this; }
ConfigBuilder& ConfigBuilder::errorBoundMode(int  value) noexcept { building.params.errorBoundMode = value; return *this; }
ConfigBuilder& ConfigBuilder::gzipMode(int value) noexcept { building.params.gzipMode = value; return *this; }
ConfigBuilder& ConfigBuilder::losslessCompressor(int value) noexcept { building.params.losslessCompressor = value; return *this; }
ConfigBuilder& ConfigBuilder::maxRangeRadius(unsigned int value) noexcept { building.params.maxRangeRadius = value; return *this; }
ConfigBuilder& ConfigBuilder::max_quant_intervals(unsigned int value) noexcept { building.params.max_quant_intervals = value; return *this; }
ConfigBuilder& ConfigBuilder::predThreshold(float value) noexcept { building.params.predThreshold = value; return *this; }
ConfigBuilder& ConfigBuilder::predictionMode(int value) noexcept { building.params.predictionMode = value; return *this; }
ConfigBuilder& ConfigBuilder::psnr(double value) noexcept { building.params.psnr = value; return *this; }
ConfigBuilder& ConfigBuilder::pw_relBoundRatio(double value) noexcept { building.params.pw_relBoundRatio = value; return *this; }
ConfigBuilder& ConfigBuilder::pwr_type(int value) noexcept { building.params.pwr_type = value; return *this; }
ConfigBuilder& ConfigBuilder::quantization_intervals(unsigned int value) noexcept { building.params.quantization_intervals = value; return *this; }
ConfigBuilder& ConfigBuilder::randomAccess(int value) noexcept { building.params.randomAccess = value; return *this; }
ConfigBuilder& ConfigBuilder::relBoundRatio(double value) noexcept { building.params.relBoundRatio = value; return *this; }
ConfigBuilder& ConfigBuilder::sampleDistance(int value) noexcept { building.params.sampleDistance = value; return *this; }
ConfigBuilder& ConfigBuilder::segment_size(int value) noexcept { building.params.segment_size = value; return *this; }
ConfigBuilder& ConfigBuilder::snapshotCmprStep(int value) noexcept { building.params.snapshotCmprStep = value; return *this; }
ConfigBuilder& ConfigBuilder::sol_ID(int value) noexcept { building.params.sol_ID = value; return *this; }
ConfigBuilder& ConfigBuilder::szMode(int value) noexcept { building.params.szMode = value; return *this; }
Config ConfigBuilder::build() { return building;} 

Compressor::Compressor(Config config) {
  int ret = SZ_Init_Params(&config.params);
  if(ret != SZ_SCES) {
    std::ostringstream ss;
    ss << "SZ Init Error: " << ret;
    throw std::runtime_error(ss.str());
  }
}


Compressor::Compressor(std::string const& configfile) {
  int ret = SZ_Init(configfile.c_str());
  if(ret != SZ_SCES) {
    std::ostringstream ss;
    ss << "SZ Init Error: " << ret;
    throw std::runtime_error(ss.str());
  }
}


Compressor::~Compressor() {
  SZ_Finalize();
}
