/**
 *
 * A set of python bindings for SZ 
 * 
 * Developed by Robert Underwood while he was at Clemson University
 * This material is based upon work supported by the National Science 
 * Foundation under Grant No. 1633608.
 * 
 * Copyright © 2019 Robert Underwood
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY Robert Underwood ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Robert Underwood BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation
 * are those of the authors and should not be interpreted as representing
 * official policies, either expressed or implied, of Robert Underwood.
 */
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
