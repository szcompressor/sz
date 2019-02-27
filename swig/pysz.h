#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "sz.h"
#include "pysz_private.h"


struct Config {
  sz_params params;
};


class ConfigBuilder {
  public:
  ConfigBuilder();
  Config build();

	ConfigBuilder& absErrBound(double value) noexcept;
	ConfigBuilder& dataType(int value) noexcept;
	ConfigBuilder& errorBoundMode(int  value) noexcept;
	ConfigBuilder& gzipMode(int value) noexcept;
	ConfigBuilder& losslessCompressor(int value) noexcept;
	ConfigBuilder& maxRangeRadius(unsigned int value) noexcept;
	ConfigBuilder& max_quant_intervals(unsigned int value) noexcept;
	ConfigBuilder& predThreshold(float value) noexcept;
	ConfigBuilder& predictionMode(int value) noexcept;
	ConfigBuilder& psnr(double value) noexcept;
	ConfigBuilder& pw_relBoundRatio(double value) noexcept;
	ConfigBuilder& pwr_type(int value) noexcept;
	ConfigBuilder& quantization_intervals(unsigned int value) noexcept;
	ConfigBuilder& randomAccess(int value) noexcept;
	ConfigBuilder& relBoundRatio(double value) noexcept;
	ConfigBuilder& sampleDistance(int value) noexcept;
	ConfigBuilder& segment_size(int value) noexcept;
	ConfigBuilder& snapshotCmprStep(int value) noexcept;
	ConfigBuilder& sol_ID(int value) noexcept;
	ConfigBuilder& szMode(int value) noexcept;
  private:
  Config building;
};

class Compressor {
  public:
  /*
   * Curently SZ has global state that needs to initialized only once; Thus
   * compressor should be a move only type for now and use RAII to manage the
   * global state
   */
  Compressor(Config config);
  Compressor(std::string const& config);
  ~Compressor();
  
  Compressor(Compressor const&)=delete;
  Compressor& operator=(Compressor&)=delete;

  //
  //Yes, this could be done using varatic arguments, but due to limitations in how SWIG handles templates
  //We need to explicitly list all of these out.
  //
  //even though sz supports 5d arrays, they aren't supported well by numpy, omit it for now
  //

  template<class T>
  std::string Compress1(T* data, size_t r1) {
    size_t outsize = 0;
    char* tmp = reinterpret_cast<char*>(SZ_compress(SZTypeToTypeID<T>::value, data, &outsize, 0, 0, 0, 0, r1));
    return std::string(tmp, outsize);
  }
  template<class T>
  std::string Compress2(T* data, size_t r1, size_t r2) {
    size_t outsize = 0;
    auto tmp = reinterpret_cast<char*>(SZ_compress(SZTypeToTypeID<T>::value, data, &outsize, 0, 0, 0, r2, r1));
    return std::string(tmp, outsize);
  }
  template<class T>
  std::string Compress3(T* data, size_t r1, size_t r2, size_t r3) {
    size_t outsize = 0;
    auto tmp = reinterpret_cast<char*>(SZ_compress(SZTypeToTypeID<T>::value, data, &outsize, 0, 0, r3, r2, r1));
    return std::string(tmp, outsize);
  }
  template<class T>
  std::string Compress4(T* data, size_t r1, size_t r2, size_t r3, size_t r4) {
    size_t outsize = 0;
    auto tmp = reinterpret_cast<char*>(SZ_compress(SZTypeToTypeID<T>::value, data, &outsize, 0, r4, r3, r2, r1));
    return std::string(tmp, outsize);

  }

  template<class T>
  std::vector<T> Decompress(std::string data, std::vector<int> r) {
    T* decompressed = nullptr;
    size_t len = 0;
    switch(r.size()) {
      case 1:
        decompressed = (T*)SZ_decompress(SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,0,0,0,r[0]);
        len = r[0];
        break;
      case 2:
        decompressed = (T*)SZ_decompress(SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,0,0,r[1],r[0]);
        len = r[0]*r[1];
        break;
      case 3:
        decompressed = (T*)SZ_decompress(SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,0,r[2],r[1],r[0]);
        len = r[0]*r[1]*r[2];
        break;
      case 4:
        decompressed = (T*)SZ_decompress(SZTypeToTypeID<T>::value, (unsigned char*)(data.c_str()), data.length(), 0,r[3],r[2],r[1],r[0]);
        len = r[0]*r[1]*r[2]*r[3];
        break;
      default:
        printf("%zu dimensional arrays not supported\n", r.size());
    }
    return std::vector<T>(decompressed, decompressed+len);
  }
};
