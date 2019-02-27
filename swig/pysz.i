%module pysz

%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_STRICT_BYTE_CHAR
#include "pysz.h"
%}

%init %{
import_array();
%}

%include <std_vector.i>
%include <std_string.i>
%include "numpy.i"

namespace std {
  %template(vectori) vector<int>;
  %template(vectorf) vector<float>;
  %template(vectord) vector<double>;
};

%apply (float* INPLACE_ARRAY1, int DIM1 ) {(float* data, size_t r1)}
%apply (float* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(float* data, size_t r1, size_t r2)}
%apply (float* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(float* data, size_t r1, size_t r2, size_t r3)}
%apply (float* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(float* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (double* INPLACE_ARRAY1, int DIM1 ) {(double* data, size_t r1)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(double* data, size_t r1, size_t r2)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(double* data, size_t r1, size_t r2, size_t r3)}
%apply (double* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(double* data, size_t r1, size_t r2, size_t r3, size_t r4)}


%include "pysz.h"

%extend Compressor {
  %template(CompressFloat1) Compress1<float>;
  %template(CompressDouble1) Compress1<double>;
  %template(CompressFloat2) Compress2<float>;
  %template(CompressDouble2) Compress2<double>;
  %template(CompressFloat3) Compress3<float>;
  %template(CompressDouble3) Compress3<double>;
  %template(CompressFloat4) Compress4<float>;
  %template(CompressDouble4) Compress4<double>;


  %template(DecompressFloat) Decompress<float>;
  %template(DecompressDouble) Decompress<double>;

  %pythoncode %{
    import numpy

    __Compress = {
      (1, numpy.dtype('float64')): CompressDouble1,
      (2, numpy.dtype('float64')): CompressDouble2,
      (3, numpy.dtype('float64')): CompressDouble3,
      (4, numpy.dtype('float64')): CompressDouble4,
      (1, numpy.dtype('float32')): CompressFloat1,
      (2, numpy.dtype('float32')): CompressFloat2,
      (3, numpy.dtype('float32')): CompressFloat3,
      (4, numpy.dtype('float32')): CompressFloat4
    }

    __Decompress = {
      numpy.dtype('float32'): DecompressFloat,
      numpy.float32: DecompressFloat,
      numpy.dtype('float64'): DecompressDouble,
      numpy.float64: DecompressDouble
    }

    def Compress(self, array):
      length = len(array.shape)
      dtype = array.dtype
      return self.__Compress[length, dtype](self, array)

    def Decompress(self, bytes, dims, dtype):
      try:
        values = self.__Decompress[dtype](self, bytes, dims)
        return self.numpy.reshape(values, dims)
      except KeyError as e:
        raise TypeError("type {} not supported".format(i.args[0]))


      

  %}

}
