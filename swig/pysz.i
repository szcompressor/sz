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

%module pysz
%feature("autodoc", 1);

%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_STRICT_BYTE_CHAR
#include "pysz.h"
#include "defines.h"
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
%include "defines.h"

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
