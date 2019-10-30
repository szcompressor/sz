#ifndef EXAFELSZ_H
#define EXAFELSZ_H

#include <sz.h>

#include <stdint.h>
#include <stdlib.h>

typedef struct exafelSZ_params{
  uint8_t *peaks;
  uint8_t *calibPanel;

  uint8_t binSize; //Binning: (pr->binSize x pr->binSize) to (1 x 1)
  double tolerance; //SZ pr->tolerance
  uint8_t szDim; //1D/2D/3D compression/decompression
  //uint8_t szBlockSize; //Currently unused
  uint8_t peakSize; //MUST BE ODD AND NOT EVEN! Each peak will have size of: (peakSize x peakSize)
  
  // string origFilename;
  // string compFilename;
  // string roiFilename;
  // string roniFilename;
   
  uint64_t nEvents;
  uint64_t panels;
  uint64_t rows;
  uint64_t cols;
  
  //CALCULATED VARIBALES:
  uint64_t binnedRows;
  uint64_t binnedCols;
  uint8_t peakRadius; //Will be calculated using peakSize

} exafelSZ_params;

unsigned char * exafelSZ_Compress(void* _origData,
                       void* _pr,
                       size_t *compressedSize);
void* exafelSZ_Decompress(unsigned char*_compressedBuffer,
                         void *_pr,
                         size_t compressedSize);


#endif