#include "exafelSZ.h"

void exafelSZ_params_process(exafelSZ_params*pr){
  pr->binnedRows=(pr->rows+pr->binSize-1)/pr->binSize;
  pr->binnedCols=(pr->cols+pr->binSize-1)/pr->binSize;
  
  pr->peakRadius=(pr->peakSize-1)/2;
}

void exafelSZ_params_checkDecomp(exafelSZ_params*pr){
  if(pr->calibPanel==NULL){
    // cout<<"ERROR: calibPanel is NULL : calibPanel="<<calibPanel<<endl;
    assert(0);
  }
  if(pr->binSize<1 || pr->tolerance<0 || pr->szDim<1 || pr->szDim>3){
    // cout<<"ERROR: Something wrong with the following:"<<endl;
    // cout<<"binSize="<<binSize<<endl;
    // cout<<"tolerance="<<tolerance<<endl;
    // cout<<"szDim="<<szDim<<endl;
  }
  if(!(pr->peakSize%2)){
    // cout<<"ERROR: peakSize = "<<peakSize<<" cannot be even. It must be odd!"<<endl;
    assert(0);
  }  
  if(pr->nEvents<1 || pr->panels<0 || pr->rows<1 || pr->cols<1){
    // cout<<"ERROR: Something wrong with the following:"<<endl;
    // cout<<"nEvents="<<nEvents<<endl;
    // cout<<"panels="<<panels<<endl;
    // cout<<"rows="<<rows<<endl;
    // cout<<"cols="<<cols<<endl;
    assert(0);
  }
}

void exafelSZ_params_checkComp(exafelSZ_params*pr){
  if(pr->peaks==NULL){
    // cout<<"ERROR: peaks is NULL : peaks="<<peaks<<endl;
    assert(0);
  }
  exafelSZ_params_checkDecomp(pr);
}

void exafelSZ_params_print(){
  // outs<<"Configuration (exafelSZ_params) : "<<endl;
  // outs<<"SMOOTHING: NO"<<"  (ROI and RONI are NOT replaced by local avg values)"<<endl;
  // outs<<"binSize:"<<binSize<<endl;
  // outs<<"tolerance:"<<tolerance<<endl;
  // outs<<"szDim:"<<szDim<<endl;
  // outs<<"peakSize:"<<peakSize<<endl;
  // outs<<"nEvents:"<<nEvents<<" (# of events per batch)"<<endl;
  // outs<<"panels:"<<panels<<" (Panels per event)"<<endl;
  // outs<<"rows:"<<rows<<" (Rows per panel)"<<endl;
  // outs<<"cols:"<<cols<<" (Columns per panel)"<<endl;
  // outs<<endl;
  // outs<<"CALCULATED VARIABLES"<<endl;
  // outs<<"binnedRows:"<<binnedRows<<" (Rows per panel after binning)"<<endl;
  // outs<<"binnedCols:"<<binnedCols<<" (Columns per panel after binning)"<<endl;
  // outs<<"peakRadius:"<<peakRadius<<" (Peak radius = (peakSize-1)/2 )"<<endl;
  // outs<<endl;
}

//*********************************************************************************
//*********************************************************************************
//*********************************************************************************

//Index Calculator
inline int calcIdx_4D(int i3, int i2, int i1, int i0, int size2, int size1, int size0){ 
  return i0+size0*(i1+size1*(i2+size2*i3));
}
inline int calcIdx_3D(int i2, int i1, int i0, int size1, int size0){ 
  return i0+size0*(i1+size1*i2);
}
inline int calcIdx_2D(int i1, int i0, int size0){ 
  return i0+size0*i1;
}

unsigned char * exafelSZ_Compress(void* _origData,
                       void* _pr,
                       size_t *compressedSize)
{
  float *origData=(float*)_origData;
  exafelSZ_params *pr=(exafelSZ_params*)_pr;
  
  exafelSZ_params_process(pr); 
  exafelSZ_params_checkDecomp(pr); 
  
  uint8_t *roiM=(uint8_t*)malloc(pr->nEvents*pr->panels*pr->rows*pr->cols) ;
  float *roiData=(float*)malloc(pr->nEvents*pr->panels*pr->rows*pr->cols*sizeof(float)) ;
  float *binnedData=(float*)malloc(pr->nEvents*pr->panels*pr->binnedRows*pr->binnedCols*sizeof(float)) ;
  //float *binnedData=(float*)malloc(pr->nEvents*pr->panels*pr->rows*pr->cols*sizeof(float)) ;
  
  int e,p,r,c,pk,ri,ci,br,bc,roii,bi;
  
  //Generate the ROI mask: NOTE: 0 means affirmative in ROI mask! This comes from the python scripts!
  //First, initialize with calibration panel:
  for(e=0;e<pr->nEvents;e++){ //Event
    for(p=0;p<pr->panels;p++){ //Panel
      for(r=0;r<pr->rows;r++){ //Row
        for(c=0;c<pr->cols;c++){ //Column
          roiM[calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)]=pr->calibPanel[calcIdx_2D(r,c,pr->cols)];
        }
      }
    }
  }
  uint64_t peaksBytePos=0; //Position in the peaks buffer
  //Now process the peaks and generate the mask:
  uint64_t nPeaksTotal=0;  //Total number of peaks
  for(e=0;e<pr->nEvents;e++){ //Event
    uint64_t nPeaks=*(uint64_t*)(&pr->peaks[peaksBytePos]);
    peaksBytePos+=8;
    //peaksBytePos+=8;//Skip the second one! This is due to the problem in Python.

    nPeaksTotal+=nPeaks;
    for(pk=0;pk<nPeaks;pk++){
      uint16_t p_=*(uint16_t*)(&pr->peaks[peaksBytePos]); //Panel for the current peak
      peaksBytePos+=2;
      uint16_t r_=*(uint16_t*)(&pr->peaks[peaksBytePos]); //Row for the current peak
      peaksBytePos+=2;
      uint16_t c_=*(uint16_t*)(&pr->peaks[peaksBytePos]); //Col for the current peak
      peaksBytePos+=2;
      
      if(p_<0 || p_>=pr->panels){
        // cout<<"ERROR: Peak coordinate out of bounds: Panel="<<p_<<", Valid range: 0,"<<pr->panels-1<<endl; 
        assert(0);
      }
      if(r_<0 || r_>=pr->rows){
        // cout<<"ERROR: Peak coordinate out of bounds: Row="<<r_<<", Valid range: 0,"<<pr->rows-1<<endl;
        assert(0);
      }
      if(c_<0 || c_>=pr->cols){
        // cout<<"ERROR: Peak coordinate out of bounds: Col="<<c_<<", Valid range: 0,"<<pr->cols-1<<endl; 
        assert(0);
      }
      
      for(ri=r_-pr->peakRadius;ri<=r_+pr->peakRadius;ri++){  //ri: row index. Just a temporary variable.
        for(ci=c_-pr->peakRadius;ci<=c_+pr->peakRadius;ci++){  //ci: column index. Just a temporary variable.
          if(ri>=0 && ri<pr->rows && ci>=0 && ci<pr->cols){  //Check whether inside bounds or not
            roiM[calcIdx_4D(e,p_,ri,ci,pr->panels,pr->rows,pr->cols)]=0;
          }
        }
      }
    }
  }
  
  //Save ROI:
  uint64_t roiSavedCount=0;
  for(e=0;e<pr->nEvents;e++){ //Event
    for(p=0;p<pr->panels;p++){ //Panel
      for(r=0;r<pr->rows;r++){ //Row
        for(c=0;c<pr->cols;c++){ //Column
          if(!roiM[calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)]){
            roiData[roiSavedCount]=origData[calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)];
            roiSavedCount++;
          }
          
          //AMG: Replace ROI and RONI pixels with avg values!
          
        }
      }
    }
  }
  
  //Binning:
  for(e=0;e<pr->nEvents;e++){ //Event
    for(p=0;p<pr->panels;p++){  //Panel
      for(r=0;r<pr->binnedRows;r++){ //Row of the binnedData
        for(c=0;c<pr->binnedCols;c++){ //Column of the binnedData
          float sum=0;
          int nPts=0;
          for(br=0;br<pr->binSize;br++) //Bin Row (from origData)
            for(bc=0;bc<pr->binSize;bc++) //Bin Column (from origData)
              if(r*pr->binSize+br<pr->rows && c*pr->binSize+bc<pr->cols){
                // cout<<p<<" "<<r<<" "<<c<<" "<<br<<" "<<bc<<" "<<r*pr->binSize+br<<" "<<c*pr->binSize+bc<<endl;
                sum+=origData[calcIdx_4D(e,p,r*pr->binSize+br,c*pr->binSize+bc,pr->panels,pr->rows,pr->cols)];
                nPts++;
              }
          // cout<<"p:"<<p<<" r:"<<r<<" c:"<<c<<" nPts:"<<nPts<<endl;
          binnedData[calcIdx_4D(e,p,r,c,pr->panels,pr->binnedRows,pr->binnedCols)]=sum/nPts;
        }
      }
    }
  }

  //Additional compression using SZ:    
  SZ_Init(NULL);
  size_t szCompressedSize=0;
  unsigned char* szComp;
   
  switch(pr->szDim){
    case 1:
      // szComp=sz_compress_3D(binnedData, 0, 0, pr->nEvents * pr->panels * pr->binnedRows * pr->binnedCols, pr->tolerance, szCompressedSize); //1D
      szComp=SZ_compress_args(SZ_FLOAT, binnedData, &szCompressedSize, ABS, pr->tolerance, 0, 0, 0, 0,0,0, pr->nEvents * pr->panels * pr->binnedRows * pr->binnedCols);
      break;
    case 2:
      // szComp=sz_compress_3D(binnedData, 0, pr->nEvents * pr->panels * pr->binnedRows, pr->binnedCols, pr->tolerance, szCompressedSize); //2D
      szComp=SZ_compress_args(SZ_FLOAT, binnedData, &szCompressedSize, ABS, pr->tolerance, 0, 0, 0, 0,0, pr->nEvents * pr->panels * pr->binnedRows, pr->binnedCols);
      break;
    case 3:
      // szComp=sz_compress_3D(binnedData, pr->nEvents * pr->panels, pr->binnedRows, pr->binnedCols, pr->tolerance, szCompressedSize); //3D
      szComp=SZ_compress_args(SZ_FLOAT, binnedData, &szCompressedSize, ABS, pr->tolerance, 0, 0, 0, 0, pr->nEvents * pr->panels, pr->binnedRows, pr->binnedCols);
      
      break;
    default:
      // cout<<"Wrong pr->szDim (SZ dimensions)"<<endl;
      // cout<<"pr->szDim:"<<pr->szDim<<endl;
      assert(0);
  }
  
  /*      
  Compressed buffer format: (Types are indicated in parenthesis)
    WRITE: nPeaksTotal(uint64_t) (Total number of peaks in this batch)
    for(e=0;e<pr->nEvents;e++){  (e for "event")
      WRITE: nPeaks[e]  (uint64_t) (Number of peaks in this event)
      for(p=0;p<nPeaks;p++){  (p for "peak")
       nPeaks{
         WRITE: peak[e][p] (uint16_t x 3)
       }
    }
    WRITE: roiSavedCount  (uint64_t) (How many pixels there are in the ROI data.)
       (roiSavedCount is the same # as # of 0's in ROI mask.) 
       (NOTE:0 means affirmative in ROI mask!)
    for(roii=0;roii<roiSavedCount;roii++){  (roii for "ROI data index")
      WRITE: ROI_data[roii]  (float, 32-bit)
    }
    WRITE: szCompressedSize  (uint64_t) (Compressed data size from SZ.)
    WRITE: szComp (unsigned char x SZ_compressed_buffer_size)  (Compressed data from SZ.)
    
    NOTE: Calibration panel is not saved. It should be handled by the user.
    
    SUMMARY:
    nPeaksTotal : 8 bytes : (1 x uint64_t)
    peaks : (8 x pr->nEvents + nPeaksTotal x 3 x 2) bytes : (pr->nEvents x (nPeaks + nPeaks x 3 x uint16_t))
    roiSavedCount : 8 Bytes : (1 x uint64_t)
    ROI_data : roiSavedCount x 4 : roiSavedCount x float 
    szCompressedSize : 8 : uint64_t
    szComp : szComp x 1 : szComp x (unsigned char)
  */
  (*compressedSize)=8+pr->nEvents*8+nPeaksTotal*(2+2+2)+8+roiSavedCount*4+8+szCompressedSize;
  //compressedBuffer=new uint8_t[(*compressedSize)];
  uint8_t * compressedBuffer=(uint8_t*)malloc(*compressedSize);
  uint64_t bytePos;
  
  bytePos=0;
  //*(uint64_t*)(&compressedBuffer[bytePos])=pr->nEvents;
  //bytePos+=8;
  *(uint64_t*)(&compressedBuffer[bytePos])=nPeaksTotal;
  bytePos+=8;
  // cout<<endl;
  // cout<<"COMPRESS:"<<endl;
  // cout<<"nPeaksTotal="<<nPeaksTotal<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("\nCOMPRESS:\n");
  printf("nPeaksTotal=%d\n",nPeaksTotal);
  printf("bytePos=%d\n",bytePos);
  
  peaksBytePos=0;
  for(e=0;e<pr->nEvents;e++){
    uint64_t nPeaks=*(uint64_t*)(&pr->peaks[peaksBytePos]);
    peaksBytePos+=8;
    //peaksBytePos+=8;//Skip the second one. This is due to the error in Python!
    
    *(uint64_t*)(&compressedBuffer[bytePos])=nPeaks;
    bytePos+=8;
    for(pk=0;pk<nPeaks;pk++){
      *(uint16_t*)(&compressedBuffer[bytePos])=*(uint16_t*)(&pr->peaks[peaksBytePos]); //Panel for the current peak
      bytePos+=2;
      peaksBytePos+=2;
      *(uint16_t*)(&compressedBuffer[bytePos])=*(uint16_t*)(&pr->peaks[peaksBytePos]); //Row for the current peak
      bytePos+=2;
      peaksBytePos+=2;      
      *(uint16_t*)(&compressedBuffer[bytePos])=*(uint16_t*)(&pr->peaks[peaksBytePos]); //Column for the current peak
      bytePos+=2;
      peaksBytePos+=2;
    }
  }
  // cout<<"peaks"<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("peaks\n");
  printf("bytePos=%d\n",bytePos);

  *(uint64_t*)(&compressedBuffer[bytePos])=roiSavedCount;
  bytePos+=8;
  // cout<<"roiSavedCount="<<roiSavedCount<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  // cout<<"roiData"<<endl;
  printf("roiSavedCount=%d\n",roiSavedCount);
  printf("bytePos=%d\n",bytePos);
  printf("roiData\n");
  for(roii=0;roii<roiSavedCount;roii++){
    *(float*)(&compressedBuffer[bytePos])=roiData[roii];
    // cout<<roiData[roii]<<",";
    bytePos+=4;
  }
  // cout<<"bytePos="<<bytePos<<endl;
  printf("bytePos=%d\n",bytePos);
  *(uint64_t*)(&compressedBuffer[bytePos])=szCompressedSize;
  bytePos+=8;
  // cout<<"szCompressedSize="<<szCompressedSize<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("szCompressedSize=%d\n",szCompressedSize);
  printf("bytePos=%d\n",bytePos);
  for(bi=0;bi<szCompressedSize;bi++){  //bi for "byte index"
    *(unsigned char*)(&compressedBuffer[bytePos])=szComp[bi];
    bytePos+=1;
  }
  // cout<<"szComp"<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("szComp\n");
  printf("bytePos=%d\n",bytePos);
  
  if(bytePos!=(*compressedSize)){
    // cout<<"ERROR: bytePos = "<<bytePos<<" != "<<(*compressedSize)<<" = compressedSize"<<endl;
    assert(0);
  }
  SZ_Finalize();
  
  free(roiM);
  free(roiData);
  free(binnedData);
  // delete [] roiM;
  // delete [] roiData;
  // delete [] binnedData;
  
  return compressedBuffer;
}

void* exafelSZ_Decompress(unsigned char*_compressedBuffer,
                         void *_pr,
                         size_t compressedSize
                         )
{ 
  uint8_t *compressedBuffer=(uint8_t *)_compressedBuffer;
  exafelSZ_params *pr=(exafelSZ_params *)_pr;
  exafelSZ_params_process(pr); 
  exafelSZ_params_checkDecomp(pr); 
  
  float *decompressedBuffer=(float*)malloc(pr->nEvents*pr->panels*pr->rows*pr->cols*sizeof(float));
  
  uint8_t *roiM=(uint8_t*)malloc(pr->nEvents*pr->panels*pr->rows*pr->cols);
  int e,p,r,c,pk,ri,ci,br,bc;
  
  
  /*
  Compressed Data Layout:
  nPeaksTotal : 8 bytes : (1 x uint64_t)
  peaks : (8 x pr->nEvents + nPeaksTotal x 3 x 2) bytes : (pr->nEvents x (nPeaks + nPeaks x 3 x uint16_t))
  roiSavedCount : 8 Bytes : (1 x uint64_t)
  ROI_data : roiSavedCount x 4 : roiSavedCount x float 
  szCompressedSize : 8 : uint64_t
  szComp : szComp x 1 : szComp x (unsigned char)
  */
  uint64_t bytePos=0;
  uint64_t nPeaksTotal=*(uint64_t*)(&compressedBuffer[bytePos]);
  bytePos += 8; 
  // cout<<endl;
  // cout<<"DECOMPRESS:"<<endl;
  // cout<<"nPeaksTotal="<<nPeaksTotal<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("\nDECOMPRESS:\n");
  printf("nPeaksTotal=%d\n",nPeaksTotal);
  printf("bytePos=%d\n",bytePos);
  
  uint8_t *peaks=(uint8_t*)(&compressedBuffer[bytePos]);
  bytePos += (8 * pr->nEvents + nPeaksTotal * 3 * 2);
  // cout<<"peaks"<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("peaks\n");
  printf("bytePos=%d\n",bytePos);
  
  uint64_t roiSavedCount=*(uint64_t*)(&compressedBuffer[bytePos]);
  bytePos+=8;
  // cout<<"roiSavedCount="<<roiSavedCount<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("roiSavedCount=%d\n",roiSavedCount);
  printf("bytePos=%d\n",bytePos);
  
  // cout<<"roiData"<<endl;
  float *roiData=(float*)(&compressedBuffer[bytePos]);
  bytePos+=(roiSavedCount*4);
  // for(uint64_t roii=0;roii<roiSavedCount;roii++){
    // cout<<roiData[roii]<<",";
  // }
  // cout<<"bytePos="<<bytePos<<endl;
  printf("bytePos=%d\n",bytePos);
  
  uint64_t szCompressedSize=*(uint64_t*)(&compressedBuffer[bytePos]);
  bytePos+=8;
  // cout<<"szCompressedSize="<<szCompressedSize<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  printf("szCompressedSize=%d\n",szCompressedSize);
  printf("bytePos=%d\n",bytePos);
  
  unsigned char *szComp=(unsigned char*)(&compressedBuffer[bytePos]);
  bytePos+=szCompressedSize;
  // cout<<"szComp"<<endl;
  // cout<<"bytePos="<<bytePos<<endl;
  // cout<<endl;
  printf("szComp\n");
  printf("bytePos=%d\n\n",bytePos);
  
  //We should have inputs ready by now. Now process them:
  
  //Generate the ROI mask: NOTE: 0 means affirmative in ROI mask! This comes from the python scripts!
  //First, initialize with calibration panel:
  for(e=0;e<pr->nEvents;e++){ //Event
    for(p=0;p<pr->panels;p++){ //Panel
      for(r=0;r<pr->rows;r++){ //Row
        for(c=0;c<pr->cols;c++){ //Column
          if(calcIdx_2D(r,c,pr->cols)<0 ||calcIdx_2D(r,c,pr->cols)>=pr->rows*pr->cols){
            // cout<<"ERROR: calcIdx(r,c,pr->cols) = calcIdx("<<r<<","<<c<<","<<pr->cols<<") = "<<calcIdx(r,c,pr->cols)<<endl;
            // cout<<"       NOT in the correct range: [0,"<<pr->rows*pr->cols-1<<"]"<<endl;
            assert(0);
          }
          if(calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)<0 ||calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)>=pr->nEvents*pr->panels*pr->rows*pr->cols){
            // cout<<"ERROR: calcIdx(e,p,r,c,pr->panels,pr->rows,pr->cols) = calcIdx("<<e<<","<<p<<","<<r<<","<<c<<","<<pr->panels<<","<<pr->rows<<","<<pr->cols<<") = "<<calcIdx(e,p,r,c,pr->panels,pr->rows,pr->cols)<<endl;
            // cout<<"       NOT in the correct range: [0,"<<pr->nEvents*pr->panels*pr->rows*pr->cols-1<<"]"<<endl;
            assert(0);
          }
          roiM[calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)]=pr->calibPanel[calcIdx_2D(r,c,pr->cols)];
        }
      }
    }
  }
  uint64_t peaksBytePos=0; //Position in the peaks buffer
  //Now process the peaks and generate the mask:
  for(e=0;e<pr->nEvents;e++){ //Event
    uint64_t nPeaks=*(uint64_t*)(&peaks[peaksBytePos]);
    peaksBytePos+=8;
    
    for(pk=0;pk<nPeaks;pk++){
      uint16_t p_=*(uint16_t*)(&peaks[peaksBytePos]); //Panel for the current peak
      peaksBytePos+=2;
      uint16_t r_=*(uint16_t*)(&peaks[peaksBytePos]); //Row for the current peak
      peaksBytePos+=2;
      uint16_t c_=*(uint16_t*)(&peaks[peaksBytePos]); //Col for the current peak
      peaksBytePos+=2;
      
      if(p_<0 || p_>=pr->panels){
        // cout<<"ERROR: Peak coordinate out of bounds: Panel="<<p_<<", Valid range: 0,"<<pr->panels-1<<endl; 
        assert(0);
      }
      if(r_<0 || r_>=pr->rows){
        // cout<<"ERROR: Peak coordinate out of bounds: Row="<<r_<<", Valid range: 0,"<<pr->rows-1<<endl; 
        assert(0);
      }
      if(c_<0 || c_>=pr->cols){
        // cout<<"ERROR: Peak coordinate out of bounds: Col="<<c_<<", Valid range: 0,"<<pr->cols-1<<endl; 
        assert(0);
      }
      
      for(ri=r_-pr->peakRadius;ri<=r_+pr->peakRadius;ri++){  //ri: row index. Just a temporary variable.
        for(ci=c_-pr->peakRadius;ci<=c_+pr->peakRadius;ci++){  //ci: column index. Just a temporary variable.
          if(ri>=0 && ri<pr->rows && ci>=0 && ci<pr->cols){  //Check whether inside bounds or not
            roiM[calcIdx_4D(e,p_,ri,ci,pr->panels,pr->rows,pr->cols)]=0;
          }
        }
      }
    }
  }
  
  //De-compress using SZ:
  SZ_Init(NULL);
  float* szDecomp;
  size_t _szCompressedSize=szCompressedSize;
  switch(pr->szDim){
    case 1:
      // szDecomp=sz_decompress_3D<float>(szComp,0,0,pr->nEvents * pr->panels * pr->binnedRows * pr->binnedCols,pr->tolerance,szCompressedSize); //1D
      szDecomp=SZ_decompress(SZ_FLOAT,szComp,_szCompressedSize,0,0,0,0, pr->nEvents * pr->panels * pr->binnedRows * pr->binnedCols);
      break;
    case 2:
      // szDecomp=sz_decompress_3D<float>(szComp,0,pr->nEvents * pr->panels * pr->binnedRows, pr->binnedCols,pr->tolerance,szCompressedSize); //2D
      szDecomp=SZ_decompress(SZ_FLOAT,szComp,_szCompressedSize,0,0,0, pr->nEvents * pr->panels * pr->binnedRows, pr->binnedCols);
      break;
    case 3:
      // szDecomp=sz_decompress_3D<float>(szComp,pr->nEvents * pr->panels, pr->binnedRows, pr->binnedCols,pr->tolerance,szCompressedSize); //3D
      szDecomp=SZ_decompress(SZ_FLOAT,szComp,_szCompressedSize,0,0,pr->nEvents * pr->panels, pr->binnedRows, pr->binnedCols);
      break;
    default:
      // cout<<"Wrong pr->szDim (SZ dimensions)"<<endl;
      // cout<<"pr->szDim:"<<pr->szDim<<endl;
      assert(0);
  }
  // double max_err = 0;
  // for(int i=0; i<pr->nEvents * pr->panels * pr->binnedRows * pr->binnedCols; i++){
    // double err = fabs(szDecomp[i]-binnedData[i]);
    // if(err > max_err) max_err = err;
  // }
  // cout << "Max err = \t\t\t" << max_err << endl;
  
  
  //De-binning:
  for(e=0;e<pr->nEvents;e++)//Event
    for(p=0;p<pr->panels;p++)  //Panel
      for(r=0;r<pr->binnedRows;r++) //Row of the binnedData
        for(c=0;c<pr->binnedCols;c++) //Column of the binnedData
            for(br=0;br<pr->binSize;br++) //Bin Row (from origData)
              for(bc=0;bc<pr->binSize;bc++) //Bin Column (from origData)
                if(r*pr->binSize+br<pr->rows && c*pr->binSize+bc<pr->cols){
                  decompressedBuffer[calcIdx_4D(e,p,r*pr->binSize+br,c*pr->binSize+bc,pr->panels,pr->rows,pr->cols)] = szDecomp[calcIdx_4D(e,p,r,c,pr->panels,pr->binnedRows,pr->binnedCols)];
                }
  //Restore ROI:
  uint64_t current=0;
  for(e=0;e<pr->nEvents;e++)//Event
    for(p=0;p<pr->panels;p++)  //Panel
      for(r=0;r<pr->rows;r++) //Row of the binnedData
        for(c=0;c<pr->cols;c++) //Column of the binnedData
          if(!roiM[calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)]){
            decompressedBuffer[calcIdx_4D(e,p,r,c,pr->panels,pr->rows,pr->cols)]=roiData[current];
            current++;
          }
  SZ_Finalize();
  // delete [] roiM;
  free(roiM);
  
  return ((void*)decompressedBuffer);
}

