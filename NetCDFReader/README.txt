#COMPILATION

1. Install NetCDF.
2. Modify Makefile.linux or Makefile.osx: NETCDFPATH 
3. Linux: make -f Makefile.linux
   OSX:   make -f Makefile.osx
4. Add SZPATH/lib and NetCDFPATH/lib to the environment variable LD_LIBRARY_PATH
