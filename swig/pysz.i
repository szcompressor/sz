%module pysz

%{
#include "sz.h"
%}

%ignore confparams_cpr;
%ignore confparams_dec;
%ignore exe_params;
%ignore sz_varset;
%ignore multisteps;
%ignore sz_tsc;
%include "sz.h"
