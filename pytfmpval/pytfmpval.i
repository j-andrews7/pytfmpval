%module pytfmpval
%{
#include "../src/Matrix.h"
#define SWIG_FILE_WITH_INIT
%}

%include "cpointer.i"
%include "std_string.i"
%include "std_vector.i"
%include "typemaps.i"
%include "../src/Matrix.h"

%pointer_class(double, doublep)
%pointer_class(int, intp)

%nodefaultdtor Matrix;
