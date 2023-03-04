#ifndef FILTERS_H
#define FILTERS_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector wavelet_filter(String filter);
NumericVector scaling_filter(String filter);

#endif
