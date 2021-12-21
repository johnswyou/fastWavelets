#ifndef FILTERS_H
#define FILTERS_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector wavelet_filter(String wavelet);
NumericVector scaling_filter(String wavelet);

#endif
