#include "filters.h"
#include <Rcpp.h>

using namespace Rcpp;

// Use C++ 11
// ==========

// [[Rcpp::plugins("cpp11")]]

// Forward declarations
// ====================

// NumericVector wavelet_filter(String wavelet);
// NumericVector scaling_filter(String wavelet);

// Main function
// =============

//' Computes wavelet coefficients
//'
//' @param X A vector
//' @param wavelet A character string indicating the wavelet filter desired
//' @param j The decomposition level
//' @return Matrix of wavelet coefficients
// [[Rcpp::export]]
NumericMatrix wavelet_coefs(NumericVector X, String wavelet, int j)
{

    int N = X.size();
    NumericVector h = wavelet_filter(wavelet);
    int L = h.size();

    //' convolve time series with scaling filter to calculate scaling coefficients
    //' at this particular decomposition level

    // preallocae space
    NumericMatrix d(N, 1);
    int pos;
    // double tmp;

    for (int k = N; k >= pow(2, j + 1) + 1; --k)
    {
        // tmp = 0;
        for (int i = 1; i <= L; ++i)
        {

            pos = k - pow(2, j - 1) * (i - 1); // time position

            //' catch cases where time position <= 0 and wrap the position based
            //' on the circularity assumption (i.e., mod(t,N) = t)

            if (pos <= 0)
            {

                pos += N; // circularity assumption
            }

            // tmp += h(i - 1) * X(pos - 1);
            d(k-1, 0) += h(i - 1) * X(pos - 1);
        }

        // d(k-1, 0) = tmp;
    }

    return d;

}
