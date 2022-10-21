#include "filters.h"
#include <Rcpp.h>

using namespace Rcpp;

// Use C++ 11
// ==========

// [[Rcpp::plugins("cpp11")]]

// Forward declarations
// ====================

// NumericVector scaling_filter(String wavelet);

// Main function
// =============

//' Compute Scaling Coefficients
//'
//' @param X A numeric vector
//' @param wavelet A character string indicating the scaling filter desired
//' @param j The decomposition level \[integer\]
//' @return Matrix of scaling coefficients
//' @references
//' Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge
//' University Press.
// [[Rcpp::export]]
NumericMatrix scaling_coefs(NumericVector X, String wavelet, int j)
{

    int N = X.size();
    NumericVector g = scaling_filter(wavelet);
    int L = g.size();

    //' convolve time series with scaling filter to calculate scaling coefficients
    //' at this particular decomposition level

    // preallocae space
    NumericMatrix c(N, 1);
    int pos;
    double tmp;

    for (int k = N; k >= pow(2, j + 1) + 1; --k)
    {
        tmp = 0;
        for (int i = 1; i <= L; ++i)
        {

            pos = k - pow(2, j - 1) * (i - 1); // time position

            //' catch cases where time position <= 0 and wrap the position based
            //' on the circularity assumption (i.e., mod(t,N) = t)

            if (pos <= 0)
            {

                pos += N; // circularity assumption
            }

            tmp += g(i - 1) * X(pos - 1);
        }

        c[k-1] = tmp;
    }

    return c;

}
