#include "filters.h"
#include <Rcpp.h>

using namespace Rcpp;

// Use C++ 11
// ==========

// [[Rcpp::plugins("cpp11")]]

// Forward declarations
// ====================

// NumericVector wavelet_filter(String filter);
// NumericVector scaling_filter(String filter);

// Main function
// =============

//' MODWT Wavelet Coefficients
//'
//' Compute the wavelet coefficients. This function should not be called by end users.
//' If you want to obtain the MODWT wavelet coefficient for a given filter and
//' decomposition level, please use `modwt_wavelet_coefs`.
//'
//' @param X An (N x 1) matrix or a vector
//' @param filter A character string indicating the desired filter
//' @param j The decomposition level \[integer\]
//' @details This is an internal function used by `atrous_dwt`.
//' The argument `filter` can take one of the following values:
//'
//' c('bl7', 'bl9', 'bl10',
//' 'beyl',
//' 'coif1', 'coif2', 'coif3', 'coif4', 'coif5',
//' 'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db11', 'db12',
//' 'db13', 'db14', 'db15', 'db16', 'db17', 'db18', 'db19', 'db20', 'db21', 'db22', 'db23',
//' 'db24', 'db25', 'db26', 'db27', 'db28', 'db29', 'db30', 'db31', 'db32', 'db33',
//' 'db34', 'db35', 'db36', 'db37', 'db38', 'db39', 'db40', 'db41', 'db42', 'db43', 'db44', 'db45',
//' 'fk4', 'fk6', 'fk8', 'fk14', 'fk18', 'fk22',
//' 'han2_3', 'han3_3', 'han4_5', 'han5_5',
//' 'dmey',
//' 'mb4_2', 'mb8_2', 'mb8_3', 'mb8_4', 'mb10_3', 'mb12_3', 'mb14_3', 'mb16_3', 'mb18_3', 'mb24_3', 'mb32_3',
//' 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'sym9', 'sym10', 'sym11', 'sym12', 'sym13', 'sym14',
//' 'sym15', 'sym16', 'sym17', 'sym18', 'sym19', 'sym20', 'sym21', 'sym22', 'sym23', 'sym24', 'sym25', 'sym26', 'sym27',
//' 'sym28', 'sym29', 'sym30', 'sym31', 'sym32', 'sym33', 'sym34', 'sym35', 'sym36', 'sym37', 'sym38', 'sym39', 'sym40',
//' 'sym41', 'sym42', 'sym43', 'sym44', 'sym45',
//' 'vaid',
//' 'la8', 'la10', 'la12', 'la14', 'la16', 'la18', 'la20')
//' @return (N x 1) matrix of wavelet coefficients
//' @references
//' Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge
//' University Press.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix wavelet_coefs(NumericVector X, String filter, int j)
{

    int N = X.size();
    NumericVector h = wavelet_filter(filter); // scaled by 1/sqrt(2)
    int L = h.size();

    // ERROR CHECKING

    // Check maximum decomposition level j
    if (j > log(((N-1)/(L-1)) +1)/log(2)) {
      stop("Decomposition level j is too large!");
    }

    // convolve time series with scaling filter to calculate scaling coefficients
    // at this particular decomposition level

    // preallocae space
    NumericMatrix d(N, 1);
    int pos;
    double tmp;

    for (int k = N; k >= pow(2, j + 1) + 1; --k)
    {
        tmp = 0;
        for (int i = 1; i <= L; ++i)
        {

            pos = k - pow(2, j - 1) * (i - 1); // time position

            // catch cases where time position <= 0 and wrap the position based
            // on the circularity assumption (i.e., mod(t,N) = t)

            if (pos <= 0)
            {

                pos += N; // circularity assumption
            }

            tmp += h(i - 1) * X(pos - 1);
            //d(k-1, 0) += h(i - 1) * X(pos - 1);
        }

        d(k-1, 0) = tmp;
    }

    return d;

}
