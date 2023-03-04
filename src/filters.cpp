#include "filters.h"
#include <Rcpp.h>
using namespace Rcpp;

// Use C++ 11
// ==========

// [[Rcpp::plugins("cpp11")]]

//' Maximal Overlap Discrete Wavelet Transform Wavelet Filter
//'
//' Compute the MODWT wavelet filter, which is just the regular wavelet
//' filter, divided by `sqrt(2)`. This function is used for the pyramid
//' algorithm used in `modwt`.
//'
//' @param filter A character string indicating the wavelet filter desired
//' @return Wavelet filter vector (a numeric vector)
//' @references
//' Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge
//' University Press.
//'
//' Wasilewski, F. (2008). Wavelet browser by pywavelets. Wavelet Properties Browser.
//' Retrieved November 17, 2022, from http://wavelets.pybytes.com/
//'
//' Gregory R. Lee, Ralf Gommers, Filip Wasilewski, Kai Wohlfahrt, Aaron O’Leary (2019).
//' PyWavelets: A Python package for wavelet analysis. Journal of Open Source Software,
//' 4(36), 1237, https://doi.org/10.21105/joss.01237.
//'
//' Olhede, S., &amp; Walden, A. T. (2004). The Hilbert spectrum via wavelet projections.
//' Proceedings of the Royal Society of London. Series A: Mathematical, Physical and
//' Engineering Sciences, 460(2044), 955–975. https://doi.org/10.1098/rspa.2003.1199
//'
//' Maheswaran, R., &amp; Khosa, R. (2012). Comparative study of different wavelets for
//' hydrologic forecasting. Computers &amp; Geosciences, 46, 284–295.
//' https://doi.org/10.1016/j.cageo.2011.12.015
//'
//' Daubechies, Ingrid. Ten Lectures on Wavelets. Society for Industrial and Applied Mathematics, 1992.
//'
//' Morris, Joel M, and Ravindra Peravali. “Minimum-Bandwidth Discrete-Time Wavelets.” Signal Processing 76,
//' no. 2 (July 1999): 181–93. https://doi.org/10.1016/S0165-1684(99)00007-9.
//'
//' Doroslovački, M.L. “On the Least Asymmetric Wavelets.” IEEE Transactions on Signal Processing 46,
//' no. 4 (April 1998): 1125–30. https://doi.org/10.1109/78.668562.
//'
//' Han, Bin. “Wavelet Filter Banks.” In Framelets and Wavelets: Algorithms, Analysis, and Applications,
//' 92–98. Applied and Numerical Harmonic Analysis. Cham, Switzerland: Birkhäuser, 2017.
//' https://doi.org/10.1007/978-3-319-68530-4_2.
//'
//' Vaidyanathan, P. P. & Hoang, P.-Q. Lattice structures for optimal design and robust implementation
//' of two-channel perfect-reconstruction QMF banks. IEEE Trans. Acoust. 36, 81–94 (1988).
//'
//' Wickerhauser, M. V. Adapted Wavelet Analysis from Theory to Software.
//' Adapted Wavelet Analysis from Theory to Software (A.K. Peters, 1994).
//' @keywords internal
// [[Rcpp::export]]
NumericVector wavelet_filter(String filter){

  // Make the self-defined R function "r_wavelet_filter" available in the
  // C++ code. It is called "cpp_wavelet_filter".
  Function cpp_wavelet_filter("r_wavelet_filter");
  // Environment pkg = Environment::namespace_env("fastWavelets");
  // Function cpp_wavelet_filter = pkg["r_wavelet_filter"];

  // Use function "cpp_wavelet_filter"
  NumericVector result = cpp_wavelet_filter(filter, Named("modwt", true));

  return result;
}

//' Maximal Overlap Discrete Wavelet Transform Scaling Filter
//'
//' Compute the MODWT scaling filter, which is just the regular wavelet
//' filter, divided by `sqrt(2)`. This function is used for the pyramid algorithm
//' used in `modwt`.
//'
//' @param filter A character string indicating the scaling filter desired
//' @return Scaling filter vector (a numeric vector)
//' @references
//' Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge
//' University Press.
//'
//' Wasilewski, F. (2008). Wavelet browser by pywavelets. Wavelet Properties Browser.
//' Retrieved November 17, 2022, from http://wavelets.pybytes.com/
//'
//' Gregory R. Lee, Ralf Gommers, Filip Wasilewski, Kai Wohlfahrt, Aaron O’Leary (2019).
//' PyWavelets: A Python package for wavelet analysis. Journal of Open Source Software,
//' 4(36), 1237, https://doi.org/10.21105/joss.01237.
//'
//' Olhede, S., &amp; Walden, A. T. (2004). The Hilbert spectrum via wavelet projections.
//' Proceedings of the Royal Society of London. Series A: Mathematical, Physical and
//' Engineering Sciences, 460(2044), 955–975. https://doi.org/10.1098/rspa.2003.1199
//'
//' Maheswaran, R., &amp; Khosa, R. (2012). Comparative study of different wavelets for
//' hydrologic forecasting. Computers &amp; Geosciences, 46, 284–295.
//' https://doi.org/10.1016/j.cageo.2011.12.015
//'
//' Daubechies, Ingrid. Ten Lectures on Wavelets. Society for Industrial and Applied Mathematics, 1992.
//'
//' Morris, Joel M, and Ravindra Peravali. “Minimum-Bandwidth Discrete-Time Wavelets.” Signal Processing 76,
//' no. 2 (July 1999): 181–93. https://doi.org/10.1016/S0165-1684(99)00007-9.
//'
//' Doroslovački, M.L. “On the Least Asymmetric Wavelets.” IEEE Transactions on Signal Processing 46,
//' no. 4 (April 1998): 1125–30. https://doi.org/10.1109/78.668562.
//'
//' Han, Bin. “Wavelet Filter Banks.” In Framelets and Wavelets: Algorithms, Analysis, and Applications,
//' 92–98. Applied and Numerical Harmonic Analysis. Cham, Switzerland: Birkhäuser, 2017.
//' https://doi.org/10.1007/978-3-319-68530-4_2.
//'
//' Vaidyanathan, P. P. & Hoang, P.-Q. Lattice structures for optimal design and robust implementation
//' of two-channel perfect-reconstruction QMF banks. IEEE Trans. Acoust. 36, 81–94 (1988).
//'
//' Wickerhauser, M. V. Adapted Wavelet Analysis from Theory to Software.
//' Adapted Wavelet Analysis from Theory to Software (A.K. Peters, 1994).
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector scaling_filter(String filter){

  // Make the self-defined R function "r_scaling_filter" available in the
  // C++ code. It is called "cpp_scaling_filter".
  Function cpp_scaling_filter("r_scaling_filter");
  // Environment pkg = Environment::namespace_env("fastWavelets");
  // Function cpp_scaling_filter = pkg["r_scaling_filter"];

  // Use function "cpp_scaling_filter"
  NumericVector result = cpp_scaling_filter(filter, Named("modwt", true));

  return result;
}

