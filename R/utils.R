.onUnload <- function (libpath) {
  library.dynam.unload("fastWavelets", libpath)
}

#' Number of Boundary Coefficients
#'
#' This function calculates the number of boundary coefficients for a
#' particular wavelet/scaling filter and decomposition level.
#'
#' @param filter Filter name (see Details below) \[string\]
#' @param J Decomposition level \[integer\]
#' @details
#' The argument `filter` can take one of the following values:
#'
#' c('bl7', 'bl9', 'bl10',
#' 'beyl',
#' 'coif1', 'coif2', 'coif3', 'coif4', 'coif5',
#' 'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db11', 'db12',
#' 'db13', 'db14', 'db15', 'db16', 'db17', 'db18', 'db19', 'db20', 'db21', 'db22', 'db23',
#' 'db24', 'db25', 'db26', 'db27', 'db28', 'db29', 'db30', 'db31', 'db32', 'db33',
#' 'db34', 'db35', 'db36', 'db37', 'db38', 'db39', 'db40', 'db41', 'db42', 'db43', 'db44', 'db45',
#' 'fk4', 'fk6', 'fk8', 'fk14', 'fk18', 'fk22',
#' 'han2_3', 'han3_3', 'han4_5', 'han5_5',
#' 'dmey',
#' 'mb4_2', 'mb8_2', 'mb8_3', 'mb8_4', 'mb10_3', 'mb12_3', 'mb14_3', 'mb16_3', 'mb18_3', 'mb24_3', 'mb32_3',
#' 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'sym9', 'sym10', 'sym11', 'sym12', 'sym13', 'sym14',
#' 'sym15', 'sym16', 'sym17', 'sym18', 'sym19', 'sym20', 'sym21', 'sym22', 'sym23', 'sym24', 'sym25', 'sym26', 'sym27',
#' 'sym28', 'sym29', 'sym30', 'sym31', 'sym32', 'sym33', 'sym34', 'sym35', 'sym36', 'sym37', 'sym38', 'sym39', 'sym40',
#' 'sym41', 'sym42', 'sym43', 'sym44', 'sym45',
#' 'vaid',
#' 'la8', 'la10', 'la12', 'la14', 'la16', 'la18', 'la20')
#' @return Number of boundary coefficients \[integer\]
#' @references
#'
#' M. Basta (2014),Additive Decomposition and Boundary Conditions in Wavelet-Based
#' Forecasting Approaches, Acta Oeconomica Pragensia, 2, pp. 48-70.
#'
#' Quilty, J., &amp; Adamowski, J. (2018). Addressing the incorrect usage of wavelet-based
#' hydrological and water resources forecasting models for real-world applications with best
#' practices and a new forecasting framework. Journal of Hydrology, 563, 336–353.
#' https://doi.org/10.1016/j.jhydrol.2018.05.003
#'
#' Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge
#' University Press.
#'
#' @examples
#' J <- 4 # decomposition level
#' nbc <- n_boundary_coefs('db1', J) # number of boundary-effected coefficients at decomp_level J
#' @export
n_boundary_coefs <- function(filter, J){

  g = scaling_filter(filter)
  L = length(g)
  nbc = ((2^J) - 1) * (L - 1)

  # n_boundary_coefs <- nbc
  return(nbc)

}

#' @title Check Matrix X is (N x 1)
#' @description `shape_check` checks whether `X` is a matrix representing
#' a column vector (i.e., a matrix with 1 column). If not, `shape_check` attempts
#' to coerce the user provided `X` to a matrix with 1 column. If this cannot be done,
#' an error is raised.
#' @param X Object to check and (if possible) coerce to a single column matrix
#' @return An (N x 1) matrix
#' @details This is a utility function written to check the input `X` for the
#' functions `atrous_dwt` and `mo_dwt`.
#' @keywords internal
shape_check <- function(X) {

  if (is.vector(X)) {
    X <- t(t(X))
  } else if (!is.matrix(X)) {
    stop("X must be an (N x 1) matrix or a vector.")
  } else if (length(X) == ncol(X)) { # row vector
    X <- t(X)
  } else if (nrow(X) > 1 && ncol(X) > 1) {
    stop("X must be an (N x 1) matrix or a vector.")
  } else if (nrow(X) == 1 && ncol(X) == 1) {
    stop("Matrix X should have more than 1 element.")
  }

  return(X)

}

#' @title Periodize a filter
#' @param filter_vector The filter vector that is outputted from r_wavelet_filter or r_scaling_filter
#' @param periodize_to The length of the time series that is to be decomposed
#' @details Output of the function is a zero padded version of the original filter, padded at the end.
#' @return Periodized filter
#' @export
periodize_filter <- function(filter_vector, periodize_to) {
  L <- length(filter_vector)
  if (L >= periodize_to) {stop("This function currently only supports L < periodize_to")}
  periodized_filter_vector <- utils::head(c(filter_vector, rep(0, periodize_to)), n = periodize_to)
  return(periodized_filter_vector)
}

#' @title Circular Convolution
#' @description Perform a circular convolution between a filter and a time series.
#' @param a Filter (typically, the wavelet or scaling filter)
#' @param b A vector
#' @return Result of circular convolution, length is the same as `length(b)`
#' @details The convolution theorem tells us that the discrete Fourier transform
#' can be leveraged to compute a circular convolution, which is what this function does.
#' @seealso
#'  \code{\link[stats]{fft}}
#' @rdname circular_convolution
#' @export
circular_convolution <- function(a, b) {
  # a <- utils::head(c(a,0*b), n = length(b))
  a = periodize_filter(a, length(b))
  A = stats::fft(a)
  B = stats::fft(b)
  C = A*B
  c = stats::fft(C, inverse = TRUE) / length(C)
  return(Re(c))
}





