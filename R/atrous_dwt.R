#' A Trous Discrete Wavelet Transform
#'
#' This function calculates the wavelet and scaling coefficients of the a trous (AT)
#' version of the Discrete Wavelet Transform (DWT).
#'
#' @param X An (N x 1) matrix or a vector
#' @param wavelet Scaling filter name (see Details below) (string)
#' @param decomp_level Decomposition level (integer, 1 < decomp_level < N/2)
#' @return Wavelet and scaling coefficients (N x J+1)
#' (wavelet coefficients in first J columns of returned matrix)
#' @references
#'
#' Benaouda, D., F. Murtagh, J. L. Starck, and O. Renaud (2006),
#' Wavelet-based nonlinear multiscale decomposition model for
#' electricity load forecasting, Neurocomputing,
#' doi:10.1016/j.neucom.2006.04.005.
#'
#' Maheswaran, R., and R. Khosa (2012), Comparative study of different
#' wavelets for hydrologic forecasting, Comput. Geosci.,
#' doi:10.1016/j.cageo.2011.12.015.
#'
#' Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge
#' University Press.
#'
#' @details
#'
#' The argument `wavelet` can take one of the following values:
#'
#' `c("haar", "d1", "sym1", "bior1.1", "rbio1.1",
#' "d2", "sym2", "d3", "sym3", "d4", "d5", "d6", "d7", "d8", "d9", "d10", "d11",
#' "sym4", "sym5", "sym6", "sym7", "sym8", "sym9", "sym10",
#' "coif1", "coif2", "coif3", "coif4", "coif5",
#' "bior1.3", "bior1.5", "bior2.2", "bior2.4", "bior2.6", "bior2.8", "bior3.1", "bior3.3",
#' "bior3.5", "bior3.7", "bior3.9", "bior4.4", "bior5.5", "bior6.8",
#' "rbio1.3", "rbio1.5", "rbio2.2", "rbio2.4", "rbio2.6", "rbio2.8", "rbio3.1", "rbio3.3",
#' "rbio3.5", "rbio3.7", "rbio3.9", "rbio4.4", "rbio5.5", "rbio6.8",
#' "la8", "la10", "la12", "la14", "la16", "la18", "la20",
#' "bl14", "bl18", "bl20",
#' "fk4", "fk6", "fk8", "fk14", "fk18", "fk22",
#' "b3spline")`
#'
#' @examples
#' N <- 1000 #  number of time series points
#' J <- 4 # decomposition level
#' wavelet <- 'coif1' # wavelet filter
#' X <- matrix(rnorm(N),N,1)
#' W <- atrous_dwt(X,wavelet,J)
#' Xr <- as.matrix(rowSums(W)) # reconstruct time series
#' mse_r <- mean( (X - Xr)^2) # confirm additive reconstruction
#' plot.ts(W) # plot wavelet and scaling coefficients
#' @export
atrous_dwt <- function(X,wavelet,decomp_level){

  X <- shape_check(X)

  N = length(X)
  # g = scaling_filter(wavelet)
  # L = length(g)
  J = decomp_level

  # use periodic extension on time series in order to calculate wavelet and scaling
  # coefficients at boundary of time series, i.e., the first ((2^J) - 1) * (L - 1) + 1
  # time series observations
  x = rbind(X,X) # periodic extension
  rm(X)

  # perform AT algorithm on extended time series

  # preallocate space
  v = matrix(0,2*N,J)
  w = matrix(0,2*N,J)

  for(j in seq(1,J)){

    # calculate scaling coefficients

    v[,j] = scaling_coefs(x,wavelet,j) # modified version of Eq. 5 in Maheswaran and Khosa [2012]

    # calculate wavelet coefficients

    w[,j] = x - v[,j] # see Eq. 6 in Maheswaran and Khosa [2012]

    # calculate residual information

    x = v[,j]

  }

  # combine wavelet and scaling coefficients in a single matrix and only store the
  # last N observations

  W = cbind(w[seq(N+1,2*N),], v[seq(N+1,2*N),J])

  # atrous_dwt <- W
  return(W)

}
