#' A Trous Discrete Wavelet Transform
#'
#' This function calculates the wavelet and scaling coefficients of the a trous (AT)
#' version of the Discrete Wavelet Transform (DWT).
#'
#' @param X Time series (N x 1)
#' @param wavelet Scaling filter name (string)
#' @param decomp_level Decomposition level (1 < integer < N/2)
#' @return W - wavelet and scaling coefficients (N x J+1) (wavelet coefs in first J columns)
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
#' @examples
#' \dontrun{
#' N = 1000 #  number of time series points
#' J = 4 # decomposition level
#' wavelet = 'coif1' # wavelet filter
#' X = matrix(rnorm(N),N,1)
#' W = atrous_dwt(X,wavelet,J)
#' Xr = as.matrix(rowSums(W)) # reconstruct time series
#' mse_r = mean( (X - Xr)^2) # confirm additive reconstruction
#' plot.ts(W) # plot wavelet and scaling coefficients
#' }
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

} # EOF
