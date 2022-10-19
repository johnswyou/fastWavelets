#' Maximal Overlap Discrete Wavelet Transform (MODWT)
#'
#' This function calculates the wavelet and scaling coefficients of the
#' MODWT.
#'
#' @param X Time series (N x 1)
#' @param wavelet Scaling filter name (string)
#' @param decomp_level Decomposition level (1 < integer < N/2)
#' @return Wavelet and scaling coefficients (N x J+1) (wavelet coefs in first J columns)
#' @references
#'
#' M. Basta (2014), Additive Decomposition and Boundary Conditions in Wavelet-Based
#' Forecasting Approaches, Acta Oeconomica Pragensia, 2, pp. 48-70.
#'
#' @examples
#' \dontrun{
#' N = 1000 #  number of time series points
#' J = 4 # decomposition level
#' wavelet = 'coif1' # wavelet filter
#' X = matrix(rnorm(N),N,1)
#' W = mo_dwt(X,wavelet,J)
#' }
#' @export
mo_dwt <- function(X,wavelet,decomp_level){

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

    v[,j] = scaling_coefs(x,wavelet,j) # Eq. 14 in Basta [2014]

    # calculate wavelet coefficients

    w[,j] = wavelet_coefs(x,wavelet,j) # Eq. 13 in Basta [2014]

    # scaling coefficients of previous step are used for calculating scaling
    # and wavelet in next step

    x = v[,j]

  }

  # combine wavelet and scaling coefficients in a single matrix and only store the
  # last N observations

  W = cbind(w[seq(N+1,2*N),], v[seq(N+1,2*N),J])

  # mo_dwt <- W
  return(W)

} # EOF
