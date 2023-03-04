#' A Trous Discrete Wavelet Transform
#'
#' This function calculates the wavelet and scaling coefficients of the a trous (AT)
#' version of the Discrete Wavelet Transform (DWT).
#'
#' @param X An (N x 1) matrix or a vector
#' @param filter Filter name (see Details below) (string)
#' @param J Decomposition level (integer, 1 < J < N/2)
#' @param remove_boundary_coefs Remove boundary affected coefficients?
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
#' @examples
#' N <- 1000 #  number of time series points
#' J <- 4 # decomposition level
#' X <- matrix(rnorm(N),N,1)
#' W <- atrous_dwt(X,'coif1',J)
#' plot.ts(W) # plot wavelet and scaling coefficients
#' @export
atrous_dwt <- function(X,filter,J,remove_boundary_coefs=FALSE){

  X <- shape_check(X)

  N = length(X)
  # g = scaling_filter(filter)
  # L = length(g)
  # J = decomp_level

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

    v[,j] = scaling_coefs(x,filter,j) # modified version of Eq. 5 in Maheswaran and Khosa [2012]

    # calculate wavelet coefficients

    w[,j] = x - v[,j] # see Eq. 6 in Maheswaran and Khosa [2012]

    # calculate residual information

    x = v[,j]

  }

  # combine wavelet and scaling coefficients in a single matrix and only store the
  # last N observations

  W = cbind(w[seq(N+1,2*N),], v[seq(N+1,2*N),J])

  if (remove_boundary_coefs) {
    L <- length(scaling_filter(filter)) # could have used wavelet_filter
    num_boundary_coefs <- ((2^J)-1)*(L-1) # length of longest equivalent wavelet/scaling filter for minus 1
    W <- utils::tail(W, -num_boundary_coefs)
  }

  # atrous_dwt <- W
  return(W)

}
