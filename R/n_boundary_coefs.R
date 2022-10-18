#' Number of Boundary Coefficients
#'
#' This function calculates the number of boundary coefficients for a
#' particular wavelet/scaling filter and decomposition level.
#'
#' @param wavelet scaling filter name \[string\]
#' @param decomp_level decomposition level \[integer\]
#' @return number of boundary coefficients \[integer\]
#'
#' @references
#'
#' M. Basta (2014),Additive Decomposition and Boundary Conditions in Wavelet-Based
#' Forecasting Approaches, Acta Oeconomica Pragensia, 2, pp. 48-70.
#'
#' @examples
#' \dontrun{
#' J = 4 # decomposition level
#' wavelet = 'b3spline' # wavelet filter
#' nbc = n_boundary_coefs(wavelet,J) # number of boundary-effected coefficients at decomp_level J
#' }
#' @export
n_boundary_coefs <- function(wavelet, decomp_level){

  g = scaling_filter(wavelet)
  L = length(g)
  nbc = ((2^decomp_level) - 1) * (L - 1) + 1; # See Eq. 10 in Basta [2014]

  # n_boundary_coefs <- nbc
  return(nbc)

} # EOF
