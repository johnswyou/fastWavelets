#' @title MODWT Wavelet Coefficients
#' @description Compute MODWT wavelet coefficients for a specified filter and
#' decomposition level. This function only gives the MODWT wavelet coefficients
#' for the user specified level. To get all MODWT wavelet coefficients up to and
#' including the user specified decomposition level (as well as the MODWT scaling
#' coefficient corresponding to the user specified decomposition level) all in one nice
#' matrix, please use `mo_dwt()`.
#' @param X A vector
#' @param filter Filter name (string)
#' @param j Decomposition level (integer)
#' @return MODWT wavelet coefficients for specified decomposition level (vector)
#' @details This function implements the first equation in equation 169a (page 169)
#' of Wavelet Methods for Time Series Analysis by Percival and Walden.
#' @rdname modwt_wavelet_coefs
#' @export
modwt_wavelet_coefs <- function(X, filter, j) {
  equiv_filters <- equivalent_filters(filter, j, modwt=TRUE)
  equiv_wavelet_filter <- equiv_filters$equivalent_wavelet_filter
  return(circular_convolution(equiv_wavelet_filter, X))
}

#' @title MODWT Scaling Coefficients
#' @description Compute MODWT scaling coefficients for a specified filter and
#' decomposition level. This function only gives the MODWT scaling coefficients
#' for the user specified level.
#' @param X A vector
#' @param filter Filter name (string)
#' @param j Decomposition level (integer)
#' @return MODWT scaling coefficients for specified decomposition level (vector)
#' @details This function implements the second equation in equation 169a (page 169)
#' of Wavelet Methods for Time Series Analysis by Percival and Walden.
#' @rdname modwt_wavelet_coefs
#' @export
modwt_scaling_coefs <- function(X, filter, j) {
  equiv_filters <- equivalent_filters(filter, j, modwt=TRUE)
  equiv_scaling_filter <- equiv_filters$equivalent_scaling_filter
  return(circular_convolution(equiv_scaling_filter, X))
}
