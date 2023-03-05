#' @title Scaling Filter
#' @description This function returns the user specified scaling filter.
#' @param filter Name of the scaling filter desired \[string\]
#' @param modwt Should the elements of the returned filter be scaled by `1/sqrt(2)`?, Default: `FALSE`
#' @return A scaling filter vector
#' @details See equation 75a (page 75) of Wavelet Methods for Time Series Analysis
#' by Percival and Walden to see how to convert a wavelet filter into a scaling filter.
#'
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
#' \dontrun{
#' if(interactive()){
#'  r_scaling_filter("db1")
#'  }
#' }
#' @rdname r_scaling_filter
#' @export
r_scaling_filter <- function(filter, modwt=FALSE) {

  supported_filters <- c('bl7', 'bl9', 'bl10', 'beyl', 'coif1', 'coif2',
                             'coif3', 'coif4', 'coif5', 'db1', 'db2', 'db3',
                             'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10',
                             'db11', 'db12', 'db13', 'db14', 'db15', 'db16',
                             'db17', 'db18', 'db19', 'db20', 'db21', 'db22',
                             'db23', 'db24', 'db25', 'db26', 'db27', 'db28',
                             'db29', 'db30', 'db31', 'db32', 'db33', 'db34',
                             'db35', 'db36', 'db37', 'db38', 'db39', 'db40',
                             'db41', 'db42', 'db43', 'db44', 'db45', 'fk4',
                             'fk6', 'fk8', 'fk14', 'fk18', 'fk22', 'han2_3',
                             'han3_3', 'han4_5', 'han5_5', 'dmey', 'mb4_2',
                             'mb8_2', 'mb8_3', 'mb8_4', 'mb10_3', 'mb12_3',
                             'mb14_3', 'mb16_3', 'mb18_3', 'mb24_3', 'mb32_3',
                             'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7',
                             'sym8', 'sym9', 'sym10', 'sym11', 'sym12', 'sym13',
                             'sym14', 'sym15', 'sym16', 'sym17', 'sym18',
                             'sym19', 'sym20', 'sym21', 'sym22', 'sym23',
                             'sym24', 'sym25', 'sym26', 'sym27', 'sym28',
                             'sym29', 'sym30', 'sym31', 'sym32', 'sym33',
                             'sym34', 'sym35', 'sym36', 'sym37', 'sym38',
                             'sym39', 'sym40', 'sym41', 'sym42', 'sym43',
                             'sym44', 'sym45', 'vaid', 'la8', 'la10', 'la12',
                             'la14', 'la16', 'la18', 'la20')

  if (!(filter %in% supported_filters)) {
    stop('The filter entered is not currently supported')
  }

  out <- scaling_list[[filter]]

  if (modwt) {
    out <- out / sqrt(2)
  }

  out <- as.vector(out)

  return(out)
}

#' @title Wavelet Filter
#' @description This function returns the user specified wavelet filter.
#' @param filter Name of the wavelet filter desired \[string\]
#' @param modwt Should the elements of the returned filter be scaled by `1/sqrt(2)`?, Default: `FALSE`
#' @return A wavelet filter vector
#' @details See section 4.2 (page 68) of Wavelet Methods for Time Series Analysis
#' by Percival and Walden for a detailed discussion of the wavelet filter.
#'
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
#' \dontrun{
#' if(interactive()){
#'  r_wavelet_filter("coif1")
#'  }
#' }
#' @rdname r_wavelet_filter
#' @export
r_wavelet_filter <- function(filter, modwt=FALSE) {

  supported_filters <- c('bl7', 'bl9', 'bl10', 'beyl', 'coif1', 'coif2',
                             'coif3', 'coif4', 'coif5', 'db1', 'db2', 'db3',
                             'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10',
                             'db11', 'db12', 'db13', 'db14', 'db15', 'db16',
                             'db17', 'db18', 'db19', 'db20', 'db21', 'db22',
                             'db23', 'db24', 'db25', 'db26', 'db27', 'db28',
                             'db29', 'db30', 'db31', 'db32', 'db33', 'db34',
                             'db35', 'db36', 'db37', 'db38', 'db39', 'db40',
                             'db41', 'db42', 'db43', 'db44', 'db45', 'fk4',
                             'fk6', 'fk8', 'fk14', 'fk18', 'fk22', 'han2_3',
                             'han3_3', 'han4_5', 'han5_5', 'dmey', 'mb4_2',
                             'mb8_2', 'mb8_3', 'mb8_4', 'mb10_3', 'mb12_3',
                             'mb14_3', 'mb16_3', 'mb18_3', 'mb24_3', 'mb32_3',
                             'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7',
                             'sym8', 'sym9', 'sym10', 'sym11', 'sym12', 'sym13',
                             'sym14', 'sym15', 'sym16', 'sym17', 'sym18',
                             'sym19', 'sym20', 'sym21', 'sym22', 'sym23',
                             'sym24', 'sym25', 'sym26', 'sym27', 'sym28',
                             'sym29', 'sym30', 'sym31', 'sym32', 'sym33',
                             'sym34', 'sym35', 'sym36', 'sym37', 'sym38',
                             'sym39', 'sym40', 'sym41', 'sym42', 'sym43',
                             'sym44', 'sym45', 'vaid', 'la8', 'la10', 'la12',
                             'la14', 'la16', 'la18', 'la20')

  if (!(filter %in% supported_filters)) {
    stop('The filter entered is not currently supported')
  }

  out <- wavelet_list[[filter]]

  if (modwt) {
    out <- out / sqrt(2)
  }

  out <- as.vector(out)

  return(out)
}

#' @title Equivalent Wavelet and Scaling Filters
#' @details Equivalent wavelet and scaling filters are useful for obtaining
#' MODWT wavelet and scaling coefficients for a given filter name and
#' decomposition level.
#' @param filter Filter name (string)
#' @param J Decomposition level (integer)
#' @param modwt Do you want the MODWT version?, Default: FALSE
#' @return A list containing the equivalent wavelet and scaling filter and their lengths.
#' @rdname equivalent_filters
#' @export
# Reference: wavelets R package
equivalent_filters <- function(filter, J, modwt=FALSE) {

  h <- r_wavelet_filter(filter, modwt=FALSE)
  g <- r_scaling_filter(filter, modwt=FALSE)
  L <- length(h)

  h.last <- h
  g.last <- g
  L.last <- L

  for(j in 2:J){
    L.new <- (2^j - 1)*(L-1) + 1
    hj <- NULL
    gj <- NULL
    for(l in 0:(L.new - 1)){
      u <- l
      ifelse(u >= L, g.mult <- 0, g.mult <- g[u+1])
      hjl <- g.mult*h.last[1]
      gjl <- g.mult*g.last[1]
      for(k in 1:(L.last-1)){
        u <- u-2
        if((u < 0) | (u >= L)) g.mult <- 0 else g.mult <- g[u+1]
        hjl <- hjl + g.mult*h.last[k+1]
        gjl <- gjl + g.mult*g.last[k+1]
      }
      hj <- c(hj,hjl)
      gj <- c(gj,gjl)
    }
    h.last <- hj
    g.last <- gj
    L.last <- L.new
  }
  if (modwt) {
    h.last <- h.last/(2^(J/2))
    g.last <- g.last/(2^(J/2))
  }
  return(list(equivalent_wavelet_filter = h.last,
              equivalent_scaling_filter = g.last,
              Lj = L.last))
}
