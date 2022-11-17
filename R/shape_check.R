#' @title Check Matrix X is (N x 1)
#' @description `shape_check` checks whether `X` is a matrix representing
#' a column vector (i.e., a matrix with 1 column). If not, `shape_check` attempts
#' to coerce the user provided `X` to a matrix with 1 column. If this cannot be done,
#' an error is raised.
#' @param X Object to check and (if possible) coerce to a single column matrix
#' @return An (N x 1) matrix
#' @details This is a utility function written to check the input `X` for the
#' functions `atrous_dwt` and `mo_dwt`.
#' @export
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
