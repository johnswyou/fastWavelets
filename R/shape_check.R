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
