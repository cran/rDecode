#' @title Implement \code{DECODE} for simple precision matrix estimation
#'
#' @description Implement \code{DECODE} to estimate a precision matrix of \code{X}. This implementation is used in Hadimaja and Pun (2018).
#'
#' @param X \eqn{n\times p}{nxp} data matrix.
#' @param lambda0 number between 0 and 1. If \code{NULL}, will use \eqn{\sqrt{2 \log{p}/n}}{\sqrt2logp/n}.
#' @param ... additional arguments to be passed to general decode function.
#'
#' @return An object of class \code{decodePM} containing:
#'   \item{Omega}{\code{DECODE} of \eqn{\Omega}.}
#'   \item{lambda0}{the \code{lambda0} used.}
#'   \item{X}{data used.}
#'   \item{theta}{final \eqn{\theta} for each column.}
#'   \item{lambda}{final \eqn{\lambda} for each column.}
#'   \item{total.iter}{number of iterations until convergence for each column.}
#'
#' @examples
#' # estimate the precision matrix of iris data
#' object <- decodePM(iris[,1:4], lambda0 = 0.01)
#'
#' object
#' summary(object)
#'
#' object$Omega
#'
#' @references
#' Hadimaja, M. Z., & Pun, C. S. (2018). A Self-Calibrated Regularized Direct Estimation for Graphical Selection and Discriminant Analysis.
#'
#' @export

decodePM <- function(X, lambda0 = NULL, ...){

  N <- dim(X)[1]
  P <- dim(X)[2]

  # sample covariance matrix
  Sigma <- stats::cov(X)

  # lambda0
  if (is.null(lambda0)) lambda0 <- sqrt(2 * log(P)/N)

  Omega <- matrix(0, P, P)
  theta <- numeric(P)
  lambda <- numeric(P)
  total.iter <- numeric(P)

  # estimate each column
  for (p in 1:P) {
    e <- numeric(P)
    e[p] <- 1
    decode.object <- decode(Sigma, e, lambda0, ...)
    Omega[,p] <- decode.object$eta
    theta[p] <- decode.object$theta
    lambda[p] <- decode.object$lambda
    total.iter[p] <- decode.object$total.iter
  }

  # symmetrization step
  Omega <- Omega * (abs(Omega) <= abs(t(Omega))) + t(Omega) * (abs(Omega) > abs(t(Omega)))

  # naming
  vars <- colnames(X)
  dimnames(Omega) <- list(vars, vars)
  names(theta) <- vars
  names(lambda) <- vars
  names(total.iter) <- vars

  object <- list()
  object$Omega = Omega
  object$theta = theta
  object$lambda = lambda
  object$lambda0 = lambda0
  object$method = decode.object$method
  object$X = X
  object$total.iter = total.iter
  object$call <- match.call()

  class(object) <- 'decodePM'
  return(object)
}

#' @export
print.decodePM <- function (x, ...)
{
  message("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nOmega:\n")
  print.default(format(x$Omega, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  invisible(x)
}

#' @export
summary.decodePM <- function (object, ...)
{
  message("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nTotal iter:\n")
  print.default(format(object$total.iter, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  message("\nTheta:\n")
  print.default(format(object$theta, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  message("\nLambda:\n")
  print.default(format(object$lambda, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  message("\nOmega:\n")
  print.default(format(object$Omega, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  invisible(object)
}
