#' @title Descent-based Calibrated Optimal Direct Estimation
#'
#' @description Implement \code{DECODE} for \code{sigma} and \code{beta} to estimate \eqn{\Sigma^{-1}\beta}{\Sigma^-1\beta} where \code{sigma} is an estimator of \eqn{\Sigma} and \code{beta} is an estimator of \eqn{\beta}.
#'
#' @param sigma \eqn{p \times p}{pxp} positive semidefinite symmetric matrix. \code{sigma} will be perturbed if needed.
#' @param beta \eqn{p}-length vector.
#' @param lambda0 number between 0 and 1.
#' @param decode.tol error tolerance for \code{DECODE}.
#' @param decode.maxit maximum iterations for \code{DECODE}
#' @param trace logical. If \code{TRUE}, will return \eqn{\eta}, \eqn{\theta}, and \eqn{\lambda} found during each iteration of \code{DECODE}
#' @param solver solver for \eqn{\ell_1}{l1}-RQP problem inside \code{DECODE}.
#' @param solver.tol tolerance for solver.
#' @param solver.maxit maximum iterations for solver (only for APG).
#' @param return.sigma logical. If \code{TRUE} the \code{sigma} entered is returned.
#' @param return.beta logical. If \code{TRUE} the \code{beta} entered is returned.
#' @param return.param logical. If \code{TRUE} the parameters used are returned.
#'
#' @return An object of class \code{decode} containing:
#'   \item{eta}{\code{DECODE} of \eqn{\Sigma^{-1}\beta}{\Sigma^-1\beta}.}
#'   \item{theta}{final \eqn{\theta} of the \code{DECODE}.}
#'   \item{lambda}{final \eqn{\lambda} of the \code{DECODE}.}
#'   \item{sigma.mult}{multiplier applied on \code{sigma} to ensure convergence.}
#'   \item{total.iter}{number of iterations until convergence.}
#'   \item{call}{the matched call.}
#'   \item{method}{the solver used, if requested.}
#'   \item{lambda0}{the \code{lambda0} entered, if requested.}
#'   \item{decode.tol}{the \code{decode.tol} used, if requested.}
#'   \item{decode.maxit}{the \code{decode.maxit} used, if requested.}
#'   \item{trace}{the \code{trace} used, if requested.}
#'   \item{solver.tol}{the \code{solver.tol} used, if requested.}
#'   \item{solver.maxit}{the \code{solver.maxit} used, if requested.}
#'   \item{eta.trace}{matrix of \eqn{\eta} used in each iteration, if requested.}
#'   \item{theta.trace}{vector of \eqn{\theta} used in each iteration, if requested.}
#'   \item{lambda.trace}{vector of \eqn{\lambda} used in each iteration, if requested.}
#'
#' @examples
#' # estimate A^(-1) b with a certain lambda0
#' X <- matrix(rnorm(100), 10, 10)
#' A <- t(X) %*% X
#' b <- rnorm(10)
#' object <- decode(A, b, lambda0 = 0.8)
#'
#' object
#' summary(object)
#'
#' coef(object)
#'
#' @references
#' Pun, C. S. (2018). A Sparse Learning Approach to Relative-Volatility-Managed Portfolio Selection.
#' Hadimaja, M. Z., & Pun, C. S. (2018). A Self-Calibrated Regularized Direct Estimation for Graphical Selection and Discriminant Analysis.
#' @export

decode <- function(sigma,
                   beta,
                   lambda0,
                   decode.tol = 1e-06,
                   decode.maxit = 100,
                   trace = FALSE,
                   solver = c("apg", "homotopy"),
                   solver.tol = 1e-08,
                   solver.maxit = 10000,
                   return.sigma = FALSE,
                   return.beta = FALSE,
                   return.param = FALSE) {

  call <- match.call()
  solver <- match.arg(solver)

  # check lambda0
  if (!is.numeric(lambda0) | length(lambda0)!=1) stop("lambda0 has to be a number")
  if (lambda0 < 0) stop("Lambda0 has to be strictly larger than 0")

  # check sigma and beta
  if (!is.numeric(sigma)) stop("sigma has to be numeric")
  if (!is.numeric(beta)) stop("beta has to be numeric")
  if (!is.matrix(sigma)) stop("sigma has to be a matrix")
  if (!is.vector(beta)) stop("beta has to be a vector")
  if (dim(sigma)[1] != dim(sigma)[2]) stop("sigma has to be a square matrix")
  if (dim(sigma)[1] != length(beta)) stop("sigma and beta have incompatible dimensions")
  if (!isSymmetric(sigma)) stop("sigma and beta have incompatible dimensions")

  # store original sigma and beta if requested
  if(return.sigma) sigma.ori <- sigma
  if(return.beta) beta.ori <- beta

  # perturb sigma if needed
  psd <-  all(eigen(sigma, TRUE, TRUE)$values > 1e-4)
  if (!psd) sigma <- sigma + lambda0 * diag(diag(sigma))

  # adjust sigma and beta
  diag.min.half <- diag(sqrt(diag(sigma)) ^ (-1))
  sigma <- diag.min.half %*% sigma %*% diag.min.half
  beta <- drop(diag.min.half %*% beta)

  # compute multiplier for sigma
  lambda.max = max(abs(beta))
  sigma.mult = (as.numeric(crossprod(beta, solve(sigma, beta))) / lambda.max) ^ 2
  if(sigma.mult * lambda0 ^ 2 < 1) {
    sigma.mult <- 1
  }
  sigma <- sigma.mult * sigma

  # run decodeIter
  decode.result <-
    decodeIter(
      sigma,
      beta,
      lambda0,
      decode.tol = decode.tol/sqrt(sigma.mult),
      decode.maxit = decode.maxit,
      trace = trace,
      solver = solver,
      solver.tol = solver.tol,
      solver.maxit = solver.maxit
    )

  decode.result$eta <- as.vector(diag.min.half %*% decode.result$eta) * sigma.mult
  decode.result$theta <- decode.result$theta * sqrt(sigma.mult)
  decode.result$lambda <- decode.result$lambda * sqrt(sigma.mult)

  if(trace) {
    decode.result$eta.trace <- decode.result$eta.trace %*% diag.min.half * sigma.mult
    decode.result$theta.trace <- decode.result$theta.trace * sqrt(sigma.mult)
    decode.result$lambda.trace <- decode.result$lambda.trace * sqrt(sigma.mult)
  }

  object <- list(call = call,
                 method = ifelse(solver == 'apg', 'APG', 'Homotopy'),
                 sigma.mult = sigma.mult)
  if(return.sigma) object$sigma <- sigma.ori
  if(return.beta) object$beta <- beta.ori
  if(return.param) object <- c(object, list(lambda0 = lambda0,
                                            decode.tol = decode.tol,
                                            decode.maxit = decode.maxit,
                                            trace = trace,
                                            solver.tol = solver.tol,
                                            solver.maxit = solver.maxit))
  object <- c(object, decode.result)


  class(object) <- 'decode'
  object
}

decodeIter <- function(sigma,
                       beta,
                       lambda0,
                       decode.tol = 1e-06,
                       decode.maxit = 100,
                       trace = FALSE,
                       solver = c("apg", "homotopy"),
                       solver.tol = 1e-08,
                       solver.maxit = 10000) {


  solver <- match.arg(solver)

  theta.old <- 0
  flag <- 0

  use.apg <- (solver == 'apg')

  if(use.apg) {
    eta = RqpAPG(sigma, beta, lambda=max(abs(beta))/2, tol=solver.tol, maxit=solver.maxit)
  } else {
    rqp.path = RqpHmtp(sigma, beta, tol=solver.tol)
    lambda.path = rqp.path$lambda.path
    eta.path = rqp.path$eta.path
    eta = eta.path[floor(dim(eta.path)[1]) / 2, ]
  }
  theta.new <- Theta(eta,sigma,beta)
  lambda <- theta.new * lambda0

  # trace variables
  if (trace) {
    P = length(beta)
    eta.trace <- matrix(0, decode.maxit, P)
    lambda.trace <- double(decode.maxit)
    theta.trace <- double(decode.maxit)
  }

  # DECODE
  while ((abs(theta.new - theta.old) > decode.tol) & (flag < decode.maxit)) {

    flag <- flag + 1

    theta.old <- theta.new

    # 1
    lambda <- theta.new * lambda0

    # 2
    if(use.apg) {
      eta <- RqpAPG(sigma, beta, lambda, tol=solver.tol, maxit=solver.maxit)
    } else {
      eta <- HmtpSolve(lambda, lambda.path, eta.path)
    }

    # store trace variables
    if (trace) {
      eta.trace[flag, ] <- eta
      theta.trace[flag] <- theta.new
      lambda.trace[flag] <- lambda
    }
    # 3
    theta.new <- Theta(eta,sigma,beta)

  }
  object <- list(
    eta = eta,
    theta = theta.new,
    lambda = lambda,
    total.iter = flag
  )
  if(trace) {
    object$eta.trace = eta.trace[1:flag,]
    object$theta.trace = theta.trace[1:flag]
    object$lambda.trace = lambda.trace[1:flag]
  }
  object
}

L = function(eta, sigma, beta) {
  drop(0.5 * t(eta) %*% sigma %*% eta - t(beta) %*% eta)
}

Theta <- function(eta, sigma, beta) {
  sqrt(-2 * L(eta, sigma, beta))
}

#' @export
print.decode <- function (x, ...)
{
  message("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nEta:\n")
  print.default(format(x$eta,...), print.gap = 2L, quote = FALSE)
  message("\n")
  invisible(x)
}


#' @export
summary.decode <- function (object, ...)
{
  message("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nTotal iter:\n", paste(object$total.iter, sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nTheta:\n", paste(format(object$theta, ...), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nLambda:\n", paste(format(object$lambda, ...), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nSigma multiplier:\n", paste(format(object$sigma.mult, ...), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nEta:\n")
  print.default(format(object$eta, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  invisible(object)
}

#' @export
coef.decode <- function (object, ...) object$eta
