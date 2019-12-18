#' @title Implement \code{DECODE} for simple LDA
#'
#' @description Implement \code{DECODE} for simple LDA. The LDA assumes both classes have equal prior probabilities. This implementation is used in Hadimaja and Pun (2018).
#'
#' @param X \eqn{n\times p}{nxp} data matrix.
#' @param y binary \eqn{n}-length vector containing the class of each observation.
#' @param lambda0 number between 0 and 1. If \code{NULL}, will use \eqn{\sqrt{2 \log{p}/n}}{\sqrt2logp/n}.
#' @param ... additional arguments to be passed to general decode function.
#'
#' @return An object of class \code{decodeLDA} containing:
#'   \item{eta}{\code{DECODE} of \eqn{\Omega\delta}}
#'   \item{X}{training data used}
#'   \item{y}{training label used}
#' and various outputs from \code{decode} function.
#'
#' @examples
#' # for efficiency, we will only use 500 variables
#'
#' # load the training data (Lung cancer data, cleaned)
#' data(lung.train) # 145 x 1578
#' X.train <- lung.train[,1:500]
#' y.train <- lung.train[,1578]
#'
#' # build the DECODE
#' object <- decodeLDA(X.train, y.train)
#'
#' object
#' summary(object)
#' coef(object)
#'
#' # test on test data
#' data(lung.test)
#' X.test <- lung.test[,1:500]
#' y.test <- lung.test[,1578]
#' y.pred <- predict(object, X.test)
#' table(y.pred, y.test)
#'
#' @references
#' Hadimaja, M. Z., & Pun, C. S. (2018). A Self-Calibrated Regularized Direct Estimation for Graphical Selection and Discriminant Analysis.
#'
#' @export

decodeLDA <- function(X, y, lambda0 = NULL, ...){

  N <- dim(X)[1]
  P <- dim(X)[2]
  if (length(y) != N) stop("size mismatch")

  y <- as.factor(y)
  c <- levels(y)
  if(length(c) != 2) stop("too many classes")

  # X1
  X1 <- X[y==c[1],]
  N1 <- dim(X1)[1]
  S1 <- stats::cov(X1)
  mu1 <- colMeans(X1)

  # X2
  X2 <- X[y==c[2],]
  N2 <- dim(X2)[1]
  S2 <- stats::cov(X2)
  mu2 <- colMeans(X2)

  # combined covariance matrix
  Sigma <- ((N1-1)*S1 + (N2-1)*S2)/(N-2)

  # delta
  delta <- mu1 - mu2

  # lambda0
  if (is.null(lambda0)) lambda0 <- sqrt(2 * log(P)/N)


  decode.object <- decode(Sigma, delta, lambda0, ...)

  vars <- colnames(X)

  object <- c(decode.object, X = list(X), y = list(y))
  names(object$eta) <- vars
  object$call <- match.call()
  object$lambda0 <- lambda0
  class(object) <- 'decodeLDA'

  return(object)
}

#' @export
print.decodeLDA <- function (x, ...)
{
  message("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nEta:\n")
  print.default(format(x$eta, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  invisible(x)
}

#' @export
summary.decodeLDA <- function (object, ...)
{
  message("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nTotal iter:\n", paste(object$total.iter, sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nTheta:\n", paste(format(object$theta, ...), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nLambda:\n", paste(format(object$lambda, ...), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nSigma multiplier:\n", paste(object$sigma.mult, sep = "\n", collapse = "\n"),
      "\n", sep = "")
  message("\nEta:\n")
  print.default(format(object$eta, ...), print.gap = 2L, quote = FALSE)
  message("\n")
  invisible(object)
}

#' @export
coef.decodeLDA <- function (object, ...) object$eta

#' @export
predict.decodeLDA <- function(object, newdata, ...)
{
  if(is.null(newdata)) newdata = object$X
  means <- colMeans(object$X)
  c <- levels(object$y)
  newdata <- scale(newdata, means, FALSE)
  y.pred <- drop(newdata %*% object$eta) >= 0
  y.pred <- as.factor(ifelse(y.pred, c[1], c[2]))
  y.pred
}
