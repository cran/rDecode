RqpAPG <- function(sigma, beta, lambda, tol=1e-8, maxit=10000) {

  # gradient of quadratic function 1/2 eta' Sigma eta - beta' eta
  # equal to sigma eta - beta
  GradQuad <- function(eta, sigma, beta) {
    return(sigma %*% eta - beta)
  }

  # proximal operator of l1-norm
  ProxL1 <- function(eta, t, lambda) {
    ## Compute the soft-thresholded operator
    thres <- t * lambda
    idx.1 <- which(eta < -thres)
    idx.2 <- which(eta > thres)
    prox <- rep(0, length(eta))
    if (length(idx.1) > 0) prox[idx.1] <- eta[idx.1] + thres
    if (length(idx.2) > 0) prox[idx.2] <- eta[idx.2] - thres
    return(prox)
  }

  # function wrappers for GradF and ProxH
  GradF <- function(eta) GradQuad(eta, sigma, beta)
  ProxH <- function(eta, t) ProxL1(eta, t, lambda)

  APG(GradF, ProxH, length(beta), tol, maxit)
}

APG <- function(GradF, ProxH, P, tol=1e-8, maxit=10000) {

  # Set default parameters
  x0 <- numeric(P) # initial starting point
  restart <- TRUE # use adaptive restart scheme
  growth <- 1.01 # step-size growth factor
  shrinkage <- 0.5 # step-size shrinkage factor
  quiet <- TRUE # if false writes out information every 100 iters
  step.size = NULL # starting step-size estimate, if not set then apg makes initial guess
  fixed.step.size <- FALSE # don't change step-size (forward or back tracking), uses initial step-size throughout, only useful if good STEP_SIZE se

  # Initialization
  x <- x0
  y <- x
  g <- GradF(y)
  if (sqrt(sum(g^2)) < tol) {
    return(list(x=x,t=0))
  }
  theta <- 1

  # Initial step size
  if (is.null(step.size)) {

    # Barzilai-Borwein step-size initialization:
    t <- 1 / sqrt(sum(g^2))
    x_hat <- x - t*g
    g_hat <- GradF(x_hat)
    t <- abs(sum( (x - x_hat)*(g - g_hat)) / (sum((g - g_hat)^2)))
  } else {
    t <- step.size
  }

  # Main loop
  for (k in seq(maxit)) {

    if (!quiet && (k %% 100==0)) {
      message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
    }

    x_old <- x
    y_old <- y

    # The proximal gradient step (update x)
    x <- ProxH(y - t*g, t)

    # The error for the stopping criterion
    err1 <- sqrt(sum((y-x)^2)) / max(1, sqrt(sum(x^2)))
    if (err1 < tol) break

    # Update theta for acceleration
    theta <- 2/(1 + sqrt(1+4/(theta^2)))

    # Update y
    if (restart && sum((y-x)*(x-x_old))>0) {
      x <- x_old
      y <- x
      theta <- 1
    } else {
      y <- x + (1-theta)*(x-x_old)
    }

    # New gradient
    g_old <- g
    g <- GradF(y)

    # Update step size by TFOCS-style backtracking
    if (!fixed.step.size) {
      t_hat <- 0.5*sum((y-y_old)^2)/abs(sum((y - y_old)*(g_old - g)))
      t <- min( growth*t, max( shrinkage*t, t_hat ))
    }
  }
  if (!quiet) {
    message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
    if (k==maxit) message(paste('Warning: maximum number of iterations reached'))
    message('Terminated')
  }

  # return solution
  return(x)
}
