RqpHmtp <- function(sigma, beta, tol = .Machine$double.eps, trace = FALSE) {

  P = length(beta)

  eta <- double(P)
  lam = max(abs(beta))
  res = 0.5 * beta

  J = NULL
  I = (1:P)

  k = 0

  Eta = eta
  Lam = lam
  Det = 0

  if (trace) {
    U = double(P)
    V = double(P)
    Res = res
  }


  while((length(J) < P) & (lam >= tol)) {
    k = k + 1
    j = intersect(I, which(abs(abs(res) - lam/2) <= tol))
    J = union(J, j)
    I = setdiff(I, j)
    j = NULL


    if (length(J) > 1) Det = c(Det, det(sigma[J, J])) else Det[1] = sigma[J, J]

    u = double(P)
    u[J] <- u[J] <- 2 * solve(sigma[J, J], sign(res[J]))
    v = drop(0.5 * sigma %*% u)

    if (length(I) == 0) {
      gamma = lam / 2
    } else {
      gamma = c((lam/2 - res[I]) / (1 - v[I]),
                (lam/2 + res[I]) / (1 + v[I]))
      gamma = min(gamma[gamma > tol])
    }

    eta = eta + gamma * u
    res <- res - gamma * v
    lam = lam - 2 * gamma

    Eta = rbind(Eta, eta)

    if (lam <= 1e-10) {
      lam = 0
      break
    }
    Lam = c(Lam, lam)

    if (trace) {
      U = rbind(U, u)
      V = rbind(V, v)
      Res = rbind(Res, res)
    }
  }

  Eta = unname(Eta)

  if (trace) {
    U = unname(U)
    V = unname(V)
    Res = unname(Res)
    return(list(seq = J,
                lambda.path = Lam,
                eta.path = Eta,
                u = U,
                v = V,
                res = Res,
                Det = Det))
  }

  return(list(lambda.path = Lam,
              eta.path = Eta))
}

HmtpSolve <- function(lambda, lambda.path, eta.path) {
  i = sum(lambda < lambda.path)
  if (i == 0) {
    return(double(dim(eta.path)[2]))
  } else if (i == length(lambda.path)) {
    stop("Lambda is too small")
  } else {
    return(eta.path[i,] + (eta.path[i+1,] - eta.path[i,]) * (lambda - lambda.path[i]) / (lambda.path[i+1] - lambda.path[i]))
  }
}
