AIC_GSS <- function(estimates, N){

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    2 * sum(estimates$m_delta > 0.5)

}


BIC_GSS <- function(estimates, N){

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    log(N) * sum(estimates$m_delta > 0.5)

}

EBIC_GSS <- function(estimates, gamma =0.5, N, P){

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    (log(N) + 4 * gamma * log(P)) * sum(estimates$m_delta > 0.5)

}
