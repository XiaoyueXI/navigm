#############

## E-steps ##

#############

#' @importFrom stats pnorm dnorm
get_P1 <- function(Omega, Alpha, tau1,  v0, v1, c = 1) {

  logdens0 <-
    c * (
      stats::dnorm(
        Omega,
        mean = 0,
        sd = v0 / sqrt(tau1),
        log = TRUE
      ) + stats::pnorm(Alpha, log.p = TRUE, lower.tail = FALSE)
    )

  logdens1 <-
    c * (
      stats::dnorm(
        Omega,
        mean = 0,
        sd = v1 / sqrt(tau1),
        log = TRUE
      ) + stats::pnorm(Alpha, log.p = TRUE, lower.tail = TRUE)
    )

  # for numerical stability. = dens1 / (dens0 + dens1)
  max_logdens <-
    pmax(logdens0, logdens1) # pmax: componentwise maximum

  P1 <- exp(logdens1 - max_logdens) / (exp(logdens0 - max_logdens) + exp(logdens1 - max_logdens))
  diag(P1) <- 0 # not used

  P1
}


#' @importFrom stats dnorm
get_P2 <- function(beta, o, tau2, s0, s1, c = 1) {

  eps <- .Machine$double.eps

  logdens0 <-
    c * (stats::dnorm(
      beta,
      mean = 0,
      sd = s0 / sqrt(tau2),
      log = TRUE
    ) + log(1 - o + eps))
  logdens1 <-
    c * (stats::dnorm(
      beta,
      mean = 0,
      sd = s1 / sqrt(tau2),
      log = TRUE
    ) + log(o + eps))

  # for numerical stability. = dens1 / (dens0 + dens1)
  max_logdens <-
    pmax(logdens0, logdens1) # pmax: componentwise maximum

  exp(logdens1 - max_logdens) / (exp(logdens0 - max_logdens) + exp(logdens1 - max_logdens))

}


get_E1 <- function(P1, v0, v1) {

  (1 - P1) / v0 ^ 2 + P1 / v1 ^ 2

}

get_E2 <- function(P1, Alpha, c = 1) {

  sqrt_c <- sqrt(c)
  log_pnorm <- pnorm(sqrt_c * Alpha, log.p = TRUE)
  log_1_pnorm <-
    pnorm(sqrt_c * Alpha, log.p = TRUE, lower.tail = FALSE)

  imr0 <-
    inv_mills_ratio_(0, sqrt_c * Alpha, log_1_pnorm, log_pnorm)
  imr1 <-
    inv_mills_ratio_(1, sqrt_c * Alpha, log_1_pnorm, log_pnorm)

  E2 <- Alpha  + imr0 / sqrt_c + P1 * (imr1 - imr0) / sqrt_c

  diag(E2) <- 0

  return(E2)

}


get_E2_2 <- function(E2, Alpha, cst = 1) {

  Alpha * E2 + 1 / cst

}

get_E5 <- function(P2, s0, s1) {

  (1 - P2) / s0 ^ 2 + P2 / s1 ^ 2

}

#############

## M-steps ##

#############

get_o <- function(P2, a_o, b_o) {

  (sum(P2) + a_o - 1) / (length(P2) + a_o + b_o - 2)

}


get_tau1 <- function(Omega, E1, a_tau, b_tau, P) {

  bool_up <- upper.tri(Omega)

  (2 * a_tau - 2 + P * (P - 1) / 2) / (sum(Omega[bool_up] ^ 2 * E1[bool_up]) + 2 * b_tau)

}


get_tau2 <- function(beta, E5, a_sigma, b_sigma, Q) {

  (2 * a_sigma - 2 + Q) / (sum(beta ^ 2 * E5) + 2 * b_sigma)

}


get_zeta <- function(E2, theta, n0, t02, P) {

  (2 * n0 + 2 * t02 * sum(E2[upper.tri(E2)]) - 2 * t02 * (P - 1) * sum(theta)) / (t02 * P * (P - 1) + 2)

}


get_beta <- function(beta,
                     zeta,
                     tau2,
                     E2,
                     E5,
                     V,
                     sV,
                     pV,
                     spV1,
                     spV2,
                     P,
                     Q) {


  (sapply(sV, function(x)
    sum(x[upper.tri(x, diag = T)] * E2[upper.tri(E2)])) -
     (P - 1) *  zeta * apply(V, 2, sum) -
     (P - 1) *  sapply(1:Q, function(x)
       sum(spV1[[x]] * beta[-x])) -
     sapply(1:Q, function(x)
       sum(spV2[[x]] * beta[-x]))) /
    ((P - 1) *  apply(V ^ 2, 2, sum) + 2 * sapply(pV, function(x)
      sum(x[upper.tri(x, diag = T)])) + tau2 * E5)

}



####################

## M-steps in GMN ##

####################

get_tau2_gmn <- function(beta, a_sigma, b_sigma, Q) {

  (2 * a_sigma - 2 + Q) / (sum(beta ^ 2) + 2 * b_sigma)

}


get_beta_gmn <- function(beta,
                         zeta,
                         tau2,
                         E2,
                         V,
                         sV,
                         pV,
                         spV1,
                         spV2,
                         P,
                         Q) {


  (sapply(sV, function(x)
    sum(x[upper.tri(x, diag = T)] * E2[upper.tri(E2)])) -
     (P - 1) *  zeta * apply(V, 2, sum) -
     (P - 1) *  sapply(1:Q, function(x)
       sum(spV1[[x]] * beta[-x])) -
     sapply(1:Q, function(x)
       sum(spV2[[x]] * beta[-x]))) /
    ((P - 1) *  apply(V ^ 2, 2, sum) + 2 * sapply(pV, function(x)
      sum(x[upper.tri(x, diag = T)])) + tau2)

}

###################

## M-steps in GM ##

###################
get_rho <- function(P1, a_rho, b_rho) {

  bool_up <- upper.tri(P1)
  (sum(P1[bool_up]) + a_rho - 1) / (sum(bool_up) + a_rho + b_rho - 2)

}


