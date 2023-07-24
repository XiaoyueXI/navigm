# This file is part of the `navigss` R package:
#     https://github.com/XiaoyueXI/navigss
#
# Internal functions gathering the updates for the core algorithms.

####################

## tau's updates ##

####################

update_alpha_tau <- function(a_tau, P, c = 1) {
  c * (P * (P - 1) / 4 + a_tau  - 1) + 1
}

update_beta_tau <- function(Omega, E1, b_tau, c = 1) {
  bool_up <- upper.tri(Omega)
  c * (b_tau + sum(Omega[bool_up] ^ 2 * E1[bool_up]) / 2)
}

get_m_tau <- function(alpha_tau, beta_tau) {
  alpha_tau / beta_tau
}

get_m_log_tau <- function(alpha_tau,beta_tau) {
  eps <- .Machine$double.eps
  digamma(alpha_tau + eps) - log(beta_tau + eps)
}

####################

## zeta's updates ##

####################

update_sig2_inv_zeta <- function(t02, P, c = 1) {
  c * (P * (P - 1) / 2 + 1 / t02)
}

update_mu_zeta <-
  function(m1_beta, sig2_inv_zeta, m_z, V, n0, t02, P, c = 1) {
    bool_up <- upper.tri(m_z)
    c * (sum(m_z[bool_up]) +
           n0 / t02 -
           (P - 1) * sum(apply(V, 2, sum) * m1_beta)) / sig2_inv_zeta
  }

get_m2_zeta <- function(mu_zeta, sig2_inv_zeta) {
  mu_zeta ^ 2 + sig2_inv_zeta ^ (-1)
}

update_mu_zeta_gm <-
  function(sig2_inv_zeta, m_z, n0, t02, c = 1) {
    bool_up <- upper.tri(m_z)
    c * (sum(m_z[bool_up]) + n0 / t02) / sig2_inv_zeta
  }

####################

## sigma's updates ##

####################

update_alpha_sigma <- function(m_gamma, a_sigma, c = 1) {
  c * (sum(m_gamma) / 2 + a_sigma -1) + 1
}

update_beta_sigma <- function(m_gamma, m2_beta, b, c = 1) {
  c * (sum(m_gamma * m2_beta) / 2 + b_sigma)
}

get_m_sig2_inv <- function(alpha_sigma, beta_sigma){
  alpha_sigma / beta_sigma
}

get_m_log_sig2_inv <- function(alpha_sigma, beta_sigma) {
  eps <- .Machine$double.eps
  digamma(alpha_sigma + eps) - log(beta_sigma + eps)
}

####################

## o's updates ##

####################

update_alpha_o <- function(m_gamma, a_o, c = 1) {
  # c * (m_gamma + a_o - 1) + 1
  c * (effective_sum(m_gamma) + a_o - 1) + 1
}

update_beta_o <- function(m_gamma, b_o, c = 1) {
  c * (effective_sum(1 - m_gamma) + b_o - 1) + 1
}

get_m_o <- function(alpha_o, beta_o) {
  alpha_o / (alpha_o + beta_o)
}

get_m_log_o <- function(alpha_o, beta_o) {
  eps <- .Machine$double.eps
  digamma(alpha_o + eps) - digamma(alpha_o + beta_o + eps)
}

get_m_log_one_minus_o <- function(alpha_o, beta_o) {
  eps <- .Machine$double.eps
  digamma(beta_o + eps) - digamma(alpha_o + beta_o + eps)
}

####################

## gamma's updates ##

####################

update_m_gamma <-
  function(m_log_o,
           m_log_one_minus_o,
           m_log_sig2_inv,
           mu_beta,
           sig2_inv_beta,
           c = 1) {

    eps <- .Machine$double.eps
    (1 + exp(
      c * (m_log_one_minus_o -  m_log_o -
             m_log_sig2_inv/2) +
        log(sig2_inv_beta + eps)/2 -
        mu_beta ^ 2 * sig2_inv_beta/2
    )) ^ (-1)
  }

####################

## beta's updates ##

####################

update_sig2_inv_beta <- function(m_sig2_inv, V, pV, p, c = 1) {
  c *  (m_sig2_inv + (p - 1) * apply(V ^ 2, 2, sum) + 2 * sapply(pV, function(x)
    sum(x[upper.tri(x, diag = T)])))
}


update_mu_beta <-
  function(sig2_inv_beta,
           m1_beta,
           mu_zeta,
           E2,
           V,
           sV,
           spV1,
           spV2,
           p,
           q,
           c = 1) {
    bool_up <- upper.tri(E2)
    c * (
      sapply(sV, function(x) {
        sum(x[upper.tri(x, diag = T)] * E2[bool_up])
      }) -
        (p - 1) * mu_zeta * apply(V, 2, sum) -
        (p - 1) * sapply(1:q, function(j) {
          sum(spV1[[j]] * m1_beta[-j])
        }) -
        sapply(1:q, function(j) {
          sum(spV2[[j]] * m1_beta[-j])
        })
    ) / sig2_inv_beta
  }

get_m1_beta <- function(mu_beta, m_gamma) {
  mu_beta * m_gamma
}

get_m2_beta <- function(mu_beta, sig2_inv_beta, m_gamma) {
  (mu_beta ^ 2 + sig2_inv_beta ^ (-1)) * m_gamma^2
}


####################

## delta's updates ##

####################

update_m_delta <- function(Omega, m_tau, m1_alpha, v0, v1, c=1) {
  1 / (1 + exp(
    c * log(v1 / v0) +
      c * m_tau * Omega ^ 2 * (1 / v1 ^ 2 - 1 / v0 ^ 2) / 2 +
      pnorm(sqrt(c) * m1_alpha, log.p = T, lower.tail = F) -
      pnorm(sqrt(c) * m1_alpha, log.p = T, lower.tail = T)
  ))
}

get_sum_var_alpha <- function(sig2_inv_zeta, m1_beta, m2_beta, V, pV, p){
  p* (p-1)/2 * sig2_inv_zeta^(-1) +
    (p-1) * sum(apply(V^2, 2, sum) * (m2_beta - m1_beta^2)) +
    2 * sum(sapply(pV, sum) *( m2_beta - m1_beta^2))
}

#####################

## omega's updates ##

#####################


get_omega <- function(E1, S, Omega, lambda, n, p) {
  for (j in 1:p) {
    IOmega_nj_nj <-
      solve(Omega[-j, -j, drop = FALSE]) # implement update based on Sigma to avoid inverting here.

    s_j_j <- S[j, j]

    Omega[-j, j] <-
      Omega[j, -j] <-
      -solve((s_j_j + lambda) * IOmega_nj_nj + diag(E1[-j, j]), S[-j, j])
    Omega[j, j] <-
      Omega[j, -j, drop = FALSE] %*% IOmega_nj_nj %*% Omega[-j, j] + n / (lambda + s_j_j)

  }

  Omega

}


#######################################################

## delta's updates in GM with beta prior on edge PIP ##

#######################################################

update_m_delta_betap <- function(Omega, m_tau, m_log_rho, m_log_one_minus_rho, v0, v1, c=1) {

  1 / (1 + exp(
    c * log(v1 / v0) +
      c * m_tau * Omega ^ 2 * (1 / v1 ^ 2 - 1 / v0 ^ 2) / 2 +
      c * m_log_one_minus_rho -
      c * m_log_rho
  ))

}

####################

## rho's updates ##

####################

update_alpha_rho <- function(m_delta, a_rho, c = 1){
  c * (sum(m_delta[upper.tri(m_delta)]) + a_rho - 1) + 1
}

update_beta_rho <- function(m_delta, b_rho, c = 1){
  c * (sum(1 - m_delta[upper.tri(m_delta)]) + b_rho - 1) + 1
}


