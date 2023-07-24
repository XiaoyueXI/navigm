# This file is part of the `navigss` R package:
#     https://github.com/XiaoyueXI/navigss
#
# Internal functions gathering the ELBO terms for the ECM algorithms


get_elbo_gm_em_v1 <- function(Omega,
                              rho,
                              tau1,
                              P1,
                              E1,
                              S,
                              lambda,
                              v0,
                              v1,
                              a_rho,
                              b_rho,
                              a_tau,
                              b_tau,
                              N,
                              P) {

  #
  eps <- .Machine$double.eps
  bool_up <- upper.tri(Omega)

  #
  N * determinant(Omega, logarithm = TRUE)$modulus[1] / 2 -
    #
    #= sum(S * t(Omega)) since Omega is symmetric tr(Y^T Y Omega)
    #
    sum(S * Omega) / 2 -
    lambda * sum(diag(Omega)) / 2 -
    sum(log(v1) * P1[bool_up] + log(v0) * (1 - P1[bool_up]))  +
    #
    (P * (P - 1) / 2) * (log(tau1 + eps) / 2) -
    tau1 * sum(Omega[bool_up] ^ 2 * E1[bool_up]) / 2 +
    (a_tau - 1) * log(tau1 + eps) - b_tau * tau1 +
    #
    sum(a_rho - 1 + P1[bool_up]) * log(rho + eps) +
    sum(b_rho - P1[bool_up]) * log(1 - rho + eps)
}


get_elbo_gm_em_v2 <- function(Omega,
                              zeta,
                              tau1,
                              P1,
                              E1,
                              E2,
                              E2_2,
                              S,
                              lambda,
                              v0,
                              v1,
                              a_rho,
                              b_rho,
                              a_tau,
                              b_tau,
                              N,
                              P) {

  #
  eps <- .Machine$double.eps
  bool_up <- upper.tri(Omega)

  #
  theta <- rep(0, P)
  mat_tmp <- get_Alpha(theta, zeta, P)


  #
  N * determinant(Omega, logarithm = TRUE)$modulus[1] / 2 -
    #
    #= sum(S * t(Omega)) since Omega is symmetric tr(Y^T Y Omega)
    #
    sum(S * Omega) / 2 -
    lambda * sum(diag(Omega)) / 2 -
    sum(log(v1) * P1[bool_up] + log(v0) * (1 - P1[bool_up]))  +
    #
    (P * (P - 1) / 2) * (log(tau1 + eps) / 2) -
    tau1 * sum(Omega[bool_up] ^ 2 * E1[bool_up]) / 2 +
    (a_tau - 1) * log(tau1 + eps) - b_tau * tau1 -
    #
    sum( E2_2[bool_up] ) / 2 +
    sum( mat_tmp[bool_up] * E2[bool_up] ) -
    sum( mat_tmp[bool_up] ^ 2 ) / 2 -
    #
    ( zeta ^ 2 - 2 * n0 * zeta ) / ( 2 * t02 )

}

get_elbo_gmn_em <- function(Omega,
                            zeta,
                            beta,
                            tau1,
                            tau2,
                            P1,
                            E1,
                            E2,
                            E2_2,
                            S,
                            V,
                            lambda,
                            v0,
                            v1,
                            n0,
                            t02,
                            a_tau,
                            b_tau,
                            a_sigma,
                            b_sigma,
                            N,
                            P,
                            Q) {

  #
  eps <- .Machine$double.eps
  bool_up <- upper.tri(Omega)

  #
  theta <- get_theta(beta, V)
  mat_tmp <- get_Alpha(theta, zeta, P)


  #
  N * determinant(Omega, logarithm = TRUE)$modulus[1] / 2 -
    #
    #= sum(S * t(Omega)) since Omega is symmetric tr(Y^T Y Omega)
    #
    sum(S * Omega) / 2 -
    lambda * sum(diag(Omega)) / 2 -
    sum(log(v1) * P1[bool_up] + log(v0) * (1 - P1[bool_up]))  +
    #
    (P * (P - 1) / 2) * (log(tau1 + eps) / 2) -
    tau1 * sum(Omega[bool_up] ^ 2 * E1[bool_up]) / 2 +
    (a_tau - 1) * log(tau1 + eps) - b_tau * tau1 -
    #
    sum( E2_2[bool_up] ) / 2 +
    sum( mat_tmp[bool_up] * E2[bool_up] ) -
    sum( mat_tmp[bool_up] ^ 2 ) / 2 -
    #
    ( zeta ^ 2 - 2 * n0 * zeta ) / ( 2 * t02 ) +
    #
    Q * log(tau2 + eps) / 2 -
    tau2 * sum( beta ^ 2 ) / 2 +
    #
    (a_sigma - 1) * log(tau2 + eps) - b_sigma * tau2

}


get_elbo_gmss_em <- function(Omega,
                             zeta,
                             beta,
                             o,
                             tau1,
                             tau2,
                             P1,
                             P2,
                             E1,
                             E2,
                             E2_2,
                             E5,
                             S,
                             V,
                             lambda,
                             v0,
                             v1,
                             s0,
                             s1,
                             n0,
                             t02,
                             a_tau,
                             b_tau,
                             a_sigma,
                             b_sigma,
                             a_o,
                             b_o,
                             N,
                             P,
                             Q) {

  #
  eps <- .Machine$double.eps
  bool_up <- upper.tri(Omega)

  #
  theta <- get_theta(beta, V)
  mat_tmp <- get_Alpha(theta, zeta, P)


  #
  N * determinant(Omega, logarithm = TRUE)$modulus[1] / 2 -
    #
    #= sum(S * t(Omega)) since Omega is symmetric tr(Y^T Y Omega)
    #
    sum(S * Omega) / 2 -
    lambda * sum(diag(Omega)) / 2 -
    sum(log(v1) * P1[bool_up] + log(v0) * (1 - P1[bool_up]))  +
    #
    (P * (P - 1) / 2) * (log(tau1 + eps) / 2) -
    tau1 * sum(Omega[bool_up] ^ 2 * E1[bool_up]) / 2 +
    (a_tau - 1) * log(tau1 + eps) - b_tau * tau1 -
    #
    sum( E2_2[bool_up] ) / 2 +
    sum( mat_tmp[bool_up] * E2[bool_up] ) -
    sum( mat_tmp[bool_up] ^ 2 ) / 2 -
    #
    ( zeta ^ 2 - 2 * n0 * zeta ) / ( 2 * t02 ) -
    #
    sum(log(s1) * P2 + log(s0) * ( 1 - P2 ))  +
    #
    Q * log(tau2 + eps) / 2 -
    tau2 * sum( beta ^ 2 * E5 ) / 2 +
    #
    (sum(P2) + a_o - 1) * log(o + eps) +
    (sum(1 - P2) + b_o - 1) * log(1 - o + eps) +
    #
    (a_sigma - 1) * log(tau2 + eps) - b_sigma * tau2

}
