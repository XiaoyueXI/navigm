# This file is part of the `navigss` R package:
#     https://github.com/XiaoyueXI/navigss
#
# Internal functions gathering the ELBO terms for the core algorithms

#####################

## E log p(y|rest) ##

#####################

e_y <- function(Omega, S, N){

  N * determinant(Omega, logarithm = TRUE)$modulus[1] / 2 -
    sum(S * Omega) / 2

}

############################################

## E log p(omega | rest) - E log q(omega) ##

############################################

e_omega <- function(Omega,
                    m_delta,
                    E1,
                    m_tau,
                    m_log_tau,
                    lambda,
                    v0,
                    v1,
                    P){

  bool_up <- upper.tri(Omega)

  - lambda * sum(diag(Omega)) / 2 -
    log(v1) * sum(m_delta[bool_up]) -
    log(v0) * sum(1-m_delta[bool_up]) -
    m_tau * sum(Omega[bool_up]^2 * E1[bool_up])/2 +
    m_log_tau * P * (P-1)/4

}

##################################################

## E log p(delta, z | rest) - E log q(delta, z) ##

##################################################

e_delta_z <- function(m_delta,
                      m1_alpha,
                      sum_var_alpha,
                      c){

  eps <- .Machine$double.eps
  bool_up <- upper.tri(Omega)

  - sum_var_alpha/2 +
    #
    sum(m_delta[bool_up] * pnorm(sqrt(c) * m1_alpha[bool_up], log.p = T, lower.tail = T)) +
    sum((1-m_delta[bool_up]) * pnorm(sqrt(c) * m1_alpha[bool_up], log.p = T, lower.tail = F)) -
    #
    sum(m_delta[bool_up] * log(m_delta[bool_up] + eps))-
    sum((1-m_delta[bool_up]) * log(1-m_delta[bool_up] + eps))

}

##########################################

##  E log p(tau | rest) - E log q(tau)  ##

##########################################

# compare with b4
e_tau <- function(alpha_tau,
                  beta_tau,
                  m_tau,
                  m_log_tau,
                  a_tau,
                  b_tau){

  eps <- .Machine$double.eps

  m_log_tau * (a_tau - alpha_tau) -
    m_tau * (b_tau - beta_tau) -
    alpha_tau * log(beta_tau + eps) + lgamma(alpha_tau)

}

##########################################

## E log p(zeta | rest) - E log q(zeta) ##

##########################################

e_zeta <- function(mu_zeta,
                   sig2_inv_zeta,
                   m2_zeta,
                   n0,
                   t02){

  eps <- .Machine$double.eps
  - m2_zeta/(2 * t02) +
    n0 * mu_zeta/t02 +
    (1 + log(2 * pi * sig2_inv_zeta^(-1) + eps))/2

}

########################################################

## E log p(beta, gamma | rest) - E log q(beta, gamma) ##

########################################################

e_beta_gamma <- function(m_gamma,
                         m_sig2_inv,
                         m_log_sig2_inv,
                         m_log_o,
                         m_log_one_minus_o,
                         m2_beta,
                         sig2_inv_beta,
                         a_sigma,
                         b_sigma){

  eps <- .Machine$double.eps

  sum(m_gamma)/2 * m_log_sig2_inv -
    sum(m_gamma * m2_beta)/2 * m_sig2_inv +
    m_log_o * sum(m_gamma) +
    m_log_one_minus_o * sum(1 - m_gamma)  +
    sum(m_gamma * (1 + log(2 * pi * sig2_inv_beta^(-1) + eps))/2 ) -
    sum(m_gamma * log(m_gamma + eps))-
    sum((1-m_gamma) * log(1 - m_gamma + eps))

}

####################################

## E log p(o | rest) - E log q(o) ##

####################################

e_o <- function(alpha_o,
                beta_o,
                m_log_o,
                m_log_one_minus_o,
                a_o,
                b_o){

  m_log_o * (a_o - alpha_o) +
    m_log_one_minus_o * ( b_o - beta_o) +
    lbeta(alpha_o, beta_o)

}

############################################

## E log p(sigma | rest) - E log q(sigma) ##

############################################

e_sigma <- function(alpha_sigma,
                    beta_sigma,
                    m_sig2_inv,
                    m_log_sig2_inv,
                    a_sigma,
                    b_sigma){

  eps <- .Machine$double.eps

  (a_sigma - alpha_sigma) * m_log_sig2_inv -
    (b_sigma - beta_sigma) * m_sig2_inv -
    alpha_sigma * log(beta_sigma + eps) + lgamma(alpha_sigma)

}
#######################################

## ELBO (select auxiliary variables) ##

#######################################

get_elbo_gmss_vbem <- function(Omega,
                               m_delta,
                               alpha_tau,
                               beta_tau,
                               m_tau,
                               m_log_tau,
                               mu_zeta,
                               sig2_inv_zeta,
                               m2_zeta,
                               m_gamma,
                               alpha_sigma,
                               beta_sigma,
                               m_log_sig2_inv,
                               m_sig2_inv,
                               m2_beta,
                               sig2_inv_beta,
                               m1_alpha,
                               sum_var_alpha,
                               alpha_o,
                               beta_o,
                               m_log_o,
                               m_log_one_minus_o,
                               E1,
                               S,
                               lambda,
                               v0,
                               v1,
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
                               c) {

  c * (
    e_y(Omega, S, N) +
      e_omega(Omega,
              m_delta,
              E1,
              m_tau,
              m_log_tau,
              lambda,
              v0,
              v1,
              P) +
      e_delta_z(m_delta,
                m1_alpha,
                sum_var_alpha,
                c) +
      e_tau(alpha_tau,
            beta_tau,
            m_tau,
            m_log_tau,
            a_tau,
            b_tau) +
      e_zeta(mu_zeta,
             sig2_inv_zeta,
             m2_zeta,
             n0,
             t02) +
      e_beta_gamma(m_gamma,
                   m_sig2_inv,
                   m_log_sig2_inv,
                   m_log_o,
                   m_log_one_minus_o,
                   sig2_inv_beta,
                   a_sigma,
                   b_sigma) +
      e_o(alpha_o,
          beta_o,
          m_log_o,
          m_log_one_minus_o,
          a_o,
          b_o) +
      e_sigma(alpha_sigma,
              beta_sigma,
              m_sig2_inv,
              m_log_sig2_inv,
              a_sigma,
              b_sigma)
  )
}
