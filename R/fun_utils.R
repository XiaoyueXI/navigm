# check the presence of character name in the list x; if not, set to default.
#
set_default <- function(x, name, default){
  if (!name %in% names(x) | (name %in% names(x) & is.null(x[[name]]))) {
    x[[name]] <- default
  }
  return(x)
}

get_theta <- function(beta, V) {
  V %*% beta
}

get_Alpha <- function(theta, zeta, p) {

  M_theta <- matrix(theta, nrow = p, ncol = p)
  Alpha <- zeta + M_theta + t(M_theta)
  diag(Alpha) <- 0 # not used

  Alpha

}

list_upper_tri_matrix_sum_Viq_Vjq <- function(V , p, q) {
  # V_{iq} + V_{jq}, i<j
  # q lists of p-1 times p-1 upper triangular matrix

  ans <- lapply(1:q, function(j) {
    sv <- matrix(V[1:(p - 1), j], nrow = p - 1, ncol = p - 1)
    sv <-
      sv + t(matrix(V[p:2, j], nrow = p - 1, ncol = p - 1))[, ncol(sv):1]
    sv[lower.tri(sv)] <- 0
    return(sv)
  })

  ans
}


list_vec_sum_prod_Viq_Vjnq_prod_Vinq_Vjq <- function(V, p, q) {
  # sum_{i<j} (ViqVjq' + Viq'Vjq)
  # q lists of (q-1) vector

  Map('+',
      lapply(1:q, function(k){
        Reduce('+',lapply(1:(p-1), function(i){
          apply(matrix(V[i,k] * V[(i+1):p, -k],ncol = q-1), 2, sum)
        }))
      }),
      lapply(1:q, function(k){
        Reduce('+',lapply(2:p, function(j){
          apply(matrix(V[(1:j-1),-k] * V[j,k],ncol = q-1), 2, sum)
        }))
      })
  )

}


list_upper_tri_matrix_prod_Viq_Vjq <- function(V , p, q) {
  # sum_{i<j} Viq Vjq
  # q lists of p-1 times p-1 upper triangular matrix

  ans <- lapply(1:q, function(j) {
    pv <- matrix(V[1:(p - 1), j], nrow = p - 1, ncol = p - 1)
    pv  <-
      pv * t(matrix(V[p:2, j], nrow = p - 1, ncol = p - 1))[, ncol(pv):1]
    pv[lower.tri(pv)] <- 0
    return(pv)
  })

  ans

}


list_vec_sum_prod_Viq_Vinq <- function(V, p, q) {
  # sum_i Viq Viq'
  # q lists of q-1 vector

  ans <- lapply(1:q, function(i) {
    apply(matrix(V[, i], nrow = p, ncol = (q - 1)) * V[,-i], 2, sum)
  })

  ans

}

inv_mills_ratio_ <- function(delta, U, log_1_pnorm_U, log_pnorm_U) {
  stopifnot(delta %in% c(0, 1))

  # writing explicitely the formula for pnorm(, log = TRUE) is faster...
  if (delta == 1) {
    m <- exp(-U ^ 2 / 2 - log(sqrt(2 * pi)) - log_pnorm_U)

    m[m < -U] <- -U[m < -U] # to do correct in other packages

  } else {
    m <- -exp(-U ^ 2 / 2 - log(sqrt(2 * pi)) - log_1_pnorm_U)

    m[m > -U] <- -U[m > -U]

  }

  m

}

get_annealing_ladder_ <- function(anneal, verbose) {
  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10, arxiv v4

  k_m <- 1 / anneal[2]
  m <- anneal[3]

  if (anneal[1] == 1) {
    type <- "geometric"

    delta_k <- k_m ^ (1 / (1 - m)) - 1

    ladder <- (1 + delta_k) ^ (1 - m:1)

  } else if (anneal[1] == 2) {
    # harmonic spacing

    type <- "harmonic"

    delta_k <- (1 / k_m - 1) / (m - 1)

    ladder <- 1 / (1 + delta_k * (m:1 - 1))

  } else if (anneal[1] == 3) {
    # linear spacing

    type <- "linear"

    delta_k <- (1 - k_m) / (m - 1)

    ladder <- k_m + delta_k * (1:m - 1)
  } else {
    type <- "fixed"
    ladder <- k_m
  }

  if (verbose != 0)
    cat(paste0("** Annealing with ", type, " spacing ** \n\n"))

  ladder

}


create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))

}

log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))

  if (max(abs(x)) > max(x)) {
    offset <- min(x)
  } else {
    offset <- max(x)
  }

  log(sum(exp(x - offset))) + offset

}

effective_sum <- function(x) {

  # sum(x)
  xp <- x[x > 0]
  xn <- x[x < 0]

  if (length(xp) != 0 & length(xn) != 0) {

    exp(log_sum_exp(log(xp))) - exp(log_sum_exp(log(-xn)))

  } else if (length(xp) == 0 & length(xn) != 0) {

    -exp(log_sum_exp(log(-xn)))

  } else if (length(xp) != 0 & length(xn) == 0) {

    exp(log_sum_exp(log(xp)))

  } else if (length(xp) == 0 & length(xn) == 0) {

    0

  }

}

# Functions for hyperparameter settings
#
E_Phi_X <- function(mu, s2, lower_tail = TRUE) {

  pnorm(mu / sqrt(1 + s2), lower.tail = lower_tail)

}

E_Phi_X_2 <- function(mu, s2) {

  pnorm(mu / sqrt(1 + s2)) -
    2 * PowerTOST::OwensT(mu / sqrt(1 + s2), 1 / sqrt(1 + 2 * s2))

}

get_V_p_t <- function(mu, s2, p) {
  p * (p - 1) * E_Phi_X_2(mu, s2) -
    p^2 * E_Phi_X(mu, s2)^2 +
    p * E_Phi_X(mu, s2)
}


get_mu <- function(E_p_t, s2, p) {

  sqrt(1 + s2) * qnorm(E_p_t / p)

}


get_n0_t02 <- function(p, p_star) {

  E_p_t <- p_star[1]
  V_p_t <- min(p_star[2], floor(2 * p / 3))

  dn <- 1e-6
  up <- 1e5

  # Get n0 and t02 similarly as for a_omega_t and b_omega_t in HESS
  # (specify expectation and variance of number of active predictors per response)
  #
  # Look at : gam_st | theta_s = 0
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of edges when no hubs appear. ",
                  "Please change their values."))
    })

  # n0 sets the level of sparsity.
  n0 <- get_mu(E_p_t, t02, p)

  create_named_list_(n0, t02)
}



# Functions for hyperparameter settings
#
E_Phi_X <- function(mu, s2, lower_tail = TRUE) {

  pnorm(mu / sqrt(1 + s2), lower.tail = lower_tail)

}

E_Phi_X_2 <- function(mu, s2) {

  pnorm(mu / sqrt(1 + s2)) -
    2 * PowerTOST::OwensT(mu / sqrt(1 + s2), 1 / sqrt(1 + 2 * s2))

}

get_V_p_t <- function(mu, s2, p) {
  p * (p - 1) * E_Phi_X_2(mu, s2) -
    p^2 * E_Phi_X(mu, s2)^2 +
    p * E_Phi_X(mu, s2)
}


get_mu <- function(E_p_t, s2, p) {

  sqrt(1 + s2) * qnorm(E_p_t / p)

}


get_n0_t02 <- function(p, p_star) {

  E_p_t <- p_star[1]
  V_p_t <- min(p_star[2], floor(2 * p / 3))

  dn <- 1e-6
  up <- 1e5

  # Get n0 and t02 similarly as for a_omega_t and b_omega_t in HESS
  # (specify expectation and variance of number of active predictors per response)
  #
  # Look at : gam_st | theta_s = 0
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of edges when no hubs appear. ",
                  "Please change their values."))
    })

  # n0 sets the level of sparsity.
  n0 <- get_mu(E_p_t, t02, p)

  create_named_list_(n0, t02)
}



