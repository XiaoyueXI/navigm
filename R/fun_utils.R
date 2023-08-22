# check the presence of character name in the list x; if not, set to default.
#
set_default <- function(x, name, default){
  if (!name %in% names(x) | (name %in% names(x) & is.null(x[[name]]))) {
    x[[name]] <- default
  }
  return(x)
}

# deduced quantities
#
get_theta <- function(beta, V) {
  V %*% beta
}

get_Alpha <- function(theta, zeta, p) {

  M_theta <- matrix(theta, nrow = p, ncol = p)
  Alpha <- zeta + M_theta + t(M_theta)
  diag(Alpha) <- 0 # not used

  Alpha

}

inv_mills_ratio_ <- function(delta, U, log_1_pnorm_U, log_pnorm_U) {
  stopifnot(delta %in% c(0, 1))

  # writing explicitly the formula for pnorm(, log = TRUE) is faster...
  if (delta == 1) {
    m <- exp(-U ^ 2 / 2 - log(sqrt(2 * pi)) - log_pnorm_U)

    m[m < -U] <- -U[m < -U] # to do correct in other packages

  } else {
    m <- -exp(-U ^ 2 / 2 - log(sqrt(2 * pi)) - log_1_pnorm_U)

    m[m > -U] <- -U[m > -U]

  }

  m

}

# V quantities
#
list_upper_tri_matrix_sum_Viq_Vjq <- function(V, p, q) {
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


# may not be used
#
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

#' @importFrom stats setNames
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


#' @importFrom stats pnorm
E_Phi_X <- function(mu, s2, lower_tail = TRUE) {

  stats::pnorm(mu / sqrt(1 + s2), lower.tail = lower_tail)

}

#' @importFrom PowerTOST OwensT
#' @importFrom stats pnorm
E_Phi_X_2 <- function(mu, s2) {

  stats::pnorm(mu / sqrt(1 + s2)) -
    2 * PowerTOST::OwensT(mu / sqrt(1 + s2), 1 / sqrt(1 + 2 * s2))

}

get_V_p_t <- function(mu, s2, p) {
  p * (p - 1) * E_Phi_X_2(mu, s2) -
    p^2 * E_Phi_X(mu, s2)^2 +
    p * E_Phi_X(mu, s2)
}

#' @importFrom stats qnorm
get_mu <- function(E_p_t, s2, p) {

  sqrt(1 + s2) * stats::qnorm(E_p_t / p)

}

#' @importFrom stats uniroot
get_n0_t02 <- function(p, p_star) {

  E_p_t <- p_star[1]
  # V_p_t <- min(p_star[2], floor(2 * p / 3))
  V_p_t <- p_star[2]

  dn <- 1e-6
  up <- 1e5

  # Get n0 and t02
  #
  tryCatch(t02 <- stats::uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of edges when no hubs appear. ",
                  "Please change their values."))
    })

  #
  n0 <- get_mu(E_p_t, t02, p)

  create_named_list_(n0, t02)
}


#' Simulate precision matrices and data given an adjacency matrix.
#'
#' This function simulates precision matrices and the observations based on
#' the pre-specified adjacency matrix,
#' inspired by Learning Graphical Models With Hubs, JMLR, 2014, P.3307.
#'
#' @param N Scalar. Number of observations.
#'
#' @param A An adjacency matrix.
#'
#' @param vec_magnitude A vector of two positive numbers indicating the range of absolute magnitudes of precision matrix entries.
#'
#' @param bool_scale Logical. If TRUE (default), scale the samples; otherwise, not scale.
#'
#' @return A list containing the simulated data:
#'  \describe{
#'  \item{A}{the input adjacency matrix.}
#'  \item{Omega}{the simulated adjacency matrix of same size and structure of the pre-specified adjacency matrix.}
#'  \item{Y}{the simulated observations with N rows and \code{ncol(A)} columns.}
#' }
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats runif
#' @export
generate_data_from_adjancency <- function(N,
                                          A,
                                          vec_magnitude = c(0.25, 0.75),
                                          bool_scale = TRUE) {

  #
  P <- nrow(A)
  nb_edges <- sum(A == 1)

  # matrix E
  E <- A

  E[A == 1] <- stats::runif(nb_edges, min = vec_magnitude[1], max = vec_magnitude[2])

  E_bar <- (E + t(E)) / 2

  msign <- matrix(1, nrow = nrow(E), ncol = ncol(E))
  msign[upper.tri(msign)] <- sample(c(-1,1), size = sum(upper.tri(msign)),  prob = c(0.5, 0.5), replace = TRUE)
  msign[lower.tri(msign)] <- t(msign)[lower.tri(msign)]
  E_bar <- E_bar * msign

  # minimum eigenvalue
  min_eigen <- min(eigen(E_bar, only.values = TRUE)$values)

  if (min_eigen < 0) {
    Omega <- E_bar + (0.1 - min_eigen) * diag(P)
  } else{
    Omega <- E_bar + 0.1 * diag(P)
  }

  #
  Y <- mvtnorm::rmvnorm(N, rep(0, P), solve(Omega))

  if (bool_scale) {
    Y <- scale(Y)
  }

  create_named_list_(A, Omega, Y)

}

#' Simulate beta distributed auxiliary variables.
#'
#' This function simulates auxiliary variable matrices that mimic the PPIs in a Bayesian spike-and-slab regression.
#'
#' @param P Scalar. Number of nodes in the graph.
#' @param Q Scalar. Number of node-level auxiliary variables.
#' @param alpha,beta Scalars. Shape parameters of beta distribution with default 0.05 and 0.2.
#' @param Sigma A correlation matrix of auxiliary variables.
#' @param min_gene Scalar. Number of influential entries per auxiliary variable.
#' @param verbose Logical. If FALSE (default), not show messages; otherwise, show the messages.
#'
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm cor qnorm runif
#' @export
generate_V <- function(P, Q,
                       alpha = 0.05, beta = 0.2,
                       Sigma = diag(Q),
                       min_gene = round(0.05 * P),
                       verbose = F) {
  #
  if(P < Q){
    warning('P < Q, set empirical = F. May not be comparable with those simulated under empirical = T.\n')
    Z <- MASS::mvrnorm(n = P,
                 mu = rep(0, nrow(Sigma)),
                 Sigma)

  }else{
    Z <- MASS::mvrnorm(n = P,
                 mu = rep(0, nrow(Sigma)),
                 Sigma,
                 empirical = T)
  }

  #
  bool_up <- upper.tri(cor(Z))
  V <- t(stats::qbeta(t(stats::pnorm(Z, 0, sqrt(
    diag(Sigma)
  ))), alpha, beta))

  if(verbose){
    cat('Range of empirical correlations:', range(stats::cor(V)[bool_up]), '\n')
    cat('Range of absolute empirical correlations:', range(abs(cor(V)[bool_up])), '\n')
  }


  # set min_gene
  #
  Vnz <- apply(V, 2, function(x)
    sum(x > 0.5))

  if (verbose)
    cat(sum(Vnz < min_gene),
        ' variant(s) do not have significant effects on gene expression.')

  if (any(Vnz < min_gene)) {
    for (j in which(Vnz < min_gene)) {
      V[sample(P, min_gene), j] <- stats::runif(min_gene, 0.9, 1)
    }
  }
  return(V)
}

