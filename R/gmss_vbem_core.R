# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm
#
# Internal functions to call the variational EM algorithm for GMSS

#' @importFrom matrixcalc is.symmetric.matrix is.positive.definite
#' @importFrom Matrix nearPD
gmss_vbem_core  <- function(Y,
                            V,
                            list_hyper = NULL,
                            list_init = NULL,
                            tol = 1e-1,
                            maxit = 1e3,
                            verbose = T,
                            track_ELBO = F,
                            debug = F
) {

  # Disabled options
  #
  anneal <- NULL

  # Anneal
  #
  if (!is.null(anneal)) {

    vec_c <- get_annealing_ladder_(anneal, verbose = verbose)
    c <- 1

  } else {

    vec_c <- NULL
    c <- 1

  }

  # Save inputs
  #
  args <-
    list(
      Y = Y,
      V = V,
      list_hyper = list_hyper,
      list_init = list_init,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      track_ELBO = track_ELBO,
      debug = debug
    )

  # Time
  #
  pt <- Sys.time()

  # Dimension
  #
  P <- ncol(Y)
  N <- nrow(Y)
  Q <- ncol(V)

  #
  #
  S <- crossprod(Y)

  if (verbose) cat("== Preparing the hyperparameters ... \n\n")

  list_hyper <- set_default(list_hyper, 'lambda', 2)
  if(list_hyper$lambda <= 0)stop("lambda must be positive.")

  list_hyper <- set_default(list_hyper, 'v0', 0.1)
  if(list_hyper$v0 <= 0)stop("v0 must be positive.")


  list_hyper <- set_default(list_hyper, 'v1', 100)
  if(list_hyper$v1 <= 0)stop("v1 must be positive.")
  if(list_hyper$v0 >= list_hyper$v1)stop("v1 should be much greater than v0.")

  list_hyper <- set_default(list_hyper, 'a_tau', 2)
  if(list_hyper$a_tau <= 0)stop("a_tau must be positive.")

  list_hyper <- set_default(list_hyper, 'b_tau', 2)
  if(list_hyper$b_tau <= 0)stop("b_tau must be positive.")

  list_hyper <- set_default(list_hyper, 'a_sigma', 2)
  if(list_hyper$a_sigma <= 0)stop("a_sigma must be positive.")

  list_hyper <- set_default(list_hyper, 'b_sigma', 2)
  if(list_hyper$b_sigma <= 0)stop("b_sigma must be positive.")

  list_hyper <- set_default(list_hyper, 'n0', -2)

  list_hyper <- set_default(list_hyper, 't02', 0.5)
  if(list_hyper$t02 <= 0)stop("t02 must be positive.")

  list_hyper <- set_default(list_hyper, 'a_o', 1)
  if(list_hyper$a_o <= 0)stop("a_o must be positive.")

  list_hyper <- set_default(list_hyper, 'b_o', Q)
  if(list_hyper$b_o <= 0)stop("b_o must be positive.")

  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Preparing the parameter initialisation ... \n\n")

  # Specify initial parameters unless provided
  #
  # list2env(list_hyper, envir = .GlobalEnv)

  # copy to global env
  lambda <- list_hyper$lambda
  v0 <- list_hyper$v0
  v1 <- list_hyper$v1
  a_tau <- list_hyper$a_tau
  b_tau <- list_hyper$b_tau
  a_sigma <- list_hyper$a_sigma
  b_sigma <- list_hyper$b_sigma
  n0 <- list_hyper$n0
  t02 <- list_hyper$t02
  a_o <- list_hyper$a_o
  b_o <- list_hyper$b_o
  #

  # if(!is.null(list_init))list2env(list_init, envir = .GlobalEnv)

  if (!'Omega' %in% names(list_init) | ('Omega' %in% names(list_init) & is.null(list_init$Omega))) {

    Sigma <- (S + v0 * diag(P)) / N
    list_init$Omega <- as.matrix(Matrix::nearPD(solve(Sigma))$mat)

  }else{

    if(nrow(list_init$Omega)!=ncol(list_init$Omega)){
      stop("Omega should be initialised as a square matrix")
    }else if(matrixcalc::is.symmetric.matrix(list_init$Omega)){
      stop("Omega should be initialised as a symmetric matrix")
    }else if(matrixcalc::is.positive.definite(list_init$Omega)){
      stop("Omega should be initialised as a positive definite matrix")
    }

  }

  if (!'mu_beta' %in% names(list_init) | ('mu_beta' %in% names(list_init) & is.null(list_init$mu_beta))) {

    list_init$mu_beta <- rep(0,Q)

  }else{

    if(length(list_init$mu_beta) != ncol(V)){
      stop("Length of mu_beta does not match the number of columns in V.")
    }

  }

  if (!'sig2_inv_beta' %in% names(list_init) | ('sig2_inv_beta' %in% names(list_init) & is.null(list_init$sig2_inv_beta))) {

    list_init$sig2_inv_beta <- rep(1,Q)

    if(length(list_init$sig2_inv_beta) != ncol(V)){
      stop("Length of sig2_inv_beta does not match the number of columns in V.")
    }else if(any(list_init$sig2_inv_beta <=0 )){
      stop("sig2_inv_beta must be positive.")
    }
  }

  if (!'mu_zeta' %in% names(list_init) | ('mu_zeta' %in% names(list_init) & is.null(list_init$mu_zeta))) {

    list_init$mu_zeta <- list_hyper$n0

  }

  if (!'sig2_inv_zeta' %in% names(list_init) | ('sig2_inv_zeta' %in% names(list_init) & is.null(list_init$sig2_inv_zeta))) {

    list_init$sig2_inv_zeta <- 1/list_hyper$t02

  }else{

    if( list_init$sig2_inv_zeta <= 0){
      stop("sig2_inv_zeta must be positive.")
    }

  }

  if (!'alpha_sigma' %in% names(list_init) | ('alpha_sigma' %in% names(list_init) & is.null(list_init$alpha_sigma))) {

    list_init$alpha_sigma <- 1

  }else{

    if( list_init$alpha_sigma <= 0){
      stop("alpha_sigma must be positive.")
    }

  }

  if (!'beta_sigma' %in% names(list_init) | ('beta_sigma' %in% names(list_init) & is.null(list_init$beta_sigma))) {

    list_init$beta_sigma <- 1

  }else{

    if( list_init$beta_sigma <= 0){
      stop("beta_sigma must be positive.")
    }

  }

  if (!'alpha_tau' %in% names(list_init) | ('alpha_tau' %in% names(list_init) & is.null(list_init$alpha_tau))) {

    list_init$alpha_tau <- 1

  }else{

    if( list_init$alpha_tau <= 0){
      stop("alpha_tau must be positive.")
    }

  }

  if (!'beta_tau' %in% names(list_init) | ('beta_tau' %in% names(list_init) & is.null(list_init$beta_tau))) {

    list_init$beta_tau <- 1

  }else{

    if( list_init$beta_tau <= 0){
      stop("beta_tau must be positive.")
    }

  }

  if (!'alpha_o' %in% names(list_init) | ('alpha_o' %in% names(list_init) & is.null(list_init$alpha_o))) {

    list_init$alpha_o <- 1

  }else{

    if( list_init$alpha_o <= 0){
      stop("alpha_o must be positive.")
    }

  }

  if (!'beta_o' %in% names(list_init) | ('beta_o' %in% names(list_init) & is.null(list_init$beta_o))) {

    list_init$beta_o <- Q

  }else{

    if( list_init$beta_o <= 0){
      stop("beta_o must be positive.")
    }

  }

  #
  Omega <- list_init$Omega
  mu_zeta <- list_init$mu_zeta
  sig2_inv_zeta <- list_init$sig2_inv_zeta
  mu_beta <- list_init$mu_beta
  sig2_inv_beta <- list_init$sig2_inv_beta
  alpha_sigma <- list_init$alpha_sigma
  beta_sigma <- list_init$beta_sigma
  alpha_tau <- list_init$alpha_tau
  beta_tau <- list_init$beta_tau
  alpha_o <- list_init$alpha_o
  beta_o <- list_init$beta_o
  #

  # Initialise deduced quantities
  #
  m_tau <- get_m_tau(alpha_tau, beta_tau)

  m_log_o <- get_m_log_o(alpha_o, beta_o)
  m_log_one_minus_o <- get_m_log_one_minus_o(alpha_o, beta_o)
  m_log_sig2_inv <- get_m_log_sig2_inv(alpha_sigma, beta_sigma)
  m_sig2_inv <- get_m_sig2_inv(alpha_sigma, beta_sigma)
  m_gamma <-
    update_m_gamma(
      m_log_o,
      m_log_one_minus_o,
      m_log_sig2_inv,
      mu_beta,
      sig2_inv_beta,
      c
    )

  m1_beta <- get_m1_beta(mu_beta, m_gamma)
  m2_beta <- get_m2_beta(mu_beta, sig2_inv_beta, m_gamma)

  theta <- get_theta(m1_beta, V)
  m1_alpha <- get_Alpha(theta, mu_zeta, P)

  if (verbose) cat("... done. == \n\n")

  # track ELBO after each maximisation step
  #
  if (track_ELBO) {

    vec_ELBO_M <- c()

  } else{

    vec_ELBO_M <- NA

  }

  # debug mode
  #
  if (debug) {

    # track ELBO within each variational step
    list_ELBO <- list()

    # record number of warnings
    n_warning <- 0
    vec_n_warning_VB <- c()

  }

  # Pre-compute
  #
  sV <- list_upper_tri_matrix_sum_Viq_Vjq(V , P, Q)
  pV <- list_upper_tri_matrix_prod_Viq_Vjq(V, P, Q)
  spV1 <- list_vec_sum_prod_Viq_Vinq(V, P, Q)
  spV2 <- list_vec_sum_prod_Viq_Vjnq_prod_Vinq_Vjq(V, P, Q)

  #
  eps <- .Machine$double.eps

  #
  it <- 0
  vec_VB_it <- c()
  ELBO_M_diff <- Inf
  ELBO_M_old <- -Inf

  #
  while ((ELBO_M_diff > tol & (it < maxit))) {

    # Iteration
    #
    it <- it + 1

    if (verbose & it %% 5 == 0)
      cat(paste0("Iteration ", format(it), "... \n"))

    # Anneal
    #
    if (!is.null(vec_c) && it <= length(vec_c)) {
      vbc <- vec_c[it]
      print(paste0("Temperature: ", format(1 / vbc, digits = 3)))
    } else {
      vbc <- 1
    }


    # VBE step :
    # ======== #

    ELBO_diff <- Inf
    ELBO_old <- -Inf
    VB_it <- 0

    if(debug){
      n_warning_VB <- 0
      tmp <- c()
    }

    while (ELBO_diff > tol & VB_it < maxit) {

      #
      VB_it <- VB_it + 1

      if (verbose != 0 & VB_it %% 5 == 0)
        cat(paste0("VBE iteration ", format(VB_it), "... \n"))

      # % # m_delta
      m_delta <- update_m_delta(Omega, m_tau, m1_alpha, v0, v1, vbc)
      E1 <- get_E1(m_delta, v0, v1)
      # % #

      # % # z
      m1_z <- get_E2(m_delta, m1_alpha, vbc)
      m2_z <- get_E2_2(m1_z, m1_alpha, vbc)
      # % #

      # % # o
      alpha_o <- update_alpha_o(m_gamma, a_o, vbc)
      beta_o <- update_beta_o(m_gamma, b_o, vbc)

      m_o <- get_m_o(alpha_o, beta_o)
      m_log_o <- get_m_log_o(alpha_o, beta_o)
      m_log_one_minus_o <- get_m_log_one_minus_o(alpha_o, beta_o)
      # % #

      # % # tau
      alpha_tau <- update_alpha_tau(a_tau, P, vbc)
      beta_tau <- update_beta_tau(Omega, E1, b_tau, vbc)

      m_tau <- get_m_tau(alpha_tau, beta_tau)
      m_log_tau <- get_m_log_tau(alpha_tau, beta_tau)
      # % #

      # % # sigma
      alpha_sigma <- update_alpha_sigma(m_gamma, a_sigma, vbc)
      beta_sigma <- update_beta_sigma(m_gamma, m2_beta, b_sigma, vbc)

      m_log_sig2_inv <- get_m_log_sig2_inv(alpha_sigma, beta_sigma)
      m_sig2_inv <- get_m_sig2_inv(alpha_sigma, beta_sigma)
      # % #

      # % # zeta
      sig2_inv_zeta <- update_sig2_inv_zeta(t02, P, vbc)
      mu_zeta <-
        update_mu_zeta(m1_beta, sig2_inv_zeta, m1_z, V, n0, t02, P, vbc)
      m2_zeta <- get_m2_zeta(mu_zeta, sig2_inv_zeta)
      # % #

      # % # beta
      sig2_inv_beta <- update_sig2_inv_beta(m_sig2_inv, V, pV, P, vbc)
      mu_beta <- update_mu_beta(sig2_inv_beta, m1_beta, mu_zeta, m1_z,
                                V, sV, spV1, spV2, P, Q, vbc)
      # % #

      # % # gamma
      m_gamma <-
        update_m_gamma(
          m_log_o,
          m_log_one_minus_o,
          m_log_sig2_inv,
          mu_beta,
          sig2_inv_beta,
          vbc
        )
      # % #

      # % #
      m1_beta <- get_m1_beta(mu_beta, m_gamma)
      m2_beta <- get_m2_beta(mu_beta, sig2_inv_beta, m_gamma)
      # % #

      # % #
      theta <- get_theta(m1_beta, V)
      m1_alpha <- get_Alpha(theta, mu_zeta, P)
      sum_var_alpha <- get_sum_var_alpha(sig2_inv_zeta, m1_beta, m2_beta, V, pV, P)
      # % #

      # % #
      ELBO <-  get_elbo_gmss_vbem(Omega,
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
                                  vbc)
      ELBO_diff <- abs(ELBO - ELBO_old)

      # check the increasing ELBO in the variational step
      #
      if (debug && ELBO + eps < ELBO_old) {

        warning(paste0(
          "Non-increasing in the VB step: ELBO_old = ",
          ELBO_old,
          ", ELBO = ",
          ELBO,
          '\n'
        ))
        n_warning_VB <- n_warning_VB + 1
      }

      if (verbose & VB_it %% 5 == 0)
        cat(paste0(
          "Difference ELBO from previous iteration: ",
          format(ELBO_diff),
          "\n"
        ))

      ELBO_old <- ELBO

      if (debug) {
        tmp <- c(tmp, ELBO)
      }
      # % #
    }


    # % #
    if (debug) {
      list_ELBO <- c(list_ELBO, list(tmp))
      vec_n_warning_VB <- c(vec_n_warning_VB, n_warning_VB)
    }

    vec_VB_it <- c(vec_VB_it, VB_it)

    if (VB_it == maxit) {
      warning('Maximal number of iterations reached before convergence in VBE step. \n')
    }
    # % #

    # M step :
    # ====== #
    #  check the increasing ELBO in the M-step
    #
    if (debug) {
      ELBO_M0 <- get_elbo_gmss_vbem(Omega,
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
                                    vbc)
    }

    # % # Omega
    #
    # bool_cpp <- F
    # if (bool_cpp) {
    #   bool_direct_solve <-  F # keep F as otherwise Omega inverted at each iteration
    #   if (bool_direct_solve) {
    #     out <- M_Omega_direct_solve(N, P, Omega, S, lambda, m_tau * E1)
    #     Omega <- out$Omega
    #   } else {
    #     out <- M_Omega(N, P, Sigma, Omega, S, lambda, m_tau * E1)
    #     Omega <- out$Omega
    #     Sigma <- out$Sigma
    #   }
    # } else {
      Omega <- get_omega(m_tau * E1, S, Omega, lambda, N, P)
    # }
    # % #

    #
    ELBO_M <- get_elbo_gmss_vbem(Omega,
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
                                 vbc)

    ELBO_M_diff <- abs(  ELBO_M -   ELBO_M_old)

    if (debug &&   ELBO_M + eps <   ELBO_M0) {

      warning(paste0("Non-increasing in the M-step :   ELBO_0 = ",   ELBO_M0, ",   ELBO = ",   ELBO_M))
      n_warning <- n_warning + 1

    }

    if (verbose & it %% 5 == 0)
      cat(paste0(
        "Difference ELBO_M from previous iteration: ",
        format(ELBO_M_diff),
        "\n"
      ))

    ELBO_M_old <- ELBO_M

    if (track_ELBO) {
      vec_ELBO_VBEM <- c(vec_ELBO_VBEM, ELBO_M)
    }
  }


  if(ELBO_M_diff <= tol){
    if(verbose)
      cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                 "Optimal marginal log-likelihood variational lower bound ",
                 "(ELBO) = ", format(ELBO_M), ". \n\n"))
  }

  if (it == maxit) {
    warning('Maximal number of iterations reached before convergence. Exit.')
  }

  pt <- Sys.time() - pt
  cat('Algorithm runtime: ',format(pt), '\n')

  estimates <- list( Omega = Omega,
                     m_delta = m_delta,
                     alpha_tau = alpha_tau,
                     beta_tau = beta_tau,
                     m_gamma =  m_gamma,
                     mu_beta = mu_beta,
                     sig2_inv_beta = sig2_inv_beta,
                     mu_zeta = mu_zeta,
                     sig2_inv_zeta = sig2_inv_zeta,
                     alpha_o = alpha_o,
                     beta_o = beta_o,
                     alpha_sigma = alpha_sigma,
                     beta_sigma = beta_sigma,
                     S = S # for model comparison
  )

  if(debug){

    debugs <- list( n_warning = n_warning,
                    vec_n_warning_VB = vec_n_warning_VB,
                    list_ELBO = list_ELBO)

  }else{

    debugs <- NA

  }

  create_named_list_(
    args,
    estimates,
    debugs,
    it,
    vec_VB_it,
    vec_ELBO_M,
    pt
  )
}


