# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm
#
# Internal functions to call the variational EM algorithm for GM

#' @importFrom matrixcalc is.symmetric.matrix is.positive.definite
#' @importFrom Matrix nearPD
gm_vbem_core  <- function(Y,
                          list_hyper = NULL,
                          list_init = NULL,
                          tol = 1e-1,
                          maxit = 1e3,
                          verbose = T,
                          debug = F,
                          version = 1
) {

  # Disabled options
  #
  anneal <- NULL

  # Anneal
  #
  if (!is.null(anneal)) {

    vec_c <- get_annealing_ladder_(anneal, verbose = verbose)

  } else {

    vec_c <- NULL

  }
  c <- 1

  # Save inputs
  #
  args <-
    list(
      Y = Y,
      list_hyper = list_hyper,
      list_init = list_init,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      debug = debug,
      version = version
    )

  if(is.null(version)){

    warning('No versions are specified. Set to 1 by default.')
    version <- 1

  }else{

    if (version!=1 & version!=2)
      stop('version takes arguments 1 (a beta prior on edge inclusion) or 2 (a normal prior on probit edge inclusion)')

  }

  # Time
  #
  pt <- Sys.time()

  # Dimension
  #
  P <- ncol(Y)
  N <- nrow(Y)

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

  if(version == 1){

    list_hyper <- set_default(list_hyper, 'a_rho', 2)
    if(list_hyper$a_rho <= 0)stop("a_rho must be positive.")

    list_hyper <- set_default(list_hyper, 'b_rho', 2)
    if(list_hyper$b_rho <= 0)stop("b_rho must be positive.")

  }else if (version == 2){

    list_hyper <- set_default(list_hyper, 'n0', -2)

    list_hyper <- set_default(list_hyper, 't02', 0.5)
    if(list_hyper$t02 <= 0)stop("t02 must be positive.")

  }

  # list2env(list_hyper, envir = .GlobalEnv)
  lambda <- list_hyper$lambda
  v0 <- list_hyper$v0
  v1 <- list_hyper$v1
  a_tau <- list_hyper$a_tau
  b_tau <- list_hyper$b_tau
  if(version == 1){
    a_rho <- list_hyper$a_rho
    b_rho <- list_hyper$b_rho
  }else if(version == 2){
    n0 <- list_hyper$n0
    t02 <- list_hyper$t02
  }
  #

  #
  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Preparing the parameter initialisation ... \n\n")

  # Specify initial parameters unless provided
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

  if(version == 1){

    if (!'alpha_rho' %in% names(list_init) | ('alpha_rho' %in% names(list_init) & is.null(list_init$alpha_rho))) {

      list_init$alpha_rho <- 1

    }else{

      if( list_init$alpha_rho <= 0){
        stop("alpha_rho must be positive.")
      }

    }

    if (!'beta_rho' %in% names(list_init) | ('beta_rho' %in% names(list_init) & is.null(list_init$beta_rho))) {

      list_init$beta_rho <- 1

    }else{

      if( list_init$beta_rho <= 0){
        stop("beta_rho must be positive.")
      }

    }

  }else if(version == 2){

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


  #
  Omega <- list_init$Omega
  alpha_tau <- list_init$alpha_tau
  beta_tau <- list_init$beta_tau

  if(version == 1 ){

    alpha_rho <- list_init$alpha_rho
    beta_rho <- list_init$beta_rho

  }else if(version == 2){

    mu_zeta <- list_init$mu_zeta
    sig2_inv_zeta <- list_init$sig2_inv_zeta

  }
  #


  # Initialise deduced quantities
  #
  m_tau <- get_m_tau(alpha_tau, beta_tau)

  if(version==1){

    m_log_rho <- get_m_log_o(alpha_rho, beta_rho)
    m_log_one_minus_rho <- get_m_log_one_minus_o(alpha_rho, beta_rho)

  }else if(version==2){

    theta <- rep(0, P)
    m1_alpha <- get_Alpha(theta, mu_zeta, P)

  }


  if (verbose) cat("... done. == \n\n")

  # debug mode
  #
  if (debug) {

    # track ELBO after each maximisation step
    vec_ELBO_VBEM <- c()

    # track ELBO within each variational step
    list_ELBO <- list()

    # record number of warnings
    n_warning <- 0
    vec_n_warning_VB <- c()

  }

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
    # ====== #

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
      if(version == 1){

        m_delta <- update_m_delta_betap(Omega, m_tau, m_log_rho, m_log_one_minus_rho, v0, v1, vbc)

      }else if(version == 2){

        m_delta <- update_m_delta(Omega, m_tau, m1_alpha, v0, v1, vbc)

      }

      E1 <- get_E1(m_delta, v0, v1)
      # % #


      # % # tau
      alpha_tau <- update_alpha_tau(a_tau, P, vbc)
      beta_tau <- update_beta_tau(Omega, E1, b_tau, vbc)

      m_tau <- get_m_tau(alpha_tau, beta_tau)
      m_log_tau <- get_m_log_tau(alpha_tau, beta_tau)
      # % #


      if(version == 1){

        # % # rho
        alpha_rho <- update_alpha_rho(m_delta, a_rho, vbc)
        beta_rho <- update_beta_rho(m_delta, b_rho, vbc)
        m_log_rho <- get_m_log_o(alpha_rho, beta_rho)
        m_log_one_minus_rho <- get_m_log_one_minus_o(alpha_rho, beta_rho)
        # % #

      }else if(version ==2){

        # % # z
        m1_z <- get_E2(m_delta, m1_alpha, vbc)
        m2_z <- get_E2_2(m1_z, m1_alpha, vbc)
        # % #

        # % # zeta
        sig2_inv_zeta <- update_sig2_inv_zeta(t02, P, vbc)
        mu_zeta <-
          update_mu_zeta_gm(sig2_inv_zeta, m1_z, n0, t02, vbc)
        m2_zeta <- get_m2_zeta(mu_zeta, sig2_inv_zeta)
        # % #

        # % #
        m1_alpha <- get_Alpha(theta, mu_zeta, P)
        sum_var_alpha <- P * (P-1) / 2 * sig2_inv_zeta^(-1)
        # % #

      }



      # % #
      if(version == 1){
        ELBO <-  get_elbo_gm_vbem_v1(Omega,
                                     m_delta,
                                     alpha_tau,
                                     beta_tau,
                                     m_tau,
                                     m_log_tau,
                                     alpha_rho,
                                     beta_rho,
                                     m_log_rho,
                                     m_log_one_minus_rho,
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
                                     P,
                                     vbc)

      }else if(version == 2){

        ELBO <-  get_elbo_gm_vbem_v2(Omega,
                                     m_delta,
                                     alpha_tau,
                                     beta_tau,
                                     m_tau,
                                     m_log_tau,
                                     mu_zeta,
                                     sig2_inv_zeta,
                                     m2_zeta,
                                     m1_alpha,
                                     sum_var_alpha,
                                     E1,
                                     S,
                                     lambda,
                                     v0,
                                     v1,
                                     n0,
                                     t02,
                                     a_tau,
                                     b_tau,
                                     N,
                                     P,
                                     vbc)

      }

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

      if(version==1){

        ELBO_M0 <-  get_elbo_gm_vbem_v1(Omega,
                                        m_delta,
                                        alpha_tau,
                                        beta_tau,
                                        m_tau,
                                        m_log_tau,
                                        alpha_rho,
                                        beta_rho,
                                        m_log_rho,
                                        m_log_one_minus_rho,
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
                                        P,
                                        vbc)

      }else if(version == 2){

        ELBO_M0 <- get_elbo_gm_vbem_v2(Omega,
                                       m_delta,
                                       alpha_tau,
                                       beta_tau,
                                       m_tau,
                                       m_log_tau,
                                       mu_zeta,
                                       sig2_inv_zeta,
                                       m2_zeta,
                                       m1_alpha,
                                       sum_var_alpha,
                                       E1,
                                       S,
                                       lambda,
                                       v0,
                                       v1,
                                       n0,
                                       t02,
                                       a_tau,
                                       b_tau,
                                       N,
                                       P,
                                       vbc)

      }
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
    if(version == 1){

      ELBO_M <-  get_elbo_gm_vbem_v1(Omega,
                                     m_delta,
                                     alpha_tau,
                                     beta_tau,
                                     m_tau,
                                     m_log_tau,
                                     alpha_rho,
                                     beta_rho,
                                     m_log_rho,
                                     m_log_one_minus_rho,
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
                                     P,
                                     vbc)

    }else if(version == 2){

      ELBO_M <- get_elbo_gm_vbem_v2(Omega,
                                    m_delta,
                                    alpha_tau,
                                    beta_tau,
                                    m_tau,
                                    m_log_tau,
                                    mu_zeta,
                                    sig2_inv_zeta,
                                    m2_zeta,
                                    m1_alpha,
                                    sum_var_alpha,
                                    E1,
                                    S,
                                    lambda,
                                    v0,
                                    v1,
                                    n0,
                                    t02,
                                    a_tau,
                                    b_tau,
                                    N,
                                    P,
                                    vbc)

    }


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

    if (debug) {
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

  if (version ==1){

    estimates <- list( Omega = Omega,
                       m_delta = m_delta,
                       alpha_tau = alpha_tau,
                       beta_tau = beta_tau,
                       alpha_rho = alpha_rho,
                       beta_rho =  beta_rho,
                       S = S # for model comparison
    )

  }else if(version == 2){

    estimates <- list( Omega = Omega,
                       m_delta = m_delta,
                       alpha_tau = alpha_tau,
                       beta_tau = beta_tau,
                       mu_zeta = mu_zeta,
                       sig2_inv_zeta = sig2_inv_zeta,
                       S = S # for model comparison
    )

  }


  if(debug){

    debugs <- list( n_warning = n_warning,
                    vec_n_warning_VB = vec_n_warning_VB,
                    list_ELBO = list_ELBO,
                    vec_ELBO_VBEM = vec_ELBO_VBEM)

  }else{

    debugs <- NA

  }

  create_named_list_(
    args,
    estimates,
    debugs,
    it,
    vec_VB_it,
    pt
  )
}


