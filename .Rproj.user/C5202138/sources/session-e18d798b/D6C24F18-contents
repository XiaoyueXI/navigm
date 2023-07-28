# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm
#
# Internal functions to call the EM algorithm for GMSS



gmss_em_core <-  function(Y,
                          V,
                          list_hyper = NULL,
                          list_init = NULL,
                          tol = 1e-1,
                          maxit = 1e3,
                          verbose = T,
                          track_ELBO = F,
                          debug = F) {


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
  N <- nrow(Y)
  P <- ncol(Y)
  Q <- ncol(V)

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

  list_hyper <- set_default(list_hyper, 's0', 0.1)
  if(list_hyper$s0 <= 0)stop("s0 must be positive.")

  list_hyper <- set_default(list_hyper, 's1', 100)
  if(list_hyper$s1 <= 0)stop("s1 must be positive.")
  if(list_hyper$s0 >= list_hyper$s1)stop("s1 should be much greater than s0.")


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
  list2env(list_hyper, envir = .GlobalEnv)
  if(!is.null(list_init))list2env(list_init, envir = .GlobalEnv)

  if (!'Omega' %in% names(list_init) | ('Omega' %in% names(list_init) & is.null(list_init$Omega))) {

    Sigma <- (S + v0 * diag(P)) / N
    Omega <- as.matrix(Matrix::nearPD(solve(Sigma))$mat)

  }else{

    if(nrow(Omega)!=ncol(Omega)){
      stop("Omega should be initialised as a square matrix")
    }else if(matrixcalc::is.symmetric.matrix(Omega)){
      stop("Omega should be initialised as a symmetric matrix")
    }else if(matrixcalc::is.positive.definite(Omega)){
      stop("Omega should be initialised as a positive definite matrix")
    }

  }

  if (!'beta' %in% names(list_init) | ('beta' %in% names(list_init) & is.null(list_init$beta))) {

    beta <- rep(0, Q)

  }else{

    if(length(beta) != ncol(V)){
      stop("Length of beta does not match the number of columns in V.")
    }

  }

  if (!'zeta' %in% names(list_init) | ('zeta' %in% names(list_init) & is.null(list_init$zeta))) {

    zeta <- list_hyper$n0

  }

  if (!'tau1' %in% names(list_init) | ('tau1' %in% names(list_init) & is.null(list_init$tau1))) {

    tau1 <- 1

  }else{

    if( tau1 <= 0){
      stop("tau1 must be positive.")
    }

  }

  if (!'tau2' %in% names(list_init) | ('tau2' %in% names(list_init) & is.null(list_init$tau2))) {

    tau2 <- 1

  }else{

    if( tau2 <= 0){
      stop("tau2 must be positive.")
    }

  }

  if (!'o' %in% names(list_init) | ('o' %in% names(list_init) & is.null(list_init$o))) {

    o <- 1/Q

  }else{

    if( o < 0 | o > 1){
      stop("o must be in [0,1].")
    }

  }

  if (verbose) cat("... done. == \n\n")

  # track ELBO
  #
  if (track_ELBO) {

    vec_ELBO_M <- c()

  } else{

    vec_ELBO_M <- NA

  }

  # debug mode
  #
  if (debug) {

    # record number of warnings
    n_warning <- 0

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
  ELBO_diff <- Inf
  ELBO_old <- -Inf

  while (ELBO_diff > tol & it < maxit) {

    # Iteration
    #
    it <- it + 1

    if (verbose & it %% 5 == 0)
      cat(paste0("Iteration ", format(it), "... \n"))

    # Anneal
    #
    if (!is.null(vec_c) && it <= length(vec_c)) {
      c <- vec_c[it]
      print(paste0("Temperature: ", format(1 / c, digits = 3)))
    } else {
      c <- 1
    }

    # E step :
    # ====== #

    theta <- get_theta(beta, V)
    Alpha <- get_Alpha(theta, zeta, P)
    P1 <- get_P1(Omega, Alpha, tau1, v0, v1, c)
    P2 <- get_P2(beta, o, tau2, s0, s1, c)
    E1 <- get_E1(P1, v0, v1)
    E2 <- get_E2(P1, Alpha, c)
    E2_2 <- get_E2_2(E2, Alpha, c)
    E5 <- get_E5(P2, s0, s1)

    # M step :
    # ====== #

    if (isTRUE(all.equal(c, 1))) {
      if (debug == T) {
        ELBO0 <- get_elbo_gmss_em(Omega,
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
                                  Q)
      }
    }

    # save order as VBEM
    #
    o <- get_o(P2, a_o, b_o)
    tau1 <- get_tau1(Omega, E1, a_tau, b_tau, P)
    tau2 <- get_tau2(beta, E5, a_sigma, b_sigma, Q)
    zeta <- get_zeta(E2, theta, n0, t02, P)
    beta <- get_beta(beta, zeta, tau2, E2, E5, sV, pV, spV1, spV2, P, Q)


    # omega
    #
    bool_cpp <- F
    bool_direct_solve <- F # keep F as otherwise Omega inverted at each iteration

    if (bool_cpp) {
      if (bool_direct_solve) {
        out <- M_Omega_direct_solve(N, P, Omega, S, lambda, tau1 * E1)
        Omega <- out$Omega
      } else {
        out <- M_Omega(N, P, Sigma, Omega, S, lambda, tau1 * E1)
        Omega <- out$Omega
        Sigma <- out$Sigma
      }
    } else {
      Omega <- get_omega(tau1 * E1, S, Omega, lambda, N, P)
    }


    #
    if (isTRUE(all.equal(c, 1))) {
      ELBO <- get_elbo_gmss_em(Omega,
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
                               Q)

      ELBO_diff <- abs(ELBO - ELBO_old)

      if (debug == T && ELBO + eps < ELBO0) {
        # not a check for the whole iteration as the M-step will increase ELBO, but the E-step may decrease it
        #
        warning(paste0("Non-increasing in the M-step : ELBO0 = ", ELBO0, ", ELBO = ", ELBO))
        n_warning <- n_warning + 1
      }

      if (verbose & it %% 5 == 0)
        cat(paste0(
          "Difference ELBO from previous iteration: ",
          format(ELBO_diff),
          "\n"
        ))

      ELBO_old <- ELBO

      if (track_ELBO) {
        vec_ELBO_M <- c(vec_ELBO_M, ELBO)
      }
    }
  }

  if(ELBO_diff <= tol){
    if(verbose)
      cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                 "Optimal objective function in the EM ",
                 "(Q) = ", format(ELBO), ". \n\n"))
  }

  if (it == maxit) {
    warning('maxiter is reached. The algorithm may not converge ... \n')
  }

  pt <- Sys.time() - pt
  cat('Algorithm runtime: ',format(pt), '\n')


  estimates <- list(Omega = Omega,
                    zeta = zeta,
                    beta = beta,
                    o = o,
                    tau1 = tau1,
                    tau2 = tau2,
                    P1 = P1,
                    P2 = P2,
                    E1 = E1,
                    E2 = E2,
                    E2_2 = E2_2,
                    E5 = E5,
                    S = S # for model comparison
  )


  if(debug){
    debugs <- list( n_warning = n_warning )
  }else{
    debugs <- NA
  }


  create_named_list_(
    args,
    estimates,
    debugs,
    it,
    vec_ELBO_M,
    pt
  )

}
