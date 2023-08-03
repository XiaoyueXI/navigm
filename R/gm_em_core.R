# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm
#
# Internal functions to call the EM algorithm for GM


#' @importFrom matrixcalc is.symmetric.matrix is.positive.definite
#' @importFrom Matrix nearPD
gm_em_core <-  function(Y,
                        list_hyper = NULL,
                        list_init = NULL,
                        tol = 1e-1,
                        maxit = 1e3,
                        seed = 123,
                        verbose = T,
                        track_ELBO = F,
                        debug = F,
                        version = 1) {


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
      list_hyper = list_hyper,
      list_init = list_init,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      track_ELBO = track_ELBO,
      debug = debug,
      version = version
    )

  if(is.null(version)){

    warning('No versions are specified. Set to 1 by default.')
    version <- 1

  }else{

    if (version!=1 & version!=2)
      stop('version takes arguments 1 (a betag prior on edge inclusion) or 2 (a normal prior on probit edge inclusion)')

  }

  # Time
  #
  pt <- Sys.time()

  # Dimension
  #
  N <- nrow(Y)
  P <- ncol(Y)

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
  #
  lambda <- list_hyper$lambda
  v0 <- list_hyper$v0
  v1 <- list_hyper$v1
  if(version == 1){
    a_rho <- list_hyper$a_rho
    b_rho <- list_hyper$b_rho
  }else if(version == 2){
    n0 <- list_hyper$n0
    t02 <- list_hyper$t02
  }
  a_tau <- list_hyper$a_tau
  b_tau <- list_hyper$b_tau
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

    if (!'rho' %in% names(list_init) | ('rho' %in% names(list_init) & is.null(list_init$rho))) {

      list_init$rho <- 1 / P

    }else{

      if( list_init$rho < 0 | list_init$rho > 1){
        stop("rho must be in [0,1].")
      }

    }

  }else if(version == 2){

    if (!'zeta' %in% names(list_init) | ('zeta' %in% names(list_init) & is.null(list_init$zeta))) {

      list_init$zeta <- list_hyper$n0

    }
  }

  if (!'tau1' %in% names(list_init) | ('tau1' %in% names(list_init) & is.null(list_init$tau1))) {

    list_init$tau1 <- 1

  }else{

    if( list_init$tau1 <= 0){
      stop("tau1 must be positive.")
    }

  }


  #
  Omega <- list_init$Omega
  if(version == 1){
    rho <- list_init$rho
  }else if(version == 2){
    zeta <-  list_init$zeta
  }
  zeta <- list_init$zeta
  tau1 <- list_init$tau1
  #

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

    if(version == 1){
      P1 <- get_P1(Omega, rho, tau1, v0, v1, c)
    }else if (version == 2){
      theta <- rep(0, P)
      Alpha <- get_Alpha(theta, zeta, P)
      P1 <- get_P1(Omega, Alpha, tau1, v0, v1, c)
      E2 <- get_E2(P1, Alpha, c)
      E2_2 <- get_E2_2(E2, Alpha, c)
    }
    E1 <- get_E1(P1, v0, v1)


    # M step :
    # ====== #

    if (isTRUE(all.equal(c, 1))) {

      if (debug == T) {

        if(version == 1){
          ELBO0 <-  get_elbo_gm_em_v1(Omega,
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
                                      P)
        }else if(version == 2){

          ELBO0 <- get_elbo_gm_em_v2(Omega,
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
                                     n0,
                                     t02,
                                     a_tau,
                                     b_tau,
                                     N,
                                     P)
        }

      }
    }

    # save order as VBEM
    #
    tau1 <- get_tau1(Omega, E1, a_tau, b_tau, P)
    if(version == 1){
      rho <- get_rho(P1, a_rho, b_rho)
    }else if (version == 2){
      zeta <- get_zeta(E2, theta, n0, t02, P)
    }


    # omega
    #
    # bool_cpp <- F
    # bool_direct_solve <- F # keep F as otherwise Omega inverted at each iteration
    #
    # if (bool_cpp) {
    #   if (bool_direct_solve) {
    #     out <- M_Omega_direct_solve(N, P, Omega, S, lambda, tau1 * E1)
    #     Omega <- out$Omega
    #   } else {
    #     out <- M_Omega(N, P, Sigma, Omega, S, lambda, tau1 * E1)
    #     Omega <- out$Omega
    #     Sigma <- out$Sigma
    #   }
    # } else {
      Omega <- get_omega(tau1 * E1, S, Omega, lambda, N, P)
    # }


    #
    if (isTRUE(all.equal(c, 1))) {
      if(version == 1){
        ELBO <-  get_elbo_gm_em_v1(Omega,
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
                                   P)
      }else if(version == 2){
        ELBO <- get_elbo_gm_em_v2(Omega,
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
                                  n0,
                                  t02,
                                  a_tau,
                                  b_tau,
                                  N,
                                  P)
      }


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

  if (version ==1){

    estimates <- list(Omega = Omega,
                      rho = rho,
                      tau1 = tau1,
                      P1 = P1,
                      E1 = E1,
                      S =S # for model comparison
    )

  }else if(version == 2){

    estimates <- list(Omega = Omega,
                      zeta = zeta,
                      tau1 = tau1,
                      P1 = P1,
                      E1 = E1,
                      E2 = E2,
                      E2_2 = E2_2,
                      S = S # for model comparison
    )
  }


  if(debug){
    debugs <- list( n_warning )
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
