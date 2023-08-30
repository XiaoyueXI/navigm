# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm


#' Node-level auxiliary variable for improved Gaussian graphical spike-and-slab model.
#'
#' This function conducts simultaneous inference of a Gaussian graphical spike-and-slab model and
#' the effects of node-level auxiliary variables.
#' The function fixes the spike-and-slab configuration, except that its scale is estimated.
#'
#' @param Y Data matrix of dimensions N x P, where N is the number of samples and P is the number of nodes in the graph.
#' \code{Y} will be centered within the \code{navigm_core} call.
#'
#' @param V Input matrix of dimensions P x Q, where Q is the number of candidate auxiliary variables.
#' There's no need to supply an intercept.
#' If \code{V = NULL} or is not specified, the default approach "GMSS" cannot be employed (see the argument \code{method}) and thus set \code{method = 'GM'}.
#'
#' @param method Character: The method to use, which includes
#' "GM" (vanilla spike-and-slab graphical model),
#' "GMN" (spike-and-slab graphical model with normal priors on the node-level auxiliary variable coefficients), and
#' "GMSS" (default; spike-and-slab graphical model with spike-and-slab priors on the node-level auxiliary variable coefficients).
#'
#' @param inference Character: The inference algorithm to use, which includes
#' "ECM" (expectation conditional maximisation algorithm),
#' "VBECM" (default; variational Bayes expectation conditional maximisation algorithm).
#'
#' @param list_hyper A list containing the model hyperparameters:
#'   Always specify \code{lambda}, \code{v0}, \code{v1}, \code{a_tau}, and \code{b_tau} in all the models.
#'   In GM (version 1), you need additional \code{a_rho} and \code{b_rho}.
#'   In GM (version 2), you need additional \code{n0} and \code{t02}.
#'   In GMN, you need additional \code{n0}, \code{t02}, \code{a_sigma}, and \code{b_sigma}.
#'   In GMSS-VBECM, you need additional \code{n0}, \code{t02}, \code{a_sigma}, and \code{b_sigma}, \code{a_o}, and \code{b_o}.
#'   In GMSS-ECM, you need additional \code{s0_v}, \code{s1}, \code{n0}, \code{t02}, \code{a_sigma}, \code{b_sigma}, \code{a_o}, and \code{b_o}.
#'   Parameter interpretations:
#'   (1) \code{lambda}: a rate parameter of exponential priors on diagonal entries.
#'   (2) \code{v0} and \code{v1}: unscaled spike and slab variances in the bottom-level continuous spike-and-slab.
#'   (3) \code{a_tau} and \code{b_tau}: shape and rate parameters of a gamma prior on the bottom-level continuous spike-and-slab scale.
#'   (4) \code{a_rho} and \code{b_rho}: shape parameters of a beta prior on the probabilities of including edges.
#'   (5) \code{n0} and \code{t02}: mean and variance of a normal prior on the overall sparsity.
#'   (6) \code{a_sigma} and  \code{b_sigma}: shape parameters of a gamma prior on the scale of the top-level continuous spike-and-slab (GMSS-ECM),
#'   or on the top-level slab precision (GMSS-VBECM), or on the regression coefficients' precision (GMN).
#'   (7) \code{a_o} and \code{b_o}: shape parameters of a beta prior on the probabilities of including node-level variables.
#'   (8) \code{s0} and \code{s1}: unscaled spike and slab variances in the top-level continuous spike-and-slab.
#'
#'   If NULL or any hyperparameter specification is missing,
#'   they will be set to the default: \code{lambda = 2, v0 = s0 = 0.1, v1 = s1 = 100},
#'   \code{a_tau = b_tau = a_sigma = b_sigma = 2},
#'   \code{n0 = -2, t02 = 0.5, a_o = 1, b_o = Q, a_rho = 1, b_rho = P}.
#'
#' @param list_init A list containing the initialisations:
#'  Specify \code{Omega} in all the models. In ECM, specify \code{tau_1}. In VBECM, specify \code{alpha_tau} and \code{beta_tau}.
#'  In GM (version 1), you need \code{alpha_rho} and \code{beta_rho} using VBECM, and \code{rho} using ECM.
#'  In GM (version 2), you need \code{mu_zeta} and \code{sig2_inv_zeta} using VBECM, and \code{zeta} using ECM.
#'  In GMN, you need additional \code{mu_zeta}, \code{sig2_inv_zeta}, \code{mu_beta}, \code{sig2_inv_beta}, \code{alpha_sigma}, and \code{beta_sigma} using VBECM,
#'  and \code{zeta}, \code{beta}, and \code{tau2} using ECM.
#'  In GMSS, you need additional \code{mu_zeta}, \code{sig2_inv_zeta}, \code{mu_beta}, \code{sig2_inv_beta},
#'  \code{alpha_o}, \code{beta_o}, \code{alpha_sigma} and \code{beta_sigma} using VBECM,
#'  and \code{zeta}, \code{beta}, \code{o}, and \code{tau2} using ECM.
#'
#'   Parameter interpretations:
#'   (1) \code{Omega}: a P x P precision matrix.
#'   VBECM:
#'   (2) \code{alpha_tau} and \code{beta_tau}: shape and rate parameters of a gamma variational distribution on the bottom-level continuous spike-and-slab scale.
#'   (3) \code{alpha_rho} and \code{beta_rho}: shape parameters of a beta variational distribution on edge inclusion probabilities.
#'   (4) \code{mu_zeta} and \code{sig2_inv_zeta}: mean and inverse variance of a normal variational distribution on the overall sparsity.
#'   (5) \code{mu_beta} and \code{sig2_inv_beta}: mean and inverse variance of a normal variational distribution on regression coefficients.
#'   (6) \code{alpha_sigma} and \code{beta_sigma}: shape and rate parameters of a gamma variational distribution on the top-level continuous spike-and-slab scale (GMSS-ECM),
#'   or on the top-level slab precision (GMSS-VBECM), or on the regression coefficients' precision (GMN).
#'   (7) \code{alpha_o} and \code{beta_o}: shape parameters of a beta variational distribution on node-level variable inclusion probabilities.
#'   Parameter interpretations (ECM):
#'   (2) \code{tau_1}: the bottom-level continuous spike-and-slab scale.
#'   (3) \code{rho}: an edge variable inclusion probability.
#'   (4) \code{zeta}: the overall sparsity.
#'   (5) \code{beta}: regression coefficients.
#'   (6) \code{tau_2}: the top-level continuous spike-and-slab scale (GMSS-ECM), or the top-level slab precision (GMSS-VBECM), or the regression coefficients' precision (GMN).
#'   (7) \code{o}: a node-level variable inclusion probability.
#'
#'  If \code{NULL} or any initialization specification is missing,
#'  they will be set to the default:
#'  \code{Omega = } the empirical precision matrix.
#'  VBECM:
#'  \code{mu_beta = rep(0, Q)}, \code{sig2_inv_beta = rep(1, Q)},
#'  \code{alpha_tau = beta_tau = alpha_sigma = beta_sigma = 1}, \code{mu_zeta = list_hyper$n0},
#'  \code{sig2_inv_zeta = 1/list_hyper$t02}, \code{a_o = 1, b_o = Q},  \code{a_rho = 1, b_rho = P},
#'  ECM:
#'  \code{beta = rep(0, Q)}, \code{tau_1 = tau_2 = 1}, \code{zeta = list_hyper$n0}, \code{o = 1 / Q}, \code{rho = 1 / P}.
#'
#' @param ne0 A vector of length 2 containing the prior expectation and variance of
#' the number of edges in the absence of hubs. This will not be used if the \code{n0} in \code{list_hyper} is non-\code{NULL}.
#'
#' @param tol Scalar: tolerance for the stopping criterion (default is 0.001).
#'
#' @param maxit Scalar: maximum number of allowed iterations (default is 100,000).
#'
#' @param transformV Logical: if \code{FALSE} (default), \code{V} will not be transformed;
#' Otherwise, if \code{V} does not range within [0,1], \code{V} will be standardised within the \code{navigm} call.
#'
#' @param verbose Logical: if \code{TRUE} (default), standard verbosity; otherwise, no messages.
#'
#' @param debug Logical: if \code{FALSE} (default), additional terms will not be recorded for debugging purposes;
#' otherwise, the number of decreasing ELBOs in the maximisation step (ECM, VBECM) and within each variational update (VBECM) will be recorded,
#' possibly due to numerical round-off. The ELBOs will be tracked after the maximisation step (ECM, VBECM) and within each variational update (VBECM).
#'
#' @param version Integer: it takes values of 1 (indicating a beta prior on edge inclusion probabilities)
#' or 2 (indicating a normal prior on probit edge inclusion probabilities).
#' This option is only valid when \code{method = 'GM'}).
#'
#'
#' @details The \code{navigm} function includes a core Gaussian graphical spike-and-slab model
#' that enables the incorporation and selection of node-level auxiliary variables,
#' thereby enhancing the detection of conditional dependence.
#' This function utilises a user-specified spike-and-slab configuration
#' and does not permit selection from multiple alternatives.
#' Inference is conducted using a scalable (variational) expectation maximisation
#' algorithm, making it suitable for high-dimensional graphs.
#'
#'
#'
#' @return A list containing the following quantities:
#'  \describe{
#'  \item{estimates}{A list containing model estimates:
#'  \code{Omega} in all the models.
#'  In VBECM, \code{a_tau, b_tau}, \code{m_delta}. In ECM, \code{tau_1}, \code{P1}.
#'  In GM (version 1), \code{a_rho, b_rho} using VBECM, \code{rho} using ECM.
#'  In GM (version 2), \code{mu_zeta, sig2_inv_zeta} using VBECM, \code{zeta} using ECM.
#'  In GMN, \code{mu_zeta, sig2_inv_zeta}, \code{mu_beta, sig2_inv_beta}, \code{a_sigma, b_sigma} using VBECM. \code{zeta}, \code{beta}, \code{tau2} using ECM.
#'  In GMSS, \code{mu_zeta, sig2_inv_zeta}, \code{mu_beta, sig2_inv_beta}, \code{a_o, b_o}, \code{a_sigma, b_sigma}using VBECM. \code{zeta}, \code{beta}, \code{o}, \code{tau2} using ECM.
#'  All the parameter interpretations are in Arguments section \code{list_init}, except that
#'  \code{m_delta, P1} refer to a P x P matrix containing posterior probabilities of including edges, and
#'  \code{m_gamma, P2} refer to Q posterior probabilities of including auxiliary variables.
#'
#' }
#'  \item{debugs}{A list containing terms for the purpose of debugging:
#'  In VBECM, return \code{n_warning}, \code{vec_n_warning_VB}, \code{list_ELBO}, \code{vec_ELBO_CM}. In ECM, return \code{n_warning}, \code{vec_ELBO_CM}.
#'  The meaning of these quantities:
#'  (1) \code{n_warning}: Scalar. Number of ELBO drops in the maximisation step.
#'  (2) \code{vec_n_warning_VB}: Vector of length \code{it}. Number of ELBO drops within each variational update.
#'  (3) \code{vec_ELBO_CM}: Vector of length \code{it}. Track ELBO after the maximisation step.
#'  (4) \code{list_ELBO}: List of length \code{it}. Track ELBO within each variational update.  }
#' \item{it}{Scalar: total number of iterations. }
#' \item{args}{A list containing the input arguments.}
#' \item{vec_VB_it}{A vector of length \code{it} containing the number of iterations within each variational step}
#' \item{pt}{Scalar: algorithm runtime in seconds. }
#' }
#'
#' @examples
#' # simulate data
#'
#' seed <- 123; set.seed(seed)
#' N <- 100; P <- 5; Q <- 3
#' V <- matrix(rnorm(P*Q), nrow = P, ncol = Q)
#' Y <- matrix(rnorm(N*P), nrow = N, ncol = P)
#'
#' # estimate precision matrix based on Y and meanwhile leverage node-level variables V
#' res_navigm_core <-navigm_core(Y, V)
#' # res_navigm_core <- navigm_core(Y, V, method = 'GMN')
#' # res_navigm_core <- navigm_core(Y, V, method = 'GM', version = 1)
#' # res_navigm_core <- navigm_core(Y, V, method = 'GM', version = 2)
#' # res_navigm_core <- navigm_core(Y, V, inference = 'ECM')
#' # res_navigm_core <- navigm_core(Y, V, method = 'GMN', inference = 'ECM')
#' # res_navigm_core <- navigm_core(Y, V, method = 'GM', inference = 'ECM', version = 1)
#' # res_navigm_core <- navigm_core(Y, V, method = 'GM', inference = 'ECM', version = 2)
#'
#' @export
#'


navigm_core <- function(Y, V =NULL,
                        method = 'GMSS',
                        inference = 'VBECM',
                        list_hyper = NULL,
                        list_init = NULL,
                        ne0 = NULL,
                        tol = 0.1,
                        maxit = 1000,
                        # transformY = T,# as mean is 0, always centre Y
                        transformV = F,
                        debug = F,
                        verbose = T,
                        version = NULL) {


  if (verbose){

    cat(paste0("\n======================= \n",
               "== Preprocess data ... == \n",
               "======================= \n\n"))

  }


  # preprocess not in navigm
  # such that navigm_core is usable standalone.
  #
  if(!is.null(V)){

    # min requirement to include node-level variables
    #
    if(ncol(Y)!=nrow(V)){
      stop('Columns of Y and rows of V must match.')
    }

    # if scale V
    #
    if(transformV){
      # if V \in [0,1] not scale
      if(min(V) < 0 | max(V) > 1){
        V <- scale(V, center = TRUE, scale = TRUE)
      }else{
        'V in [0,1] and thus not transformed.'
      }
    }

    # no node-level variables
    #
    if(!is.null(V) & method == 'GM' ){
      warning("node-level auxiliary variables (V) will not be used in GM. Please consider GMN and GMSS to leverage V.")
    }

    # low number of node-level variables
    #
    if(ncol(V) <= 10 & method =='GMSS'){
      warning("Number of node-level auxiliary variables <= 10, not advised to select and change to GMN.")
    }

    # high number of node-level variables
    #
    if(ncol(V) > 10 & method =='GMN'){
      warning("Number of node-level auxiliary variables > 10, advised to select and change to GMSS.")
    }

    #
    if(is.null(V) & (method == 'GMN' | method == 'GMSS')){
      warning('node-level auxiliary variables (V) must be provided in GMN and GMSS. Change to GM.')
      method <- 'GM'
    }

    # variable names
    #
    if (is.null(colnames(Y)) & is.null(rownames(V)))
      colnames(Y) <- rownames(V) <- paste0("Node_", 1:ncol(Y))
    else if (is.null(colnames(Y))) colnames(Y) <- rownames(V)
    else if (is.null(rownames(V))) rownames(V) <- colnames(Y)
    else if (any(colnames(Y) != rownames(V)))
      stop("The provided column names of Y and row names of V must be the same.")

    if (is.null(rownames(Y))) rownames(Y) <- paste0("Ind_", 1:nrow(Y))
    if (is.null(colnames(V))) colnames(V) <- paste0("Var_", 1:ncol(V))

  }


  # as mean is 0, always centre Y
  #
  # if(transformY){
  Y <- scale(Y, center = TRUE, scale = FALSE) # center the data
  # }


  # version
  #
  if(!is.null(version) & (method == 'GMN' | method == 'GMSS')){
    warning("version is not used in GMN & GMSS.")
  }

  # ne0
  #
  if(!is.null(list_hyper) & 'n0' %in% names(list_hyper) & 't02' %in% names(list_hyper)){

    warning(paste0("Provided argument ne0 not used, as n0 and t02 were provided in list_hyper."))

  }else if(!is.null(ne0)){

    if(length(ne0)!=2){

      warning('ne0 should have length 2, specifying expected number of edges and variance respectively. \n
              As no valid list_hyper and ne0 are provided, use default hyperparameters n0 = -2 and t02 = 0.5.')

    }else{
      P <- ncol(Y)
      tmp <- get_n0_t02(P * (P-1) / 2, ne0)
      list_hyper$n0 <-  tmp$n0
      list_hyper$t02 <-  tmp$t02

    }
  }else{

    warning('As no valid list_hyper and ne0 are provided, use default hyperparameters n0 = -2 and t02 = 0.5.')

  }


  if (verbose) cat("... done. == \n\n")


  cat("**************************************************** \n",
      "Number of samples: ", nrow(Y), "\n",
      "Number of nodes: ", ncol(Y), "\n")


  if(method == 'GMN' | method == 'GMSS'){

    cat("Number of candidate node-level auxiliary variables: ", ncol(V), "\n",
        "**************************************************** \n\n")

  }

  if (verbose){
    cat(paste0("\n=========================================================\n",
               "== Perform inference using the ",inference, " algorithm ... == \n",
               "======================================================== \n\n"))
  }


  if(inference == 'ECM'){

    res_navigm <- navigm_ecm_core(Y = Y,
                                  V = V,
                                  method = method,
                                  list_hyper = list_hyper,
                                  list_init = list_init,
                                  tol = tol,
                                  maxit = maxit,
                                  verbose = verbose,
                                  debug = debug,
                                  version = version)

  }else if(inference == 'VBECM'){

    res_navigm <- navigm_vbecm_core(Y = Y,
                                    V = V,
                                    method = method,
                                    list_hyper = list_hyper,
                                    list_init = list_init,
                                    tol = tol,
                                    maxit = maxit,
                                    verbose = verbose,
                                    debug = debug,
                                    version = version)

  }

  return(res_navigm)

}


