# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm


#' Node-level auxiliary variable for improved Gaussian spike-and-slab.
#'
#' This function performs simultaneous inference of Gaussian spike-and-slab and node-level auxiliary variable.
#'
#' @param Y Data matrix of dimension N x P, where N is the number of samples and P is the number of nodes in the graph.
#'   \code{Y} will be centred within \code{navigm_core} call.
#'
#' @param V Input matrix of dimension P x Q, where Q is the number of candidate auxiliary variables. No need to supply intercept.
#'   If \code{V} is not specified, set \code{method = 'GM'}.
#'
#' @param method Character: the method to use, including
#'  "GM" (vanilla spike-and-slab graphical model),
#'  "GMN" (spike-and-slab graphical model with normal priors on the node-level auxiliary variable coefficients), and
#'  "GMSS" (default; spike-and-slab graphical model with spike-and-slab priors on the node-level auxiliary variable coefficients).
#'
#' @param inference Character: the inference algorithm to use, including
#'  "EM" (expectation maximisation algorithm),
#'  "VBEM" (default; variational Bayes expectation maximisation algorithm),
#'
#' @param list_hyper A list containing the model hyperparameters:
#'   always specify \code{lambda}, \code{v0}, \code{v1}, \code{a_tau} and \code{b_tau} in all the models.
#'   In GM (version 1), need additional \code{a_rho} and \code{b_rho}.
#'   In GM (version 2), need additional \code{n0} and \code{t02}.
#'   In GMN, need additional \code{n0}, \code{t02}, \code{a_sigma} and \code{b_sigma}.
#'   In GMSS-VBEM, need additional \code{n0}, \code{t02}, \code{a_sigma} and \code{b_sigma}.
#'   In GMSS-EM, need additional \code{s0}, \code{s1}, \code{n0}, \code{t02}, \code{a_sigma} and \code{b_sigma}.
#'   Parameter interpretations: (1) \code{lambda}: a rate parameter of exponential priors on diagonal entries.
#'   (2) \code{v0} and \code{v1}: unscaled spike and slab variances in the bottom-level continuous spike-and-slab.
#'   (3) \code{a_tau} and \code{b_tau}: shape and rate parameters of a gamma prior on the continuous spike-and-slab unscaled prercision.
#'   (4) \code{a_rho} and \code{b_rho}: shape parameters of a beta prior on edge inclusion probabilities.
#'   (5) \code{n0} and \code{t02}: mean and variance of a normal prior on the overall sparsity.
#'   (6) \code{a_sigma} and  \code{b_sigma}: shape parameters of a gamma prior on (slab) precision.
#'   (7) \code{a_o} and \code{b_o}: shape parameters of a beta prior on node-level variable inclusion probabilities.
#'   (8) \code{s0} and \code{s1}: unscaled spike and slab variances in the top-level continuous spike-and-slab.
#'
#'   if NULL or any hyperparameter specification is missing,
#'   set to the default, \code{lambda = 2, v0 = 0.1, v1 = 100, a_tau = b_tau = a_sigma = b_sigma = 2},
#'   \code{n0 = -2, t02 = 0.5, a_o = 1, b_o = Q, a_rho = 1, b_rho = P}.
#'
#' @param list_init A list containing the initialisations:
#'  specify \code{Omega} in all the models. In EM, specify \code{tau_1}. In VBEM, specify \code{alpha_tau} and \code{beta_tau}.
#'  In GM (version 1), need \code{alpha_rho} and \code{beta_rho} using VBEM, and \code{rho} using EM.
#'  In GM (version 2), need \code{mu_zeta} and \code{sig2_inv_zeta} using VBEM, and \code{zeta} using EM.
#'  In GMN, need additional \code{mu_zeta}, \code{sig2_inv_zeta}, \code{mu_beta}, \code{sig2_inv_beta}, \code{alpha_sigma}, and \code{beta_sigma} using VBEM, and \code{zeta}, \code{beta} and \code{tau2} using EM.
#'  In GMSS, need additional \code{alpha_o}, \code{beta_o}, \code{mu_zeta}, \code{sig2_inv_zeta}, \code{mu_beta}, \code{sig2_inv_beta}, \code{alpha_sigma},
#'  and \code{beta_sigma} using VBEM, and \code{zeta}, \code{beta}, \code{o} and \code{tau2} using EM.
#'
#'   Parameter interpretations:
#'   (1) \code{alpha_tau} and \code{beta_tau}: shape and rate parameters of a gamma variational distribution on the continuous spike-and-slab prercision scale.
#'   (2) \code{alpha_rho} and \code{beta_rho}: shape parameters of a beta variational distribution on edge inclusion probabilities.
#'   (3) \code{mu_zeta} and \code{sig2_inv_zeta}: mean and inverse variance of a normal variational distribution on the overall sparsity.
#'   (4) \code{mu_beta} and \code{sig2_inv_beta}: mean and inverse variance of a normal variational distribution on the Q regression coefficients.
#'   (5) \code{alpha_sigma} and  \code{beta_sigma}: shape parameters of a gamma variational distribution on (slab) precision.
#'   (6) \code{alpha_o} and \code{beta_o}: shape parameters of a beta variational distribution on node-level variable inclusion probabilities.
#'   (7) \code{tau_1}: the bottom-level continuous spike-and-slab precision scale.
#'   (8) \code{zeta}: the overall sparsity.
#'   (9) \code{beta}: Q regression coefficients.
#'   (10) \code{tau_2}: the top-level continuous spike-and-slab precision scale.
#'   (11) \code{o}: node-level variable inclusion probability.
#'   (12) \code{Omega}: P x P precision matrix.

#'  if \code{NULL} or any initialisation specification is missing,
#'  set to default, \code{mu_beta = rep(0,Q)}, \code{sig2_inv_beta = rep(1,Q)},
#'  \code{alpha_tau = beta_tau = alpha_sigma = beta_sigma = 1}, \code{mu_zeta = list_hyper$n0},
#'  \code{sig2_inv_zeta = 1/list_hyper$t02}, \code{a_o = 1, b_o = Q},  \code{a_rho = 1, b_rho = P},
#'  \code{beta = rep(0,Q)}, \code{tau_1 = tau_2 = 1}, \code{zeta = list_hyper$n0}, \code{o = 1 / Q}, \code{rho = 1 / P}.
#'
#' @param ne0 Vector of length 2 whose entries are the prior expectation and variance of
#'  the number of edges in absence of hubs. Will not be used if \code{n0} in \code{list_hyper} is non-\code{NULL}.
#'
#' @param tol Scalar: tolerance for the stopping criterion (default is 0.1).
#'
#' @param maxit Scalar: maximum number of iterations allowed (default is 1000).
#'
#' @param transformV Logical; if \code{FALSE} (default), \code{V} will not be transformed;
#'  Otherwise and if \code{V} does not range within [0,1],  \code{V} will be standardised within \code{navigm} call.
#'
#' @param verbose Logical; if \code{TRUE} (default), standard verbosity; otherwise, no messages.
#'
#' @param debug Logical; if \code{FALSE} (default), not record additional terms for the purpose of debugging;
#'  otherwise, record the number of decreasing ELBOs in the maximisation step (EM, VBEM) and within each variational update (VBEM),
#'  possibly due to numerical round-off; track ELBOs after the maximisation step (EM, VBEM) and within each variational update (VBEM).
#'
#' @param version Integer; take values of 1 (a beta prior on edge inclusion probabilities) or 2 (a normal prior on probit edge inclusion probabilities).
#'  (a valid option when \code{method = 'GM'}).
#'
#' @details \code{navigm_core} implements a Gaussian graphical model
#'   that allows incorporating and selecting node-level auxiliary variables,
#'   thereby enhancing the detection of conditional dependence. Inference is
#'   carried out using a scalable (variational) expectation maximisation
#'   algorithm, which is applicable to high-dimension graphs.
#'
#'
#' @return A list containing the following quantities:
#'  \describe{
#'  \item{estimates}{A list containing model estimates:
#'  \code{Omega} in all the models.
#'  In VBEM, \code{a_tau, b_tau}, \code{m_delta}. In EM, \code{tau_1}, \code{P1}.
#'  In GM (version 1), \code{a_rho, b_rho} using VBEM, \code{rho} using EM.
#'  In GM (version 2), \code{mu_zeta, sig2_inv_zeta} using VBEM, \code{zeta} using EM.
#'  In GMN, \code{mu_zeta, sig2_inv_zeta}, \code{mu_beta, sig2_inv_beta}, \code{a_sigma, b_sigma} using VBEM. \code{zeta}, \code{beta}, \code{tau2} using EM.
#'  In GMSS, \code{mu_zeta, sig2_inv_zeta}, \code{mu_beta, sig2_inv_beta}, \code{a_o, b_o}, \code{a_sigma, b_sigma}using VBEM. \code{zeta}, \code{beta}, \code{o}, \code{tau2} using EM.
#'  All the parameter interpretations are in Arguments section \code{list_init}, except that
#'  \code{m_delta, P1} refer to a P x P matrix containing posterior probabilities of including edges, and
#'  \code{m_gamma, P2} refer to Q posterior probabilities of including auxiliary variables.
#'
#' }
#'  \item{debugs}{A list containing terms for the purpose of debugging:
#'  In VBEM, return \code{n_warning}, \code{vec_n_warning_VB}, \code{list_ELBO}, \code{vec_ELBO_M}. In EM, return \code{n_warning}, \code{vec_ELBO_M}.
#'  The meaning of these quantities:
#'  (1) \code{n_warning}: Scalar. Number of ELBO drops in the maximisation step.
#'  (2) \code{vec_n_warning_VB}: Vector of length \code{it}. Number of ELBO drops within each variational update.
#'  (3) \code{vec_ELBO_M}: Vector of length \code{it}. Track ELBO after the maximisation step.
#'  (4) \code{list_ELBO}: List of length \code{it}. Track ELBO within each variational update.  }
#' \item{it}{Scalar. Total number of iterations. }
#' \item{vec_VB_it}{Vector of length \code{it}. Number of iterations within each variational update.}
#' \item{pt}{Scalar. Algorithm runtime in seconds. }
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
#' # res_navigm_core <- navigm_core(Y, V, inference = 'EM')
#' # res_navigm_core <- navigm_core(Y, V, method = 'GMN', inference = 'EM')
#' # res_navigm_core <- navigm_core(Y, V, method = 'GM', inference = 'EM', version = 1)
#' # res_navigm_core <- navigm_core(Y, V, method = 'GM', inference = 'EM', version = 2)
#'
#' @export
#'


navigm_core <- function(Y, V =NULL,
                         method = 'GMSS',
                         inference = 'VBEM',
                         list_hyper = NULL, list_init = NULL,
                         ne0 = NULL,
                         tol = 0.1, maxit = 1000,
                         verbose = T,
                         debug = F,
                         version = NULL,
                         # transformY = T,# as mean is 0, always centre Y
                         transformV = F) {


  if (verbose){
    cat(paste0("\n======================= \n",
               "== Preprocess data ... == \n",
               "======================= \n\n"))
  }

  if(!is.null(V)){

    if(ncol(Y)!=nrow(V)){
      stop('Columns of Y and rows of V must match.')
    }

    if(transformV){
      # if V \in [0,1] not scale
      if(min(V) < 0 | max(V) > 1){
        V <- scale(V, center = TRUE, scale = TRUE)
      }else{
        'V in [0,1] and thus not transformed.'
      }
    }

    if(!is.null(V) & method == 'GM' ){
      warning("node-level auxiliary variables (V) will not be used in GM. Please consider GMN and GMSS to leverage V.")
    }

    if(ncol(V) <= 10 & method =='GMSS'){
      warning("Number of node-level auxiliary variables <= 10, not advised to select and change to GMN.")
    }

    if(ncol(V) > 10 & method =='GMN'){
      warning("Number of node-level auxiliary variables > 10, advised to select and change to GMSS.")
    }

    if(is.null(V) & (method == 'GMN' | method == 'GMSS')){
      warning('node-level auxiliary variables (V) must be provided in GMN and GMSS. Change to GM.')
      method <- 'GM'
    }

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
  # if(transformY){
  Y <- scale(Y, center = TRUE, scale = FALSE) # center the data
  # }


  if(!is.null(version) & (method == 'GMN' | method == 'GMSS')){
    warning("version is not used in GMN & GMSS.")
  }

  if(!is.null(list_hyper) & 'n0' %in% names(list_hyper) & 't02' %in% names(list_hyper)){

    warning(paste0("Provided argument ne0 not used, as n0 and t02 were provided in list_hyper."))

  }else if(!is.null(ne0)){

    if(length(ne0)!=2){

      warning('ne0 should have length 2, specifying expected number of edges and variance respectively. \n
              As no valid list_hyper and ne0 are provided, use default hyperparameters n0 = -2 and t02 = 0.5.')

    }else{
      P <- ncol(Y)
      tmp <- get_n0_t02(P*(P-1)/2, ne0)
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


  if(inference == 'EM'){

    res_navigm <- navigm_em_core(Y = Y,
                                   V = V,
                                   method = method,
                                   list_hyper = list_hyper,
                                   list_init = list_init,
                                   tol = tol,
                                   maxit = maxit,
                                   verbose = verbose,
                                   track_ELBO = track_ELBO,
                                   debug = debug,
                                   version = version)
  }else if(inference == 'VBEM'){
    res_navigm <- navigm_vbem_core(Y = Y,
                                     V = V,
                                     method = method,
                                     list_hyper = list_hyper,
                                     list_init = list_init,
                                     tol = tol,
                                     maxit = maxit,
                                     verbose = verbose,
                                     track_ELBO = track_ELBO,
                                     debug = debug,
                                     version = version)
  }

  return(res_navigm)
}


