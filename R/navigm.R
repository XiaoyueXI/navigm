# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm



#' Node-level auxiliary variable for improved Gaussian graphical spike-and-slab model with automatic selection of spike variances.
#'
#' This function conducts simultaneous inference of a Gaussian graphical spike-and-slab model and
#' the effects of node-level auxiliary variables.
#' The function also enables the automatic selection of spike variance from a grid of candidates.
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
#' @param criterion Character: The model selection criterion to use, which includes
#' "AIC" (default; Akaike information criterion),
#' "BIC" (Bayesian information criterion),
#' and "EBIC" (extended Bayesian information criterion).
#'
#' @param list_hyper A list containing the model hyperparameters:
#'   Always specify \code{lambda}, \code{v0_v}, \code{v1}, \code{a_tau}, and \code{b_tau} in all the models.
#'   In GM (version 1), you need additional \code{a_rho} and \code{b_rho}.
#'   In GM (version 2), you need additional \code{n0} and \code{t02}.
#'   In GMN, you need additional \code{n0}, \code{t02}, \code{a_sigma}, and \code{b_sigma}.
#'   In GMSS-VBECM, you need additional \code{n0}, \code{t02}, \code{a_sigma}, and \code{b_sigma}, \code{a_o}, and \code{b_o}.
#'   In GMSS-ECM, you need additional \code{s0_v}, \code{s1}, \code{n0}, \code{t02}, \code{a_sigma}, \code{b_sigma}, \code{a_o}, and \code{b_o}.
#'   Parameter interpretations:
#'   (1) \code{lambda}: a rate parameter of exponential priors on diagonal entries.
#'   (2) \code{v0_v} and \code{v1}: a grid of unscaled spike variances and an unscaled slab variance in the bottom-level continuous spike-and-slab.
#'   (3) \code{a_tau} and \code{b_tau}: shape and rate parameters of a gamma prior on the bottom-level continuous spike-and-slab scale.
#'   (4) \code{a_rho} and \code{b_rho}: shape parameters of a beta prior on the probabilities of including edges.
#'   (5) \code{n0} and \code{t02}: mean and variance of a normal prior on the overall sparsity.
#'   (6) \code{a_sigma} and  \code{b_sigma}: shape parameters of a gamma prior on the scale of the top-level continuous spike-and-slab (GMSS-ECM),
#'   or on the top-level slab precision (GMSS-VBECM), or on the regression coefficients' precision (GMN).
#'   (7) \code{a_o} and \code{b_o}: shape parameters of a beta prior on the probabilities of including node-level variables.
#'   (8) \code{s0_v} and \code{s1}: a grid of unscaled spike variances and an unscaled slab variance in the top-level continuous spike-and-slab.
#'   If NULL or any hyperparameter specification is missing,
#'   they will be set to the default: \code{lambda = 2, v0_v = s0_v = seq(1e-4,  1, length.out = 16), v1 = s1 = 100},
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
#' @param full_output Logical: if \code{FALSE}, only model parameter estimates under the selected spike-and-slab configuration will be returned.
#' Otherwise, estimates for all the explored configurations will be returned.
#'
#' @param numCores Integer: number of cores for parallel computation. The default is to use all available cores on the machine.
#'
#' @details The \code{navigm} function implements a Gaussian graphical spike-and-slab model
#' that enables the incorporation and selection of node-level auxiliary variables,
#' thereby enhancing the detection of conditional dependence. Inference is
#' conducted using a scalable (variational) expectation maximisation
#' algorithm, making it suitable for high-dimensional graphs and exploration of multiple spike variances.
#'
#' @return A list containing the following quantities:
#' \describe{
#' \item{estimates}{A list containing model estimates in the model with selected spike variance:
#' \code{Omega} in all the models. In VBECM, \code{a_tau, b_tau}, \code{m_delta}. In ECM, \code{tau_1}, \code{P1}.
#' In GM-VBECM (version 1), \code{a_rho, b_rho}. In GM-ECM (version 1), \code{rho}.
#' In GM-VBECM (version 2), \code{mu_zeta, sig2_inv_zeta} In GM-ECM (version 2), \code{zeta}.
#' In GMN-VBECM, \code{mu_zeta, sig2_inv_zeta}, \code{mu_beta, sig2_inv_beta}, \code{a_sigma, b_sigma}.
#' In GMN-ECM, \code{zeta}, \code{beta}, \code{tau2}.
#' In GMSS-VBECM, \code{mu_zeta, sig2_inv_zeta}, \code{mu_beta, sig2_inv_beta}, \code{a_o, b_o}, \code{a_sigma, b_sigma}.
#' In GMSS-ECM, \code{zeta}, \code{beta}, \code{o}, \code{tau2}.
#' See \code{list_init} for parameter interpretations, except that
#' \code{m_delta, P1} refer to a P x P matrix containing posterior probabilities of including edges, and
#' \code{m_gamma, P2} refer to a vector of length Q containing posterior probabilities of including auxiliary variables.
#'
#' }
#' \item{debugs}{A list containing terms for the purpose of debugging in the model with selected spike variance:
#' In VBECM, \code{n_warning}, \code{vec_n_warning_VB}, \code{list_ELBO}, \code{vec_ELBO_CM}.
#' In ECM, \code{n_warning}, \code{vec_ELBO_CM}.
#' Interpretations:
#' (1) \code{n_warning}: a scalar showing the number of objective function drops in the conditional maximisation step.
#' (2) \code{vec_n_warning_VB}: a vector of length \code{it} showing the number of objective function drops within each variational step.
#' (3) \code{vec_ELBO_CM}: a vector of length \code{it} that tracks objective functions after each conditional maximisation step.
#' (4) \code{list_ELBO}: a list of length \code{it} of which each tracks objective function within each variational step.}
#' \item{it}{Scalar: total number of iterations in the model with selected spike variance. }
#' \item{vec_VB_it}{A vector of length \code{it} containing the number of iterations within each variational step in the model with selected spike variance.}
#' \item{pt}{Scalar: algorithm runtime (seconds) in the model with selected spike variance. }
#' \item{args}{A list containing the inputs to the model with selected spike variance.}
#' \item{index}{Integer(s): an index (or indices in GMSS-ECM) for the selected spike variance.}
#' \item{vec_criterion}{A vector of length ``length(v0_v)`` (or a matrix of dimensions ``length(v0_v)`` X ``length(s0_v)`` in the GMSS-ECM)
#' containing the evaluations of the model selection criterion on a grid (or two grids) of spike variances.}
#' \item{total_pt}{Scalar: total algorithm runtime (seconds). }
#' \item{full_outputs}{A list containing the outputs from \code{navigm_core} for all the explored models. It is non-``NULL`` if ``full_output = T``.}

#' }
#'
#'
#' @examples
#' # simulate data
#'
#' seed <- 123; set.seed(seed)
#' N <- 100; P <- 5; Q <- 3
#' V <- matrix(rnorm(P*Q), nrow = P, ncol = Q)
#' Y <- matrix(rnorm(N*P), nrow = N, ncol = P)
#'
#' # estimate a precision matrix based on Y and effects of node-level variables V at the same time.
#'
#' res_navigm <- navigm(Y, V, numCores = 2)
#' # res_navigm <- navigm(Y, V, method = 'GMN')
#' # res_navigm <- navigm(Y, V, method = 'GM', version = 1)
#' # res_navigm <- navigm(Y, V, method = 'GM', version = 2)
#' # res_navigm <- navigm(Y, V, inference = 'ECM')
#' # res_navigm <- navigm(Y, V, method = 'GMN', inference = 'ECM')
#' # res_navigm <- navigm(Y, V, method = 'GM', inference = 'ECM', version = 1)
#' # res_navigm <- navigm(Y, V, method = 'GM', inference = 'ECM', version = 2)
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @export
#'


navigm <- function(Y, V =NULL,
                   method = 'GMSS',
                   inference = 'VBECM',
                   criterion = 'AIC',
                   list_hyper = NULL,
                   list_init = NULL,
                   ne0 = NULL,
                   tol = 1e-3,
                   maxit = 1e5,
                   # transformY = T,# as mean is 0, always centre Y
                   transformV = F,
                   verbose = T,
                   debug = F,
                   version = NULL,
                   full_output = F,
                   numCores = NULL) {

  # Time
  #
  pt <- Sys.time()

  # Set up cores
  #
  if(is.null(numCores)){
    numCores <- parallel::detectCores()
  }
  doParallel::registerDoParallel(numCores)

  #
  # Set up the grid of spike variances
  #
  list_hyper <- set_default(list_hyper, 'v0_v', seq(1e-4,  1, length.out = 16))
  if(any(list_hyper$vec_v0 <= 0))stop("all the entries of vec_v0 must be positive.")

  if(method =='GMSS' & inference == 'ECM'){
    list_hyper <- set_default(list_hyper, 's0_v', seq(1e-4,  1, length.out = 16))
    if(any(list_hyper$vec_s0 <= 0))stop("all the entries of vec_s0 must be positive.")
  }

  #
  if (verbose) cat("== Parallel exploration of a", ifelse(method =='GMSS' & inference == 'ECM', "double", ""),
                   "grid of spike standard deviations \n v0 = ", list_hyper$v0_v, "\n",
                   ifelse(method =='GMSS' & inference == 'ECM',paste0("and s0 = ", paste(list_hyper$s0_v, collapse = ' '), '\n') , ""),
                   "on ",numCores," cores ... \n\n")

  s0 <- v0 <- NULL
  # keep all the checks & steps inside navigm_core
  # to make it work standalone
  #
  if(!(method =='GMSS' & inference == 'ECM')){
    out <- foreach (v0 = list_hyper$v0_v) %dopar% {
      list_hyper$v0 <- v0
      navigm_core(Y = Y,
                  V = V,
                  method = method,
                  inference = inference,
                  list_hyper = list_hyper,
                  list_init = list_init,
                  ne0 = ne0,
                  tol = tol,
                  maxit = maxit,
                  verbose = verbose,
                  debug = debug,
                  version = version,
                  transformV = transformV)
    }
  }else{
    # double grid search GMSS-ECM
    out <- foreach (v0 = list_hyper$v0_v) %:%
      foreach(s0 = list_hyper$s0_v)%dopar% {
        list_hyper$v0 <- v0
        list_hyper$s0 <- s0
        navigm_core(Y = Y,
                    V = V,
                    method = method,
                    inference = inference,
                    list_hyper = list_hyper,
                    list_init = list_init,
                    ne0 = ne0,
                    tol = tol,
                    maxit = maxit,
                    verbose = verbose,
                    debug = debug,
                    version = version,
                    transformV = transformV)
      }
  }


  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Select from a grid of spike variance using", criterion," ... \n\n")

  if(criterion == 'AIC'){

    if(!(method =='GMSS' & inference == 'ECM')){
      vec_criterion <- sapply(out, function(x){AIC_GSS(x$estimates,N = nrow(x$args$Y))})
    }else{
      vec_criterion <- do.call('rbind',lapply(out, function(x){
        sapply(x, function(y){
          AIC_GSS(y$estimates,N = nrow(y$args$Y))
        })
      }))
      # each row represents v0
    }


  }else if(criterion == 'BIC'){

    if(!(method =='GMSS' & inference == 'ECM')){
      vec_criterion <- sapply(out, function(x){BIC_GSS(x$estimates,N = nrow(x$args$Y))})
    }else{
      vec_criterion <- do.call('rbind',lapply(out, function(x){
        sapply(x, function(y){
          BIC_GSS(y$estimates,N = nrow(y$args$Y))
        })
      }))
    }

  }else if(criterion == 'EBIC'){


    if(!(method =='GMSS' & inference == 'ECM')){
      vec_criterion <- sapply(out, function(x){EBIC_GSS(x$estimates,
                                                        N = nrow(x$args$Y))})
    }else{
      vec_criterion <- do.call('rbind',lapply(out, function(x){
        sapply(x, function(y){
          EBIC_GSS(y$estimates,
                   N = nrow(y$args$Y))
        })
      }))
    }

  }
  if(!(method =='GMSS' & inference == 'ECM')){

    index <- which.min(vec_criterion)

    if(index == 1){
      warning('The selected v0 reaches the lower bound in the grid. Consider extend the grid to lower values. ')
    }else if(index == length(list_hyper$v0_v)){
      warning('The selected v0 reaches the upper bound in the grid. Consider extend the grid to higher values. ')
    }

  }else{

    index <- which(vec_criterion == min(vec_criterion), arr.ind = TRUE)

    if(length(index) != 2){
      warning('More than one best index choices \n')
      if(verbose){
        cat('These indices are all equally the best: \n')
        print(index)
      }
      # avoid selecting on two ends
      tmp <- which(apply(index, 1, function(x)sum(x==1 | x==length(list_hyper$v0_v))) == 0)
      if(length(tmp) == 0){
        # all rows reach boundary
        index <- index[sample(nrow(index),1),]
      }else{
        index <- index[sample(tmp,1),]
      }
    }

    if(index[1] == 1){
      warning('The selected v0 reaches the lower bound in the grid. Consider extend the grid to lower values. ')
    }else if(index[1] == length(list_hyper$v0_v)){
      warning('The selected v0 reaches the upper bound in the grid. Consider extend the grid to higher values. ')
    }
    if(index[2] == 1){
      warning('The selected s0 reaches the lower bound in the grid. Consider extend the grid to lower values. ')
    }else if(index[2] == length(list_hyper$s0_v)){
      warning('The selected s0 reaches the upper bound in the grid. Consider extend the grid to higher values. ')
    }
  }



  if (verbose) cat("... done. == \n\n")

  if(!(method =='GMSS' & inference == 'ECM')){
    if(verbose) cat("Select the index", index, " i.e., v0 = ", list_hyper$v0_v[index], ", the best",criterion, " = ", vec_criterion[index], '.\n\n')
  }else{
    if(verbose) cat("Select the index", index, " i.e., v0 = ", list_hyper$v0_v[index[1]] ," and s0 = ", list_hyper$s0_v[index[2]], ",the best",criterion, " = ", vec_criterion[index[1], index[2]], '.\n\n')
  }


  pt <- Sys.time() - pt
  cat('Total runtime: ',format(pt), '.\n')

  if(!(method =='GMSS' & inference == 'ECM')){

    ans <- out[[index]]

  }else{

    ans <- out[[index[1]]][[index[2]]]

  }


  ans <- c(ans,
           list(
             index = index,
             model_criterion = list(type = criterion,
                                    value = vec_criterion),
             total_pt = pt
           ))


  if(full_output){

    ans$full_output <- out

  }

  return(ans)
}


