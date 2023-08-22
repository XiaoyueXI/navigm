# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm


#' Node-level auxiliary variable for improved Gaussian spike-and-slab with model selection.
#'
#' This function performs simultaneous inference of Gaussian spike-and-slab and node-level auxiliary variable.
#' The graph-level continuous spike-and-slab priors require the exploration of a grid of spike variances and model selections.
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
#' @param criterion Character: the model selection criterion to use, including
#'  'AIC' (default), 'BIC' and 'EBIC'.
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
#' @param full_output Logical; if \code{FALSE}, only return final estimates of model parameters.
#'  Otherwise, return estimates for all the explored models.
#'
#' @param numCores number of cores for the parallel computation. The default is to use all the available cores in the machine.
#'
#' @details \code{navigm} implements a Gaussian graphical model
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
#' \item{args}{A list containing the input arguments in the final selected model.}
#' \item{full_outputs}{A list contains the outputs from \code{navigm_core} for all the explored models.}
#' \item{index}{Index (or indices for GMSS-EM) of the selected spike variance.}
#' \item{vec_criterion}{A list contains the model seletion criterion used and its evaluations on a grid of spike variances. }
#' \item{total_pt}{Total algorithm runtime in seconds. }
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
#' res_navigm <- navigm(Y, V, numCores = 2)
#' # res_navigm <- navigm(Y, V, method = 'GMN')
#' # res_navigm <- navigm(Y, V, method = 'GM', version = 1)
#' # res_navigm <- navigm(Y, V, method = 'GM', version = 2)
#' # res_navigm <- navigm(Y, V, inference = 'EM')
#' # res_navigm <- navigm(Y, V, method = 'GMN', inference = 'EM')
#' # res_navigm <- navigm(Y, V, method = 'GM', inference = 'EM', version = 1)
#' # res_navigm <- navigm(Y, V, method = 'GM', inference = 'EM', version = 2)
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @export
#'


navigm <- function(Y, V =NULL,
                   method = 'GMSS',
                   inference = 'VBEM',
                   criterion = 'AIC',
                   list_hyper = NULL,
                   list_init = NULL,
                   ne0 = NULL,
                   tol = 0.1,
                   maxit = 1000,
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

  if(method =='GMSS' & inference == 'EM'){
    list_hyper <- set_default(list_hyper, 's0_v', seq(1e-4,  1, length.out = 16))
    if(any(list_hyper$vec_s0 <= 0))stop("all the entries of vec_s0 must be positive.")
  }

  #
  if (verbose) cat("== Parallel exploration of a", ifelse(method =='GMSS' & inference == 'EM', "double", ""),
                   "grid of spike standard deviations \n v0 = ", list_hyper$v0_v, "\n",
                   ifelse(method =='GMSS' & inference == 'EM',paste0("and s0 = ", paste(list_hyper$s0_v, collapse = ' '), '\n') , ""),
                   "on ",numCores," cores ... \n\n")

  s0 <- v0 <- NULL
  # keep all the checks & steps inside navigm_core
  # to make it work standalone
  #
  if(!(method =='GMSS' & inference == 'EM')){
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
    # double grid search GMSS-EM
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

    if(!(method =='GMSS' & inference == 'EM')){
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

    if(!(method =='GMSS' & inference == 'EM')){
      vec_criterion <- sapply(out, function(x){BIC_GSS(x$estimates,N = nrow(x$args$Y))})
    }else{
      vec_criterion <- do.call('rbind',lapply(out, function(x){
        sapply(x, function(y){
          BIC_GSS(y$estimates,N = nrow(y$args$Y))
        })
      }))
    }

  }else if(criterion == 'EBIC'){


    if(!(method =='GMSS' & inference == 'EM')){
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
  if(!(method =='GMSS' & inference == 'EM')){

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

  if(!(method =='GMSS' & inference == 'EM')){
    if(verbose) cat("Select the index", index, " i.e., v0 = ", list_hyper$v0_v[index], ", the best",criterion, " = ", vec_criterion[index], '.\n\n')
  }else{
    if(verbose) cat("Select the index", index, " i.e., v0 = ", list_hyper$v0_v[index[1]] ," and s0 = ", list_hyper$s0_v[index[2]], ",the best",criterion, " = ", vec_criterion[index[1], index[2]], '.\n\n')
  }


  pt <- Sys.time() - pt
  cat('Total runtime: ',format(pt), '.\n')

  if(!(method =='GMSS' & inference == 'EM')){

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


