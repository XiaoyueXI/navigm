# This file is part of the `navigss` R package:
#     https://github.com/XiaoyueXI/navigss


#' Node-level auxiliary variable informed Gaussian spike-and-slab using EM.
#'
#' This function performs expectation maximization for
#' simultaneous Gaussian spike-and-slab inference and node-level auxiliary variable selections.
#'
#' @param Y Data matrix of dimension N x P, where N is the number of samples and P is the number of nodes in the graph.
#'   \code{Y} will be centred within \code{navigss_em_core} call.
#' @param V Input matrix of dimension P x Q, where Q is the number of candidate auxiliary variables. No need to supply intercept.
#'   \code{V} will be centred within \code{navigss_em_core} call if not within [0,1] by default, unless disabled by users via \code{transformV = F}.
#'   If \code{V} is not specified, will use \code{method = 'GM'}.
#' @param list_hyper A list containing the model hyperparameters:
#'   the gamma prior on the continuous spike-and-slab scale requires shape parameter \code{a_tau} and rate parameter \code{b_tau};
#'   the unscaled spike and slab variances in the continuous spike-and-slab are v0 and v1;
#'   the exponential prior on diagonal entries requires rate parameter lambda;
#'   the normal prior on the overall sparsity requires mean n0 and variance t02;
#'   the gamma prior on the slab precision in the discrete spike-and-slab requires shape parameter \code{a_sigma} and rate parameter \code{b_sigma} (needed in GMN & GMSS);
#'   the beta prior on auxiliary variable inclusion probability requires shape parameters \code{a_o} and \code{b_o} (only needed in GMSS);
#'   the beta prior on edge inclusion probability requires shape parameters \code{a_rho} and \code{b_rho} (only needed in \code{method = "GM"} and \code{version = 1});
#'
#'   if NULL or any hyperparameter specification is missing,
#'   set to default, \code{lambda = 2, v0 = 0.1, v1 = 100, a_tau = b_tau = a_sigma = b_sigma = 2, n0 = -2, t02 = 0.5, a_o = 1, b_o = Q, a_rho = 1, b_rho = P.}.
#' @param list_init A list containing the initial parameters:
#'  the continuous spike-and-slab scale \code{tau_1};
#'  the overall sparsity \code{zeta};
#'  the regression coefficients \code{beta} (needed in GMN & GMSS);
#'  the slab precision in the discrete spike-and-slab \code{tau_2} (needed in GMN & GMSS);
#'  the auxiliary variable inclusion probability \code{o} (only needed in GMSS);
#'  if \code{NULL} or any initialisation is missing,
#'  set to default, \code{beta}: vector of zeros, \code{tau_1 = tau_2 = 1}, \code{zeta = list_hyper$n0}, \code{o = 1 / Q}, \code{rho = 1 / P}.
#' @param tol Scalar: tolerance for the stopping criterion (default is 0.1).
#' @param maxit Scalar: maximum number of iterations allowed (default is 1000).
#' @param transformV Logical; if \code{TRUE} (default) and not ranged within [0,1],
#' \code{V} will be standardised within \code{navigss_vbem_core} call.
#' Otherwise, \code{V} will not be transformed, and can be processed by users.
#' @param seed Scalar: seed for reproducible default choices of initialisation (default is 123).
#' @param verbose Logical; if \code{TRUE} (default), standard verbosity; otherwise, no messages.
#' @param track_ELBO Logical; if \code{TRUE} (default), track objective function in the EM; otherwise, not track.
#' @param debug Logical; if \code{FALSE} (default), not record additional terms for the purpose of debug;
#'  otherwise, track number of warnings of decreasing objective functions, possibly due to numerical round-off.
#' @param version Integer; take values of 1 (a beta prior on edge inclusion) or 2 (a normal prior on probit edge inclusion) only. Only valid when \code{method = 'GM'}.
#' @details \code{navigss_em_core} implements a Gaussian graphical model
#'   that allows incorporating and selecting node-level auxiliary variables,
#'   thereby enhancing the detection of conditional dependence. Inference is
#'   carried out using a scalable expectation maximisation
#'   algorithm, which is applicable to high-dimension graphs.
#'
#'
#' @return A list containing the following estimates and settings:
#'  \describe{
#'  \item{estimates}{A list containing parameter estimation, which includes a selection of following terms depending on method:
#'  \code{Omega} (estimated precision matrix of dimension P x P),
#'  \code{P1} (posterior inclusion probability matrix of dimension P x P. off-diagonal entry (i, j) corresponds to
#'  the variational posterior probability of edge (i,j) being included in the graph),
#'  \code{P2} (posterior inclusion probability vector of dimension Q. entry Q corresponds to the variational posterior probability of
#'           qth auxiliary variable explaining the node centrality),
#'  \code{beta} (vector of length Q containing the estimated regression coefficients of
#'           node-level auxiliary variables in columns of matrix \code{V}),
#' \code{tau_1} (estimated continuous spike-and-slab scale),
#' \code{tau_2} (estimated slab precision in discrete spike-and-slab,
#' \code{o} (estimated auxiliary variable inclusion probability),
#' \code{zeta} (estimated overall sparsity parameter).}
#'  \item{debugs}{A list containing terms for debugging:
#'  \code{n_warning} (scalar. number of ELBO drops in the maximisation step).}
#'  \item{args}{A list containing the algorithm inputs including \code{list_hyper} and \code{list_init}.}
#' \item{vec_ELBO}{Vector of length \code{it}. ELBOs after each EM update. }
#' \item{it}{Scalar. Final number of iterations. \code{it == maxit} indicates the algorithm dose not converge. }
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
#' res_navigss_em_core <-navigss_em_core(Y, V, transformV = F)
#' # res_navigss_em_core <-navigss_em_core(Y, V, method = 'GMN',transformV = F)
#' # res_navigss_em_core <-navigss_em_core(Y, V, method = 'GM', version = 1,transformV = F)
#' # res_navigss_em_core <-navigss_em_core(Y, V, method = 'GM', version = 2,transformV = F)
#'
#' @import matrixcalc
#' @export
#'

navigss_em_core <- function(Y, V =NULL,
                              ne0 = NULL,
                              method = 'GMSS',
                              list_hyper = NULL, list_init = NULL,
                              tol = 0.1, maxit = 1000,
                              seed = NULL, verbose = T,
                              track_ELBO = F, debug = F,
                              version = NULL,
                              # transformY = T,# as mean is 0, always centre Y
                              transformV = T) {


  if (verbose != 0){
    cat(paste0("\n======================= \n",
               "== Preprocess data ... == \n",
               "======================= \n\n"))
  }
  if(ncol(Y)!=nrow(V)){
    stop('Columns of Y and rows of V must match.')
  }

  # if(transformY){
  Y <- scale(Y, center = TRUE, scale = FALSE) # center the data
  # }

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
      tmp <- get_n0_t02(ne0[1], ne0[2])
      list_hyper$n0 <-  tmp$n0
      list_hyper$t02 <-  tmp$t02
    }
  }else{
    warning('As no valid list_hyper and ne0 are provided, use default hyperparameters n0 = -2 and t02 = 0.5.')
  }


  if (verbose != 0) cat("... done. == \n\n")


  cat("**************************************************** \n",
      "Number of samples: ", nrow(Y), "\n",
      "Number of nodes: ", ncol(Y), "\n")

  if(method == 'GM'){

    cat("**************************************************** \n\n")

    cat(paste0("======================================== \n",
               "== GM: spike-and-slab graphical model == \n",
               "======================================== \n\n"))

    res_navigss <- gm_em_core(Y,
                                list_hyper = list_hyper,
                                list_init =  list_init,
                                tol = tol,
                                maxit = maxit,
                                verbose = verbose,
                                track_ELBO = track_ELBO,
                                debug = debug,
                                version = version)

  }else if(method == 'GMN'){

    cat("Number of candidate node-level auxiliary variables: ", ncol(V), "\n",
        "**************************************************** \n\n")

    cat(paste0("============================================================================================================== \n",
               "== GMN: spike-and-slab graphical model with normal prior for the node-level auxiliary variable coefficients == \n",
               "============================================================================================================== \n\n"))

    res_navigss <- gmn_em_core(Y,
                                 V,
                                 list_hyper = list_hyper,
                                 list_init =  list_init,
                                 tol = tol,
                                 maxit = maxit,
                                 verbose = verbose,
                                 track_ELBO = track_ELBO,
                                 debug = debug)


  }else if(method == 'GMSS'){

    cat("Number of candidate node-level auxiliary variables: ", ncol(V), "\n",
        "**************************************************** \n\n")

    cat(paste0("====================================================================================================================== \n",
               "== GMSS: spike-and-slab graphical model with spike-and-slab prior for the node-level auxiliary variable coefficients == \n",
               "====================================================================================================================== \n\n"))

    res_navigss <- gmss_em_core(Y,
                                  V,
                                  list_hyper = list_hyper,
                                  list_init =  list_init,
                                  tol = tol,
                                  maxit = maxit,
                                  verbose = verbose,
                                  track_ELBO = track_ELBO,
                                  debug = debug)

  }

  if (verbose) cat("... done. == \n\n")

  return(res_navigss)
}

