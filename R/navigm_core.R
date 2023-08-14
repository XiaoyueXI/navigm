# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm


#' Node-level auxiliary variable informed Gaussian spike-and-slab.
#'
#' This function performs simultaneous Gaussian spike-and-slab inference and node-level auxiliary variable selections.
#'
#' @param Y Data matrix of dimension N x P, where N is the number of samples and P is the number of nodes in the graph.
#'   \code{Y} will be centred within \code{navigm_core} call.
#'
#' @param V Input matrix of dimension P x Q, where Q is the number of candidate auxiliary variables. No need to supply intercept.
#'   If \code{V} is not specified, set \code{method = 'GM'}.
#'
#' @param list_hyper A list containing the model hyperparameters:
#'   the gamma prior on the continuous spike-and-slab scale requires shape parameter \code{a_tau} and rate parameter \code{b_tau};
#'   the unscaled spike and slab variances in the continuous spike-and-slab are \code{v0} and \code{v1};
#'   the exponential prior on diagonal entries requires rate parameter \code{lambda};
#'   the normal prior on the overall sparsity requires mean \code{n0} and variance \code{t02};
#'   the gamma prior on the slab precision in the discrete spike-and-slab requires shape parameter \code{a_sigma} and rate parameter \code{b_sigma} (needed in GMN & GMSS);
#'   the beta prior on auxiliary variable inclusion probability requires shape parameters \code{a_o} and \code{b_o} (only needed in GMSS);
#'   the beta prior on edge inclusion probability requires shape parameters \code{a_rho} and \code{b_rho} (only needed in \code{method = "GM"} and \code{version = 1});
#'
#'   if NULL or any hyperparameter specification is missing,
#'   set to default, \code{lambda = 2, v0 = 0.1, v1 = 100, a_tau = b_tau = a_sigma = b_sigma = 2},
#'   \code{n0 = -2, t02 = 0.5, a_o = 1, b_o = Q, a_rho = 1, b_rho = P}.
#'
#' @param ne0 Vector of length 2 whose entries are the prior expectation and variance of
#'  the number of edges in absence of hubs. Will not be used if \code{n0} in \code{list_hyper} is non-\code{NULL}.
#'
#' @param list_init A list containing the initial variational parameters if \code{inference = 'VBEM'}:
#'  the gamma variational distribution of the continuous spike-and-slab scale includes shape parameter \code{alpha_tau} and rate parameter \code{beta_tau};
#'  the normal variational distribution of the overall sparsity includes mean \code{mu_zeta} and variance \code{sig2_inv_zeta};
#'  the normal variational distribution of regression coefficients includes mean \code{mu_beta} and variance \code{sig2_inv_beta} (needed in GMN & GMSS);
#'  the gamma variational distribution of slab precision in the discrete spike-and-slab includes \code{alpha_sigma, beta_sigma} (needed in GMN & GMSS);
#'  the beta variational distribution of auxiliary variable inclusion probability includes \code{alpha_o} and \code{beta_o} (only needed in GMSS);
#'  the beta variational distribution of edge inclusion probability includes \code{alpha_rho} and \code{beta_rho} (only needed in \code{method = "GM"} and \code{version = 1});
#'  if \code{NULL} or any initialisations are missing,
#'  set to default, \code{mu_beta}: vector of zeros, \code{sig2_inv_beta}: vector of ones,
#'  \code{alpha_tau = beta_tau = alpha_sigma = beta_sigma = 1}, \code{mu_zeta = list_hyper$n0},
#'  \code{sig2_inv_zeta = 1/list_hyper$t02}, \code{a_o = 1, b_o = Q},  \code{a_rho = 1, b_rho = P}.
#'  A list containing the initial parameters if \code{inference = 'EM'}:
#'  the continuous spike-and-slab scale in the level of edge selection \code{tau_1};
#'  the overall sparsity \code{zeta};
#'  the regression coefficients \code{beta} (needed in GMN & GMSS);
#'  the continuous spike-and-slab scale in the level of auxiliary variable selection \code{tau_2} (needed in GMN & GMSS);
#'  the auxiliary variable inclusion probability \code{o} (only needed in GMSS);
#'  if \code{NULL} or any initialisation is missing,
#'  set to default, \code{beta}: vector of zeros, \code{tau_1 = tau_2 = 1}, \code{zeta = list_hyper$n0}, \code{o = 1 / Q}, \code{rho = 1 / P}.
#'
#' @param method Character: the method to be used, including
#'  'GM' (vanilla spike-and-slab graphical model),
#'  'GMN' (spike-and-slab graphical model with normal prior for the node-level auxiliary variable coefficients), and
#'  'GMSS' (default; spike-and-slab graphical model with spike-and-slab prior for the node-level auxiliary variable coefficients).
#'
#' @param inference Character: the inference algorithm to be used, including
#'  'EM' (expectation maximisation algorithm),
#'  'VBEM' (default; variational Bayes expectation maximisation algorithm),
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
#' @param track_ELBO Logical; if \code{TRUE} (default), track ELBO after each maximisation step; otherwise, not track.
#'
#' @param debug Logical; if \code{FALSE} (default), not record additional terms for the purpose of debug;
#'  otherwise, record the number of decreasing ELBOs in the maximisation step, possibly due to numerical round-off;
#'  and ELBOs and number of decreasing ELBOs within each variational update (only in the VBEM).
#'
#' @param version Integer; take values of 1 (a beta prior on edge inclusion) or 2 (a normal prior on probit edge inclusion) only.
#'  Only valid when \code{method = 'GM'}.
#'
#' @details \code{navigm_core} implements a Gaussian graphical model
#'   that allows incorporating and selecting node-level auxiliary variables,
#'   thereby enhancing the detection of conditional dependence. Inference is
#'   carried out using a scalable (variational) expectation maximisation
#'   algorithm, which is applicable to high-dimension graphs.
#'
#'
#' @return A list containing the following estimates and settings:
#'  \describe{
#'  \item{estimates}{A list containing estimates of variational distribution parameters, which includes a selection of following terms depending on method:
#'  \code{Omega} (estimated precision matrix of dimension P x P),
#'  \code{m_delta} (VBEM), \code{P1} (EM) (posterior inclusion probability matrix of dimension P x P. off-diagonal entry (i, j) corresponds to
#'  the variational posterior probability of edge (i,j) being included in the graph),
#'  \code{m_gamma} (VBEM), \code{P2} (EM) (posterior inclusion probability vector of dimension Q. entry Q corresponds to the variational posterior probability of
#'           qth auxiliary variable explaining the node centrality in GMSS),
#'  \code{mu_beta, sig2_inv_beta} (VBEM) (vector of length Q containing the variational posterior mean and inverse variance of regression coefficients of
#'           node-level auxiliary variables in columns of matrix \code{V} in GMN and GMSS),
#' \code{a_tau, b_tau} (VBEM) (estimated parameters of Gamma distributed scale in continuous spike-and-slab),
#' \code{a_sigma, b_sigma} (VBEM) (estimated parameters of Gamma distributed slab variance in discrete spike-and-slab in GMN and GMSS,
#' \code{a_o, b_o} (VBEM) (estimated parameters of Beta distributed auxiliary variable inclusion probability in GMSS),
#' \code{mu_zeta, sig2_inv_zeta} (VBEM) (variational posterior mean and inverse variance of overall sparsity parameter except for \code{GM version = 1} ).
#' \code{beta} (EM) (vector of length Q containing the estimated regression coefficients of
#'           node-level auxiliary variables in columns of matrix \code{V} in GMN and GMSS),
#' \code{tau_1} (EM) (estimated edge-level continuous spike-and-slab scale),
#' \code{tau_2} (EM) (estimated variable-level continuous spike-and-slab scale in GMN and GMSS),
#' \code{o} (EM) (estimated auxiliary variable inclusion probability in GMSS),
#' \code{zeta} (EM) (estimated overall sparsity parameter except for \code{GM version = 1}).}
#}
#'  \item{debugs}{A list containing terms for debugging:
#'  \code{n_warning} (scalar. number of ELBO drops in the maximisation step),
#'  \code{vec_n_warning_VB} (vector of length \code{it}. number of ELBO drops in each variational step. Only in VBEM),
#'  \code{list_elbo} (a list of length \code{it} containing the ELBO in each variational step. Only in VBEM).  }
#'  \item{args}{A list containing the algorithm inputs including \code{list_hyper} and \code{list_init}.}
#' \item{it}{Scalar. Final number of iterations. \code{it == maxit} indicates the algorithm dose not converge. }
#' \item{vec_ELBO_M}{Vector of length \code{it}. ELBOs after the maximisation step. }
#' \item{vec_VB_it}{Vector of length \code{it}. Number of iterations in each variational step.
#' \code{any(vec_VB_it == maxit)} indicates the variational step dose not converge. Only in VBEM.}
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
                         track_ELBO = F, debug = F,
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


