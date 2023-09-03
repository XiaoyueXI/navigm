# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm
#
# Internal functions for plotting and summarising model outputs


#' Compute AIC for a Gaussian graphical spike-and-slab model.
#'
#' This function evaluates one of the model selection criteria, Akaike information criterion (AIC),
#' for a Gaussian graphical spike-and-slab model.
#'
#' @param estimates A list of parameter estimates using the \code{navigm} function,
#' which includes the precision matrix (\code{Omega}), edge posterior inclusion probability (\code{m_delta}),
#' and a matrix \code{S (=Y^T Y)} evaluated based on data Y.
#' @param N number of samples, i.e., number of rows in the input data matrix Y of \code{navigm}.
#'
#' @return Scalar: Akaike information criterion for a Gaussian graphical spike-and-slab model on data with N samples.
#' @examples
#' # simulate data
#'
#' seed <- 123; set.seed(seed)
#' N <- 100; P <- 5; Q <- 3
#' V <- matrix(rnorm(P*Q), nrow = P, ncol = Q)
#' Y <- matrix(rnorm(N*P), nrow = N, ncol = P)
#'
#' # estimate precision matrix based on Y and meanwhile leverage node-level variables V
#' ans <- navigm_core(Y, V)
#' AIC_GSS(ans$estimates,N)
#'
#' @export
AIC_GSS <- function(estimates, N){

  stopifnot('P1'%in% names(estimates) | 'm_delta'%in% names(estimates))
  if('P1'%in% names(estimates)){
    estimates$m_delta <- estimates$P1
  }

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    2 * sum(estimates$m_delta > 0.5)

}

#' Compute BIC for a Gaussian graphical spike-and-slab model.
#'
#' This function evaluates one of the model selection criteria, Bayesian information criterion (BIC),
#' for a Gaussian graphical spike-and-slab model.
#'
#' @param estimates A list of parameter estimates using the \code{navigm} function,
#' which includes the precision matrix (\code{Omega}), edge posterior inclusion probability (\code{m_delta}),
#' and a matrix \code{S (=Y^T Y)} evaluated based on data Y.
#' @param N number of samples, i.e., number of rows in the input data matrix Y of \code{navigm}.
#'
#' @return Scalar: Bayesian information criterion for a Gaussian graphical spike-and-slab model on data with N samples.
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
#' ans <- navigm_core(Y, V)
#' BIC_GSS(ans$estimates,N)
#'
#' @export
BIC_GSS <- function(estimates, N){

  stopifnot('P1'%in% names(estimates) | 'm_delta'%in% names(estimates))
  if('P1'%in% names(estimates)){
    estimates$m_delta <- estimates$P1
  }

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    log(N) * sum(estimates$m_delta > 0.5)

}

#' Compute EBIC for a Gaussian graphical spike-and-slab model.
#'
#' This function evaluates one of the model selection criteria, extended Bayesian information criterion (EBIC),
#' for a Gaussian graphical spike-and-slab model.
#'
#' @param estimates A list of parameter estimates using the \code{navigm} function,
#' which includes the precision matrix (\code{Omega}), edge posterior inclusion probability (\code{m_delta}),
#' and a matrix \code{S (=Y^T Y)} evaluated based on data Y.
#' @param N Number of samples, i.e., number of rows in the input data matrix Y of \code{navigm}.
#' @param gamma A tuning parameter in EBIC valued within [0, 1].
#' \code{gamma = 0} corresponds to recovering BIC.
#' A positive gamma leads to stronger penalization of large graphs.
#' https://arxiv.org/pdf/1011.6640.pdf suggests that
#' \code{gamma = 0.5} (default) achieves a good compromise between positive selection rates and false discovery rates.
#' @return Scalar: extended Bayesian information criterion for a Gaussian graphical spike-and-slab model on data with N samples.
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
#' ans <- navigm_core(Y, V)
#' EBIC_GSS(ans$estimates,N)
#'
#' @export
EBIC_GSS <- function(estimates, N, gamma =0.5){

  stopifnot('P1'%in% names(estimates) | 'm_delta'%in% names(estimates))
  if('P1'%in% names(estimates)){
    estimates$m_delta <- estimates$P1
  }

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  P <- nrow(Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    (log(N) + 4 * gamma * log(P)) * sum(estimates$m_delta > 0.5)

}

#' Plot pROC curve with 95\% confidence interval.
#'
#' This function plots an average receiver operating characteristic (ROC) curve with 95% confidence interval.
#' A partial ROC (pROC) curve can be specified with the argument \code{fpr_stop}.
#'
#' @param ppi A vector of continuous prediction scores, such as posterior inclusion probability, or a list containing such vectors.
#' @param pat A vector of binary true outcomes, or a list containing such vectors.
#' @param fpr_stop Scalar: false positive rate at which the ROC curve is truncated.
#' @param nci Number of confidence bars. The 95\% confidence interval is only displayed if \code{ppi} and \code{pat} are lists.
#' @param ... Other plotting arguments, such as \code{col} and \code{add}. For full details, please refer to \code{ROCR::plot}.
#'
#' @examples
#' set.seed(123); pat <- sample(c(0,1), 10, replace = TRUE)
#' ppi <- sapply(pat,function(x){if(x==0){runif(1,0,0.1)}else{runif(1,0.9,1)}})
#' plot_roc(ppi, pat)
#'
#' @import ROCR
#' @export
plot_roc <- function(ppi, pat, fpr_stop = 1, nci = 11, ...) {

  pred <- ROCR::prediction(ppi,pat)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")

  plot(perf,
       type = "l",
       xlim = c(0,fpr_stop),
       ylim = c(0, 1),
       avg= "vertical",
       spread.estimate="stderror",
       spread.scale = 2,
       show.spread.at = seq(0,fpr_stop,length.out= nci),
       ...)

}

#' Compute AUC.
#'
#' This function computes the area under the ROC curve (AUC).
#' A partial AUC can be specified with the argument \code{fpr_stop}.
#'
#' @param ppi A vector of continuous prediction scores, such as posterior inclusion probability, or a list containing such vectors.
#' @param pat A vector of binary true outcomes, or a list containing such vectors.
#' @param fpr_stop Scalar: false positive rate at which the ROC curve is truncated.
#' @param standardise Logical: if set to FALSE (default), standardisation will not be applied;
#' if set to TRUE, standardisation will be applied following https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3744586/.
#'
#' @return Scalar: (standardised partial) area under the ROC curve.
#'
#' @examples
#' set.seed(123); pat <- sample(c(0,1), 10, replace = TRUE)
#' ppi <- sapply(pat,function(x){if(x==0){runif(1,0,0.1)}else{runif(1,0.9,1)}})
#' compute_pauc(ppi, pat)
#'
#' @import ROCR
#' @export
compute_pauc <- function(ppi, pat,  fpr_stop = 1, standardise = F) {

  pred <- ROCR::prediction(ppi,pat)
  perf <- performance(pred, measure = 'auc',fpr.stop = fpr_stop)
  pauc <- do.call('c',perf@y.values)

  if(standardise){
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3744586/
    pauc <- (1 + (pauc - fpr_stop^2/2)/(fpr_stop - fpr_stop^2/2))/2
  }

  return(pauc)

}

#' Plot auxiliary-variable-associated quantities.
#'
#' This function generates a scatter plot illustrating the values associated with auxiliary variables \code{ppi} against their corresponding indices \code{ppi_names}.
#'
#' @param ppi A vector of the values associated with auxiliary variables, such as effect sizes and posterior probability of inclusion.
#' @param ppi_names A vector of variable names. If set to \code{NULL} (default), the indices \code{1:length(ppi)} will be used.
#' @param col A character indicating the color of all the points, or a vector containing such characters with the same length as the variables.
#' @param condition A vector of logical values. If \code{condition == TRUE}, vertical lines will be drawn from points to the x-axis.
#' @param xlab Character: a title for the x-axis (the default is "auxiliary variables").
#' @param ylab Character: a title for the y axis (the default is "effect sizes").
#' @param ... Other plotting arguments. Refer to \code{base::plot} for more details.
#'
#' @examples
#' set.seed(123); pat <- sample(c(0,1), 10, replace = TRUE)
#' ppi <- sapply(pat,function(x){if(x==0){runif(1,0,0.1)}else{runif(1,0.9,1)}})
#' plot_ppi(ppi, ylab = 'inclusion probability')
#'
#' @importFrom graphics segments
#' @export
plot_ppi <- function(ppi, ppi_names = NULL,
                     col ='black', condition = (ppi >=0.5),
                     xlab = 'auxiliary variables', ylab = 'effect sizes',...){

  if(is.null(ppi_names)){
    ppi_names <- 1:length(ppi)
  }

  plot(ppi_names, ppi, col = col,
       pch = 19, xlab = xlab, ylab = ylab, ...)

  graphics::segments(which(condition), rep(0,sum(condition)), x1 = which(condition), y1 = ppi[condition])

}


#' Assess performance using a threshold.
#'
#' This function assesses estimation performance using a predefined threshold,
#' including metrics such as true positive rate (TPR), false positive rate (FPR), true negative rate (TNR),
#' false negative rate (FNR), recall, precision, and F1-scores (F1).
#'
#' @param ppi A vector of continuous prediction scores, such as posterior inclusion probability, or a list containing such vectors.
#' @param pat A vector of binary true outcomes, or a list containing such vectors.
#' @param threshold Scalar: a threshold for classifying \code{ppi}.
#'
#' @return A list of performance metrics calculated using the specified threshold.
#'
#' @examples
#' set.seed(123); pat <- sample(c(0,1), 10, replace = TRUE)
#' ppi <- sapply(pat,function(x){if(x==0){runif(1,0,0.1)}else{runif(1,0.9,1)}})
#' compute_perf(ppi, pat)
#'
#' @export
compute_perf <- function(ppi, pat, threshold = 0.5){

  pat <- as.logical(pat)
  tpr <- sum(ppi >= threshold & pat)/sum(pat)
  fpr <- sum(ppi >= threshold & !pat)/sum(!pat)
  fnr <- sum(ppi < threshold & pat)/sum(pat)
  tnr <- sum(ppi < threshold & !pat)/sum(!pat)

  recall <- tpr
  precision <-  sum(ppi >= threshold & pat)/sum(ppi >= threshold)
  f1 <- 2 * (precision * recall) / (precision + recall)

  create_named_list_(tpr, fpr, tnr, fnr, recall, precision, f1)

}

#' Plot a graph.
#'
#' This function generates plots of undirected and unweighted networks based on an adjacency matrix.
#' This visualisation is suitable for small networks.
#'
#' @param x  An adjacency matrix.
#' @param node_names A vector of the same length as the number of rows and columns of \code{x}, containing vertex names.
#' @param cex Scalar: the factor by which vertex label sizes should be enlarged relative to the default size.
#' @param vertex.color Character: vertex color, for example, \code{vertex.color = "yellow"} (default).
#' @param ... Other plotting arguments. For complete details, refer to \code{igraph::plot}.
#'
#' @import igraph
#'
#' @examples
#' A <- matrix(0, 5, 5); A[1,2] <- A[2,1] <- A[2,5] <- A[5,2] <- 1;
#' plot_network(A)
#'
#' @export
plot_network <- function(x, cex = 0.5,
                         node_names = NULL,vertex.color = "yellow" , ...){

  x <- x == 1
  diag(x) <- F

  if(is.null(node_names)){
    if(!is.null(rownames(x))){
      node_names <- rownames(x)
    }else{
      node_names <- seq_len(ncol(x))
    }

  }

  g <- graph_from_adjacency_matrix(x, mode="undirected")
  lay <- layout_in_circle(g)
  V(g)$label.cex <- cex
  pg <- plot(g, layout=lay,
             vertex.color=vertex.color,
             vertex.label = node_names, ...)

}
