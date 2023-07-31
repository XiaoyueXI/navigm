# trunc_perf <-  function(perf, fpr.stop) { # perf = list with multiple entry per replication
#
#   for (iperf in seq_along(perf@x.values)){
#     ind = which(perf@x.values[[iperf]] <= fpr.stop)
#     perf@y.values[[iperf]] = perf@y.values[[iperf]][ind]
#     perf@x.values[[iperf]] = perf@x.values[[iperf]][ind]
#   }
#
#   perf
#
# }

#' Evaluate the Akaike information criterion (AIC) for spike-and-slab Gaussian graphical models.
#'
#' This function evaluates one of the model selection criteria, Akaike information criterion,
#' for spike-and-slab Gaussian graphical models.
#'
#' @param estimates A list of parameter estimates using the \code{navigm} function,
#' which includes precision matrix (Omega), edge posterior inclusion probability (m_delta)
#' and data induced matrix S (=Y^T Y).
#' @param N number of samples (number of rows in Y).
#'
#' @return Scalar. Akaike information criterion of a spike-and-slab Gaussian graphical model on N-sample data Y.
#'
#' @export
AIC_GSS <- function(estimates, N){

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    2 * sum(estimates$m_delta > 0.5)

}

#' Evaluate the Bayesian information criterion (BIC) for spike-and-slab Gaussian graphical models.
#'
#' This function evaluates one of the model selection criteria, Bayesian information criterion,
#' for spike-and-slab Gaussian graphical models.
#'
#' @param estimates A list of parameter estimates using the \code{navigm} function,
#' which includes precision matrix (Omega), edge posterior inclusion probability (m_delta)
#' and data induced matrix S (=Y^T Y).
#' @param N number of samples (number of rows in Y).
#'
#' @return Scalar. Bayesian information criterion of a spike-and-slab Gaussian graphical model on N-sample data Y.
#'
#' @export
BIC_GSS <- function(estimates, N){

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    log(N) * sum(estimates$m_delta > 0.5)

}

#' Evaluate the extended Bayesian information criterion (BIC) for spike-and-slab Gaussian graphical models.
#'
#' This function evaluates one of the model selection criteria, extended Bayesian information criterion,
#' for spike-and-slab Gaussian graphical models.
#'
#' @param estimates A list of parameter estimates using the \code{navigm} function,
#' which includes precision matrix (Omega), edge posterior inclusion probability (m_delta)
#' and data induced matrix S (=Y^T Y).
#' @param gamma EBIC parameter in [0,1] (default = 0.5). \code{gamma = 0} recovers BIC.
#' Positive gamma leads to stronger penalization of large graphs.
#' https://arxiv.org/pdf/1011.6640.pdf shows \code{gamma=0.5} achieves a good compromise between positive selection rates and false discovery rates.
#' @param N number of samples (number of rows in Y).
#' @param P number of nodes (number of columns in Y).
#'
#' @return Scalar. Extended Bayesian information criterion of a spike-and-slab Gaussian graphical model on N-sample data Y.
#'
#'
#' @export
EBIC_GSS <- function(estimates, gamma =0.5, N, P){

  Omega <- estimates$Omega
  Omega[estimates$m_delta <= 0.5] <- 0
  diag(Omega) <- diag(estimates$Omega)

  sum(diag(estimates$S %*% Omega)) -
    N * determinant(Omega, logarithm = TRUE)$modulus[1] +
    (log(N) + 4 * gamma * log(P)) * sum(estimates$m_delta > 0.5)

}

#' Plot truncated receiver operating characteristic curves.
#'
#' This function plots the truncated receiver operating characteristic (ROC) curves at
#' the false positive rate \code{fpr_stop}. ROC curves describe the
#' performance of a classification model at all classification thresholds.
#'
#' @param ppi A vector of continuous prediction scores such as posterior inclusion probability or a list of such vectors.
#' @param pat A vector of binary true outcomes or a list of such vectors.
#' @param fpr_stop Scalar. False positive rate at which the ROC curve is truncated.
#' @param nci number of confidence bars if \code{ppi} and \code{pat} are lists.
#' @param ... other plotting arguments, such as \code{col}, \code{add}. Full details see \code{ROCR}.
#'
#' @import ROCR
#' @export
plot_roc <- function(ppi, pat, fpr_stop = 1, nci = 11, ...) {

  pred <- ROCR::prediction(ppi,pat)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")

  plot(perf,
       type = "l", xlim = c(0,fpr_stop), ylim = c(0, 1),
       avg= "vertical", spread.estimate="stderror", spread.scale = 2,
       show.spread.at = seq(0,fpr_stop,length.out= nci),
       ...)

}

#' Compute standardised partial area under the ROC curve.
#'
#' This function compute the area under the ROC curve (AUC) truncated at
#' the false positive rate \code{fpr_stop}. AUC is a threshold-free performance measure.
#'
#' @param ppi A vector of continuous prediction scores such as posterior inclusion probability or a list of such vectors.
#' @param pat A vector of binary true outcomes or a list of such vectors.
#' @param fpr_stop Scalar. False positive rate at which the ROC curve is truncated.
#' @param standardise Logical. If FALSE (default), not standardise the partial AUC; otherwise, standardise.
#'
#' @return A scalar measuring the classification performance.

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

#' Plot variable-specific quantities.
#'
#' This function plots variable-specific quantity \code{ppi} against their indices \code{ppi_names},
#' for instance posterior inclusion probabilities or effect sizes of auxiliary variables
#'
#' @param ppi A vector of the variable-specific quantity.
#' @param ppi_names A vector of variable indices. If \code{NULL} (default), use the numbered indices.
#' @param col A character indicating the point color or a vector of the same length of variables.
#' @param condition A vector of logical. If satisfying \code{condition == T}, draw vertical lines from points to the x-axis.
#' @param xlab Character. A title for the x axis.
#' @param ylab Character. A title for the y axis.
#' @param ... other plotting arguments. See \code{base::plot}.
#'
#' @export
plot_ppi <- function(ppi, ppi_names = NULL, col ='black', condition = (ppi >=0.5),
                     xlab = 'auxiliary variables', ylab = 'effect sizes',...){

  if(is.null(ppi_names)){
    ppi_names <- 1:length(ppi)
  }
  plot(ppi_names, ppi, col = col,
       pch = 19, xlab = xlab, ylab = ylab, ...)

  segments(which(condition), rep(0,sum(condition)), x1 = which(condition), y1 = ppi[condition])

}


#' Compute thresholds-based performance measures.
#'
#' This function compute the thresholds-based performance measures
#' including true positive rate (tpr), false positive rate (fpr), true negative rate (tnr),
#' false negative rate (fnr), recall, precision, F1-scores (f1).
#'
#' @param ppi A vector of continuous prediction scores such as posterior inclusion probability or a list of such vectors.
#' @param pat A vector of binary true outcomes or a list of such vectors.
#' @param threshold Scalar. A threshold for classification.
#'
#' @return A list of thresholds-based performance measures.

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


#' This function plots the undirected and unweighted networks based on an adjacency matrix (typically useful when the network is small).
#'
#' @param x A adjacency matrix.
#' @param cex Scalar. The amount by which vertex label sizes should be magnified relative to the default.
#' @param node_names A vector of same length as the number of rows & columns of \code{x} and containing vertex names.
#' @import igraph
#' @export
plot_network <- function(x, cex = 0.5, node_names = NULL, ...){

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
  pg <- plot(g, layout=lay, vertex.color="yellow", vertex.label = node_names)

  return(pg)

}
