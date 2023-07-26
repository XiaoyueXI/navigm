trunc_perf <-  function(perf, fpr.stop) { # perf = list with multiple entry per replication

  for (iperf in seq_along(perf@x.values)){
    ind = which(perf@x.values[[iperf]] <= fpr.stop)
    perf@y.values[[iperf]] = perf@y.values[[iperf]][ind]
    perf@x.values[[iperf]] = perf@x.values[[iperf]][ind]
  }

  perf

}


plot_roc <- function(ppi, pat, col, lty = 1,
                     add = FALSE, fpr_stop = 1,
                     main = "ROC curves") {

  pred <- ROCR::prediction(ppi,pat)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")

  plot(perf, col = col, type = "l", lwd = 2, lty = lty,
       main = main, add = add, xlim = c(0,fpr_stop), ylim = c(0, 1),
       avg= "vertical", spread.estimate="stderror", spread.scale = 2,
       show.spread.at = seq(0,fpr_stop,length.out= 11))

}


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

# plot_ppi(gmss_vbem$estimates$m_gamma, ylab = "PIP")
# plot_ppi(beta_true_gmss, condition = (beta_true_gmss!=0))
plot_ppi <- function(ppi, ppi_names = NULL, col ='black', condition = (ppi >=0.5),
                     xlab = 'auxiliary variables', ylab = 'effect sizes',...){

  if(is.null(ppi_names)){
    ppi_names <- 1:length(ppi)
  }
  plot(ppi_names, ppi, col = col,
       pch = 19, xlab = xlab, ylab = ylab, ...)

  segments(which(condition), rep(0,sum(condition)), x1 = which(condition), y1 = ppi[condition])

}

