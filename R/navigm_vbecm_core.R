# This file is part of the `navigm` R package:
#     https://github.com/XiaoyueXI/navigm


navigm_vbecm_core <- function(Y, V =NULL,
                              method = 'GMSS',
                              list_hyper = NULL,
                              list_init = NULL,
                              tol = 1e-3,
                              maxit = 1e5,
                              verbose = T,
                              debug = F,
                              version = NULL) {


  if(method == 'GM'){

    cat("**************************************************** \n\n")

    cat(paste0("======================================== \n",
               "== GM: spike-and-slab graphical model == \n",
               "======================================== \n\n"))

    ans <- gm_vbecm_core(Y,
                         list_hyper = list_hyper,
                         list_init =  list_init,
                         tol = tol,
                         maxit = maxit,
                         verbose = verbose,
                         debug = debug,
                         version = version)

  }else if(method == 'GMN'){

    cat(paste0("============================================================================================================== \n",
               "== GMN: spike-and-slab graphical model with normal priors for the node-level auxiliary variable coefficients == \n",
               "============================================================================================================== \n\n"))

    ans <- gmn_vbecm_core(Y,
                          V,
                          list_hyper = list_hyper,
                          list_init =  list_init,
                          tol = tol,
                          maxit = maxit,
                          verbose = verbose,
                          debug = debug)


  }else if(method == 'GMSS'){

    cat(paste0("====================================================================================================================== \n",
               "== GMSS: spike-and-slab graphical model with spike-and-slab priors for the node-level auxiliary variable coefficients == \n",
               "====================================================================================================================== \n\n"))

    ans <- gmss_vbecm_core(Y,
                           V,
                           list_hyper = list_hyper,
                           list_init =  list_init,
                           tol = tol,
                           maxit = maxit,
                           verbose = verbose,
                           debug = debug)

  }

  if (verbose) cat("... done. == \n\n")

  return(ans)
}

