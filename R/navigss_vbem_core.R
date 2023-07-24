# This file is part of the `navigss` R package:
#     https://github.com/XiaoyueXI/navigss


navigss_vbem_core <- function(Y, V =NULL,
                              method = 'GMSS',
                              list_hyper = NULL, list_init = NULL,
                              tol = 0.1, maxit = 1000,
                              verbose = T,
                              track_ELBO = F, debug = F,
                              version = NULL) {


  if(method == 'GM'){

    cat("**************************************************** \n\n")

    cat(paste0("======================================== \n",
               "== GM: spike-and-slab graphical model == \n",
               "======================================== \n\n"))

    ans <- gm_vbem_core(Y,
                                list_hyper = list_hyper,
                                list_init =  list_init,
                                tol = tol,
                                maxit = maxit,
                                verbose = verbose,
                                track_ELBO = track_ELBO,
                                debug = debug,
                                version = version)

  }else if(method == 'GMN'){

    cat(paste0("============================================================================================================== \n",
               "== GMN: spike-and-slab graphical model with normal prior for the node-level auxiliary variable coefficients == \n",
               "============================================================================================================== \n\n"))

    ans <- gmn_vbem_core(Y,
                                 V,
                                 list_hyper = list_hyper,
                                 list_init =  list_init,
                                 tol = tol,
                                 maxit = maxit,
                                 verbose = verbose,
                                 track_ELBO = track_ELBO,
                                 debug = debug)


  }else if(method == 'GMSS'){

    cat(paste0("====================================================================================================================== \n",
               "== GMSS: spike-and-slab graphical model with spike-and-slab prior for the node-level auxiliary variable coefficients == \n",
               "====================================================================================================================== \n\n"))

    ans <- gmss_vbem_core(Y,
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

  return(ans)
}

