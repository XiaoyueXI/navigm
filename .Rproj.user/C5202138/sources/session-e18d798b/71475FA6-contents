rm(list = ls())

set.seed(123)

############################
## simulate basic dataset ##
############################

#
N <- 100
P <- 20
Q <- 10
Q0 <- 1

#
V <- generate_V(P, Q, alpha = 0.05, beta = 0.2, Sigma = diag(Q), min_gene = 1)

#
beta0 <- 0.5
sig2_beta0 <- 0.1
beta_true <- rlnorm(Q0, log(beta0), sig2_beta0)
beta_true_gmss <- rep(0, Q)
beta_true_gmss[sample(1:Q, Q0)] <- beta_true

#
theta <- V %*% matrix(beta_true_gmss, ncol = 1)
zeta <-  - 1.2
pe <- matrix(theta, nrow = P, ncol = P)
pe <- pe + t(pe) + zeta
pe <- pnorm(pe)
A <- 0 + (pe >= 0.5)
diag(A) <- 0

#
net <-  generate_data_from_adjancency(N = N, A = A)

########################
## navigm inference ##
########################


gmss_vbem <- navigm(net$Y, V, numCores = 2,
                    list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))
gmn_vbem <- navigm(net$Y, V, method = 'GMN', numCores = 2,
                   list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))
gm_v1_vbem <- navigm(net$Y, V, method = 'GM', version = 1, numCores = 2,
                     list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))
gm_v2_vbem <- navigm(net$Y, V, method = 'GM', version = 2, numCores = 2,
                     list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))

gmss_em <- navigm(net$Y, V, inference = 'EM',numCores = 2,
                  list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))
gmn_em <- navigm(net$Y, V, method = 'GMN', inference = 'EM', numCores = 2,
                 list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))
gm_v1_em <- navigm(net$Y, V, method = 'GM', version = 1, inference = 'EM', numCores = 2,
                   list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))
gm_v2_em <- navigm(net$Y, V, method = 'GM', version = 2, inference = 'EM', numCores = 2,
                   list_hyper = list(v0_v = seq(1e-4, 1, length.out = 16)))

