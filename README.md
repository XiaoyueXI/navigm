
# navigss

The `navigss` is an R package that implements the graphical
spike-and-slab in Gaussian graphical models. Besides the primary
measurements of variables, this package allows for incorporating
node-level auxiliary variables. The method is tailored to improve the
detection of conditional independence by leveraging node-level
information while selecting useful node-level variables that influence
the graph structure. The package implements a scalable variational Bayes
expectation-maximization algorithm to obtain the approximated posterior
distribution of the model parameters and the posterior mode of precision
matrices.

The package offers several modeling options, including whether to
incorporate node-level variables and whether to make the selection of
variables. Additionally, it provides deterministic inference options,
such as full expectation maximization and variational Bayes
expectation-maximization algorithms.

Reference: Xiaoyue Xi, Sylvia Richardson, Helene Ruffieux, 2023. A
hierarchical framework for inferring and leveraging node-level
information in Bayesian networks.

## Installation

The package can be installed in R by

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigss")
```

## Example

The main function `navigss` takes a data matrix of dimension N x P,
where N is the sample size and P is the graph size. It allows for
another input matrix of dimension P x Q, which comprises node-level
auxiliary variables and may inform the centrality of nodes, where Q is
the number of auxiliary variables. We illustrate three modeling choices
as follows. In each model, we provide two choices of inference
algorithms.

### Graphical spike-and-slab

We start with a simple example when no node-level variables are
provided.

#### Data simulation

We simulate a random network of size $P = 50$, of which edges are
included with probability 0.01. We then simulate the precision matrix
following the adjacency structure. $N = 100$ samples are then generated
from a multivariate normal distribution with mean **0** and the
simulated precision matrix.

``` r
set.seed(123)

P <- 50
A <- matrix(rbinom(P ^ 2, 1, 0.01), nrow = P, ncol = P)
diag(A) <- 0
A[lower.tri(A)] <- t(A[upper.tri(A)])

N <- 100
net <-  generate_data_from_adjancency(N = N, A = A)

print(str(net$Y))
```

    ##  num [1:100, 1:50] 1.84 -0.539 0.371 -0.459 -1.048 ...
    ##  - attr(*, "scaled:center")= num [1:50] 0.0181 0.0234 -0.2788 0.0782 -0.1783 ...
    ##  - attr(*, "scaled:scale")= num [1:50] 1.16 1.11 1.19 1.16 1 ...
    ## NULL

#### Inference using expectation maximization

We then input the simulated data `net$Y` into the `navigss` and infer
the precision matrix by spike-and-slab Gaussian graphical models using
expectation maximization.

``` r
gm_em <- navigss(Y = net$Y, method = 'GM', inference = 'EM', version = 1)
```

    ## == Parallel exploration of a  grid of spike standard deviations 
    ##  v0 =  1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1 
    ##   on  8  cores ... 
    ## 
    ## ... done. == 
    ## 
    ## == Select from a grid of spike variance using AIC  ...

    ## Warning in navigss(Y = net$Y, method = "GM", inference = "EM", version = 1):
    ## The selected v0 reaches the lower bound in the grid. Consider extend the grid
    ## to lower values.

    ## ... done. == 
    ## 
    ## Select the index 1  i.e., v0 =  1e-04 , the best AIC  =  2622.522 .
    ## 
    ## Total runtime:  2.894116 secs .

#### Inference using variational Bayes expectation maximization

The following command uses the variational expectation maximization
algorithm for graphical spike-and-slab inference.

``` r
gm_vbem <- navigss(Y = net$Y, method = 'GM', version = 1)
```

    ## == Parallel exploration of a  grid of spike standard deviations 
    ##  v0 =  1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1 
    ##   on  8  cores ... 
    ## 
    ## ... done. == 
    ## 
    ## == Select from a grid of spike variance using AIC  ... 
    ## 
    ## ... done. == 
    ## 
    ## Select the index 9  i.e., v0 =  0.53338 , the best AIC  =  4726.792 .
    ## 
    ## Total runtime:  3.09415 secs .

### Graphical spike-and-slab with normal priors on auxiliary variable coefficients

We next simulate a network triggered by three node-level auxiliary
variables.

#### Data simulation

We generate right-skewed beta-distributed node-level variables **V**
(for instance, posterior inclusion probabilities from previous studies),
and their effect size from a log-normal distribution centered at 0.5
$\beta$, which gives rise to hub propensities **$V\beta$**. Set
$\zeta=-1.5$ and we obtain the probability of including the edge
$(i,j)$, i.e., $\zeta + \bm{v_{i} \beta} + \bm{v_{j} \beta}$. $N = 100$
samples are then generated from a multivariate normal distribution with
mean **0** and the adjacency-guided precision matrix.

``` r
set.seed(123)

Q <- 3
Q0 <- 3

#
V <- generate_V(P, Q, alpha = 0.05, beta = 0.2, Sigma = diag(Q), min_gene = 5) 
```

    ## Range of empirical correlations: -0.1386623 0.1388341 
    ## Range of absolute empirical correlations: 0.01920697 0.1388341

``` r
#
beta0 <- 0.5
sig2_beta0 <- 0.1
beta_true <- rlnorm(Q0, log(beta0), sig2_beta0)

#
theta <- V %*% matrix(beta_true, ncol = 1)
pe <- matrix(theta, nrow = P, ncol = P)
zeta <- -1.5
pe <- pe + t(pe) + zeta
pe <- pnorm(pe)
A <- 0 + (pe >= 0.5)
diag(A) <- 0

#
net <-  generate_data_from_adjancency(N = N, A = A)

print(str(net$Y))
```

    ##  num [1:100, 1:50] 0.673 -1.701 0.265 0.173 -0.698 ...
    ##  - attr(*, "scaled:center")= num [1:50] 0.00521 -0.00102 -0.01538 0.02928 -0.03965 ...
    ##  - attr(*, "scaled:scale")= num [1:50] 0.582 0.697 0.644 0.633 0.74 ...
    ## NULL

#### Inference using expectation maximization

We then input the simulated data `net$Y`and `V` into the `navigss` and
infer the precision matrix by spike-and-slab Gaussian graphical models
with normal priors on auxiliary variable coefficients using expectation
maximization.

``` r
gmn_em <- navigss(Y = net$Y, V = V, method = 'GMN', inference = 'EM')
```

    ## == Parallel exploration of a  grid of spike standard deviations 
    ##  v0 =  1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1 
    ##   on  8  cores ... 
    ## 
    ## ... done. == 
    ## 
    ## == Select from a grid of spike variance using AIC  ...

    ## Warning in navigss(Y = net$Y, V = V, method = "GMN", inference = "EM"): The
    ## selected v0 reaches the lower bound in the grid. Consider extend the grid to
    ## lower values.

    ## ... done. == 
    ## 
    ## Select the index 1  i.e., v0 =  1e-04 , the best AIC  =  2544.074 .
    ## 
    ## Total runtime:  9.021525 secs .

#### Inference using variational Bayes expectation maximization

To use the variational expectation maximization, we implement the
command as follows.

``` r
gmn_vbem <- navigss(Y = net$Y, V = V, method = 'GMN')
```

    ## == Parallel exploration of a  grid of spike standard deviations 
    ##  v0 =  1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1 
    ##   on  8  cores ... 
    ## 
    ## ... done. == 
    ## 
    ## == Select from a grid of spike variance using AIC  ... 
    ## 
    ## ... done. == 
    ## 
    ## Select the index 6  i.e., v0 =  0.3334 , the best AIC  =  4405.267 .
    ## 
    ## Total runtime:  4.110883 secs .

### Graphical spike-and-slab with spike-and-slab priors on auxiliary variable coefficients

We next consider a scenario where a large number of candidate node-level
auxiliary variables are available, and three of them contribute to the
node centrality. The selection of “active” node-level auxiliary
variables is thus of interest.

#### Data simulation

We generate the node-level variable matrix of dimension P x Q where
$P = 100$ and $Q = 50$, and set most of the effect sizes to zeros, of
which three non-null effect sizes are generated from log normal
distribution centered at 0.5 **$\beta$**. These gives rise to hub
propensities **$V\beta$**. Set $\zeta=-1.5$, we obtain the probability
of including the edge $(i,j)$, i.e.,
$\zeta + \bm{v_{i} \beta} + \bm{v_{j} \beta}$. $N = 100$ samples are
then generated from a multivariate normal distribution with mean **0**
and the adjacency-guided precision matrix.

``` r
set.seed(123)

Q <- 50 
Q0 <- 3

#
V <- generate_V(P, Q, alpha = 0.05, beta = 0.2, Sigma = diag(Q), min_gene = 5) 
```

    ## Range of empirical correlations: -0.2958226 0.4269854 
    ## Range of absolute empirical correlations: 0.0001344982 0.4269854

``` r
beta_true_gmss <- rep(0, Q)
beta_true_gmss[sample(1:Q, Q0)] <- beta_true

#
theta <- V %*% matrix(beta_true_gmss, ncol = 1)
pe <- matrix(theta, nrow = P, ncol = P)
pe <- pe + t(pe) + zeta
pe <- pnorm(pe)
A <- 0 + (pe >= 0.5)
diag(A) <- 0

#
net <-  generate_data_from_adjancency(N = N, A = A)

print(str(net$Y))
```

    ##  num [1:100, 1:50] -0.109 -0.215 1.888 1.099 0.649 ...
    ##  - attr(*, "scaled:center")= num [1:50] 0.0943 -0.0335 -0.0197 -0.1027 -0.0152 ...
    ##  - attr(*, "scaled:scale")= num [1:50] 0.694 0.592 0.533 0.616 0.589 ...
    ## NULL

#### Inference using expectation maximization

We then input the simulated data `net$Y`and `V` into the `navigss` and
infer the precision matrix by spike-and-slab Gaussian graphical models
with spike-and-slab priors on auxiliary variable coefficients using
expectation maximization.

``` r
gmss_em <- navigss(Y = net$Y, V = V, method = 'GMSS', inference = 'EM')
```

    ## == Parallel exploration of a double grid of spike standard deviations 
    ##  v0 =  1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1 
    ##  and s0 = 1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1
    ##  on  8  cores ... 
    ## 
    ## ... done. == 
    ## 
    ## == Select from a grid of spike variance using AIC  ...

    ## Warning in navigss(Y = net$Y, V = V, method = "GMSS", inference = "EM"): The
    ## selected v0 reaches the lower bound in the grid. Consider extend the grid to
    ## lower values.

    ## ... done. == 
    ## 
    ## Select the index 1 8  i.e., v0 =  1e-04  and s0 =  0.46672 ,the best AIC  =  2143.57 .
    ## 
    ## Total runtime:  2.460255 mins .

#### Inference using variational Bayes expectation maximization

To use the variational expectation maximization algorithm, we implement
the command

``` r
gmss_vbem <- navigss(Y = net$Y, V = V, method = 'GMSS')
```

    ## == Parallel exploration of a  grid of spike standard deviations 
    ##  v0 =  1e-04 0.06676 0.13342 0.20008 0.26674 0.3334 0.40006 0.46672 0.53338 0.60004 0.6667 0.73336 0.80002 0.86668 0.93334 1 
    ##   on  8  cores ... 
    ## 
    ## ... done. == 
    ## 
    ## == Select from a grid of spike variance using AIC  ... 
    ## 
    ## ... done. == 
    ## 
    ## Select the index 3  i.e., v0 =  0.13342 , the best AIC  =  4241.623 .
    ## 
    ## Total runtime:  11.27738 secs .

This shows less computational time for the variational expectation
maximization algorithm, especially when we leverage and select useful
node-level auxiliary variables, as there is no need for double grid
search in two-level continuous spike-and-slabs as in the expectation
maximization.
