
# navigm

The `navigm` is an R package that implements the graphical
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
remotes::install_github("XiaoyueXI/navigm")
```
