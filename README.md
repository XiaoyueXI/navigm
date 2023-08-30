
# navigm

The `navigm` is an R package that implements a Gaussian graphical
spike-and-slab model for estimating associations between variables
(nodes). In addition to the primary measurements of node variables, this
package allows for incorporating node-level auxiliary variables. The
method is tailored to enhance the detection of conditional dependence
structure by leveraging and selecting useful node-level information. The
package employs a scalable variational Bayes expectation conditional
maximisation algorithm to obtain posterior mode of the precision matrix
that encodes conditional dependence between nodes and approximated
posterior distributions of the remaining model parameters.

The package offers several modeling options, including the choice of
incorporating node-level variables and the option to select influential
variables. Furthermore, it provides an alternative deterministic
inference algorithm, full expectation conditional maximisation, which
has been widely used in previous studies.

Reference: Xiaoyue Xi, Helene Ruffieux, 2023. A hierarchical framework
for inferring and leveraging node-level information in Bayesian
networks.

## Installation

The package can be installed in R by

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigm")
```

For further instructions, please check the (tutorial)\[link\].
