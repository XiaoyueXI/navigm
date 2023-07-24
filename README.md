
## NAVIGSS

The `navigss` is an R package that implements the graphical
spike-and-slab in the Gaussian graphical models. In addition to primary
measurements of variables, this package allows for incorporating
node-level auxiliary variables. The method is tailored to improving the
detection of conditional independence by leveraging node-level
information and, at the same time, selecting useful node-level variables
that influence the graph structure. The package implements a scalable
variational Bayes expectation maximization algorithm for obtaining the
approximated posterior distribution of the model parameters and the
posterior mode of precision matrices.

The package offers several modeling options, including whether to
incorporate node-level variables and whether to make the selection of
variables. It also provides deterministic inference options, including
full expectation maximization and variational Bayes expectation
maximization algorithms.

Reference:

## Installation

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigss")
```
