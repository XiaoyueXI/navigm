
# navigm

`navigm` is an R package that implements variable-guided network
inference using Bayesian graphical spike-and-slab modelling. Alongside
the primary node measurements, our framework encodes node-level
auxiliary variables that may be informative on the network structure.
For instance, gene network inference may be informed by the use of
publicly available summary statistics on the regulation of genes by
genetic variants.

Our approach relies on a fully joint hierarchical model to
simultaneously infer (i) sparse precision matrices and (ii) the
relevance of node-level information for uncovering the sought-after
network structure. Inference is carried out using an efficient
variational expectation conditional maximisation algorithm that scales
to hundreds of samples, nodes and auxiliary variables, and approximates
full posterior distributions for parameters of interest.

Reference: Xiaoyue Xi, Hélène Ruffieux, 2023. A modelling framework for
detecting and leveraging node-level information in Bayesian network
inference.

## Installation

The package can be installed in R using the following command:

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("XiaoyueXI/navigm")
```

For further instructions, please check the
[tutorial](https://xiaoyuexi.github.io/navigm.github.io/).
