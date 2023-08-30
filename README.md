
# navigm

The `navigm` is an R package that implements a Gaussian graphical
spike-and-slab model for estimating associations between variables
(nodes). Alongside the primary node measurements, this package permits
the inclusion of supplementary node-level auxiliary variables. This
method is tailored to enhance the identification of conditional
dependence structure by leveraging these variables and selecting useful
information from them. Employing a scalable variational Bayes
expectation conditional maximisation algorithm, the package attains the
posterior mode of the precision matrix that pinpoints the conditional
dependence between nodes, along with approximated posterior
distributions of the remaining model parameters.

The package offers a variety of modelling options, including the
decision to integrate node-level variables as well as the flexibility to
select influential variables. Furthermore, it provides an alternative
deterministic inference algorithm, full expectation conditional
maximisation, which has been widely used in previous studies.

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
