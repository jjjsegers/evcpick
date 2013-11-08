evcpick
=======

This is a small R package providing functions for rank-based versions of three estimators for the Pickands dependence function of a bivariate extreme value copula: the Pickands estimator, the Caperaa-Fougeres-Genest estimator, and the Buecher-Dette-Volgushev estimator.

Installation
------------

To install the package in R, you must first install the `devtools` package:

```
install.package("devtools")
```
  
Then, to install the `effcopgauss` package, simply do

```
library(devtools)
install_github("evcpick", "jjjsegers")
```

Usage
-----

If the package is installed, you can do for instance

```
library(evcpick)
? evcPickands
? evcCFG
? evcBDV
```
  
to see the documentation.