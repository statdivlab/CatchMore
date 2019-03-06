
CatchMore <img src="docs/catchmore.png" align="right" width="165px"/>
=====================================================================

CatchMore is an R package to estimate species richness with mixed-Poisson and non-mixed-Poisson models. It extends the historically popular methodology CatchAll to include a more diverse catalogue of models, and provides an easy-to-use R implementation that runs on Unix machines (unlike CatchAll).

CatchMore is written by [Kendrick Qijin Li](http://students.washington.edu/qijunl2/) and [Amy Willis](http://statisticaldiversitylab.com/team/amy-willis), and is maintained by the [Statistical Diversity Lab](http://statisticaldiversitylab.com/) at the University of Washington.

Installation
------------

``` r
library(devtools)
install_github("statdivlab/CatchMore")
library(CatchMore)
```

Basic Usage
-----------

``` r
library(breakaway)
data(apples)
```

Poisson Model:

``` r
Poisson_model(apples, cutoff = 20)
```

Geometric model (or single exponential-mixed Poisson model):

``` r
geometric_model(apples, cutoff = 20)
```

Two- or three-component exponential mixture mixed Poisson model:

``` r
two_geometric_model(apples, cutoff = 50)
three_geometric_model(apples, cutoff = 70)
```

An integrated function for viewing the results of "best" models.

``` r
apples_results <- catch_all(apples)
apple_results
```

### Citing CatchMore

Writing software is time-consuming and often thankless. If you use this software, please cite our work in the following way:

Li, Q.K. and Willis, A.D. (2019). CatchMore: an R package for richness estimation. In Preparation.

`CatchMore` implements a number of different richness estimates. Please also cite the following if you use them:

-   `breakaway()` and `kemp()`: Willis, A. & Bunge, J. (2015). Estimating diversity via frequency ratios. Biometrics.
-   `catch_all()`: Bunge, J. (2011). Estimating the number of species with CatchAll. Proceedings of the Pacific Symposium on Biocomputing.
