
```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# CatchAll


CatchAll: an R package to estimate species richness with finite exponential-mixed Poisson model.



This package contains functions that estimate the species richness with finite exponential-mixed Poisson models (up to the three-component exponential)
mixed-Poisson model).
## Installation ##


```r
library(devtools)
install_github("/statdivlab/CatchMore/CatchAll")
library(CatchAll)
```

## Basic Usage ##


```r
library(breakaway)
data(apples)
```

Poisson Model:
```r
Poisson_model(apples, cutoff = 20)
```

Geometric model (or single exponential-mixed Poisson model):
```r
geometric_model(apples, cutoff = 20)
```

Two- or three-component exponential mixture mixed Poisson model:
```r
two_geometric_model(apples, cutoff = 50)
three_geometric_model(apples, cutoff = 50)
```

An integrated function for simutaneously viewing the results of all 
above models for different cut-off of species frequency counts.

```r 
apples_results <- all_parametric_model(apples)
head(apples_results, 20)
```
