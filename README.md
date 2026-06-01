
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LIsing

The goal of LIsing is to estimate the co-occurence network for binary
repeated measures data.

## Installation

You can install the development version of LIsing like so:

``` r
library(LIsing)
```

### Network Estimation

The following chunk of code simulates a dataset containing
presence/absence information for 50 subjects taken over 4 time points.

``` r
n <- 200; p <- 10
X <- matrix(data = NA, nrow = n, ncol = p + 3)
colnames(X) <- c("ID", "Measure", "Time", paste0("M", 1:p))

X[,"ID"] <- rep(1:50, each = 4)
X[,"Measure"] <- rep(c(1:4), 50)
X[,"Time"] <- X[,"Measure"]

X[,-c(1:3)] <- rbernoulli(n = n*p, p = 0.5)
#> Warning: `rbernoulli()` was deprecated in purrr 1.0.0.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
head(X)
#>      ID Measure Time M1 M2 M3 M4 M5 M6 M7 M8 M9 M10
#> [1,]  1       1    1  0  1  0  0  0  0  1  0  1   1
#> [2,]  1       2    2  0  0  0  0  1  0  1  0  1   1
#> [3,]  1       3    3  1  1  0  1  0  1  1  0  0   1
#> [4,]  1       4    4  0  0  0  0  1  0  1  1  0   1
#> [5,]  2       1    1  1  0  1  0  0  0  1  0  1   1
#> [6,]  2       2    2  0  0  1  0  1  0  1  0  0   1
```

The following chunk of code estimates the network from this using the
nodewise regression approach.

``` r
theta.res <- estimate.parameters.nw(X = X, theta.init = NULL, tau.init = 0.5, penalty = "1se", symm = "min")
print(theta.res)
#>       [,1]        [,2]        [,3] [,4] [,5] [,6]        [,7] [,8] [,9] [,10]
#>  [1,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
#>  [2,]    0 0.004266741  0.00000000    0    0    0  0.00000000    0    0     0
#>  [3,]    0 0.000000000  0.00000000    0    0    0 -0.04674539    0    0     0
#>  [4,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
#>  [5,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
#>  [6,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
#>  [7,]    0 0.000000000 -0.04674539    0    0    0  0.00000000    0    0     0
#>  [8,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
#>  [9,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
#> [10,]    0 0.000000000  0.00000000    0    0    0  0.00000000    0    0     0
```

The following chunk of code estimates the network from this using the
joint approach.

``` r
theta.res <- estimate.parameters.j(X = X, theta.init = NULL, tau.init = 0.5, penalty = "1se")
print(theta.res)
#>       [,1] [,2]        [,3] [,4] [,5] [,6]        [,7] [,8] [,9] [,10]
#>  [1,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#>  [2,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#>  [3,]    0    0  0.00000000    0    0    0 -0.05220203    0    0     0
#>  [4,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#>  [5,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#>  [6,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#>  [7,]    0    0 -0.05220203    0    0    0  0.00000000    0    0     0
#>  [8,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#>  [9,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
#> [10,]    0    0  0.00000000    0    0    0  0.00000000    0    0     0
```
