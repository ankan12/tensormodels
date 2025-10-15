
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tensor_models

Tensor operations, decompositions, and distributions in R.

## Installation

You can install the development version of regexcite from
[GitHub](https://github.com/) with:

``` r
#install.packages("devtools")
devtools::install_github("ankan12/tensormodels")
```

# Operations

# n-mode prod

The n-mode product between a tensor and a matrix can be computed with
`n_mode_prod`.

``` r
a <- array(1:3, dim = c(3, 1, 1))
b <- matrix(4:9, nrow = 2, ncol = 3)
n_mode_prod(a, b, 1)
#> , , 1
#> 
#>      [,1]
#> [1,]   40
#> [2,]   46
```

# nm-mode prod

The nm-mode product between a tensor and a tensor can be computed with
`nm-mode-prod`.

``` r
A <- matrix(c(1, 2, 3, 4), nrow = 2)
b <- matrix(c(5, 6), nrow = 2)
nm_prod(A, b, 1, 1)
#>      [,1]
#> [1,]   17
#> [2,]   39
```

# tensor prod

The tensor product between a tensor and a tensor can be computed with
`tensor_prod`.

``` r
A <- matrix(c(1, 2, 3, 4), nrow = 2)
b <- matrix(c(5, 6), nrow = 2)
tensor_prod(A, b)
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]    5   15
#> [2,]   10   20
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> [1,]    6   18
#> [2,]   12   24
```
