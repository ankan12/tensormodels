
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tensor_models

Tensor operations, decompositions, and distributions in R.

## Installation

You can install the development version of tensormodels from
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

# Simulating Tensor Variate Normal Draws

The function `rtnorm` simulates from the tensor variate normal with a
specified mean array mu and a list of covariance matrices called
list_sigmas.

Since the univariate normal and matrix variate normal are simpler cases
of the tensor variate normal, this function can simulate from them as
well. Here, we simulate from a univariate normal with mean $-2$ and
variance $4$.

``` r
univar_draws <- rtnorm(n = 10000, mu = -2, list_sigmas = 4)
mean(univar_draws)
#> [1] -1.984154
sd(univar_draws)
#> [1] 2.007343
```

And now we simulate from the matrix variate normal of size $2 \times 3$.
By default, the covariance matrices will be the identity.

``` r
matrix_draws <- rtnorm(n = 10000, mu = matrix(1:6, nrow = 2, ncol = 3))

matrix_draws[1, , ]
#>          [,1]     [,2]     [,3]
#> [1,] 1.570040 3.090437 3.874926
#> [2,] 1.756116 3.012314 8.072248
```

Below is a simulation from a tensor variate normal of size
$3 \times 4 \times 2$.

``` r
tensor_draws <- rtnorm(n = 1e3, mu = array(1:24, dim = c(3, 4, 2)))

tensor_draws[1, , , ]
#> , , 1
#> 
#>          [,1]     [,2]     [,3]     [,4]
#> [1,] 1.636847 2.465285 6.012573 11.69208
#> [2,] 2.363168 4.173519 6.918542 10.61718
#> [3,] 1.577817 5.582948 8.627753 12.00810
#> 
#> , , 2
#> 
#>          [,1]     [,2]     [,3]     [,4]
#> [1,] 12.89172 15.57493 18.23921 22.75339
#> [2,] 14.21812 18.37944 21.06593 22.43802
#> [3,] 14.67745 18.00590 20.13545 23.00758
```

# MLE Estimation for the Tensor Variate Normal

The package supports MLE estimation for the mean and list of covariance
matrices.

``` r
univar_est <- mle_est(univar_draws)

univar_est$mu
#> [1] -1.984154
univar_est$sigmas
#> [[1]]
#>          [,1]
#> [1,] 4.029426
```

Here, the mean matrix from the draws above has values 1 through 6. The
covariance matrices are the identity matrices.

``` r
matrix_est <- mle_est(matrix_draws)
#> Iteration 1: max relative change = 1.228e-02
#> Iteration 2: max relative change = 9.964e-05
#> Iteration 3: max relative change = 4.082e-08
#> Converged at iteration 3

matrix_est$mu
#>           [,1]     [,2]     [,3]
#> [1,] 0.9889842 2.993124 5.002082
#> [2,] 1.9874096 3.980462 5.988387
matrix_est$sigmas
#> [[1]]
#>             [,1]        [,2]
#> [1,]  1.00308406 -0.01180287
#> [2,] -0.01180287  0.99706430
#> 
#> [[2]]
#>             [,1]         [,2]         [,3]
#> [1,] 1.007360206  0.007453601  0.008568699
#> [2,] 0.007453601  0.990977958 -0.002157685
#> [3,] 0.008568699 -0.002157685  1.001864848
```

This function works on array objects of dimensions above 2.

``` r
tensor_est <- mle_est(tensor_draws)
#> Iteration 1: max relative change = 3.653e-02
#> Iteration 2: max relative change = 7.607e-04
#> Iteration 3: max relative change = 2.000e-06
#> Iteration 4: max relative change = 9.828e-09
#> Converged at iteration 4

tensor_est$mu
#> , , 1
#> 
#>           [,1]     [,2]     [,3]      [,4]
#> [1,] 0.9950717 4.004322 6.964651  9.999908
#> [2,] 1.9865597 4.991675 7.935281 10.977554
#> [3,] 3.0780787 6.037056 9.005040 12.045153
#> 
#> , , 2
#> 
#>          [,1]     [,2]     [,3]     [,4]
#> [1,] 12.98197 16.04675 19.03532 22.01601
#> [2,] 13.99639 16.97609 19.99714 23.00097
#> [3,] 14.99562 18.02654 20.98631 24.05676
tensor_est$sigmas
#> [[1]]
#>              [,1]        [,2]        [,3]
#> [1,]  0.970626129 0.006219061 -0.01527996
#> [2,]  0.006219061 1.015100266  0.02040243
#> [3,] -0.015279957 0.020402435  1.01563140
#> 
#> [[2]]
#>              [,1]       [,2]         [,3]         [,4]
#> [1,]  1.037256620 0.01974396 -0.005775892  0.002004411
#> [2,]  0.019743960 0.97462151  0.019342349  0.015668624
#> [3,] -0.005775892 0.01934235  1.002474777 -0.023123390
#> [4,]  0.002004411 0.01566862 -0.023123390  0.988340795
#> 
#> [[3]]
#>              [,1]         [,2]
#> [1,]  1.000059904 -0.006475362
#> [2,] -0.006475362  0.999982027
```
