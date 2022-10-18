# fastWavelets

Computes Maximal Overlap Discrete Wavelet Transform (MODWT) and A Trous DWT.

## Installation

You can install the development version of fastWavelets as follows:

``` r
# install.packages("devtools")
devtools::install_github("johnswyou/fastWavelets")
```

## Loading

``` r
library(fastWavelets)
```
## Usage

```r
N = 1000                # number of time series points
J = 4                   # decomposition level
wavelet = 'coif1'       # scaling filter
X = matrix(rnorm(N),N,1)
W = mo_dwt(X,wavelet,J)
```

The set of possible value for the argument `wavelet` is as follows:

```r
c("haar", "d1", "sym1", "bior1.1", "rbio1.1",
"d2", "sym2", "d3", "sym3", "d4", "d5", "d6", "d7", "d8", "d9", "d10", "d11",
"sym4", "sym5", "sym6", "sym7", "sym8", "sym9", "sym10",
"coif1", "coif2", "coif3", "coif4", "coif5",
"bior1.3", "bior1.5", "bior2.2", "bior2.4", "bior2.6", "bior2.8", "bior3.1", "bior3.3",
"bior3.5", "bior3.7", "bior3.9", "bior4.4", "bior5.5", "bior6.8",
" rbio1.3", " rbio1.5", " rbio2.2", " rbio2.4", " rbio2.6", " rbio2.8", " rbio3.1", " rbio3.3",
" rbio3.5", " rbio3.7", " rbio3.9", " rbio4.4", " rbio5.5", " rbio6.8",
"la8", "la10", "la12", "la14", "la16", "la18", "la20",
"bl14", "bl18", "bl20",
"fk4", "fk6", "fk8", "fk14", "fk18", "fk22",
"b3spline")
```
