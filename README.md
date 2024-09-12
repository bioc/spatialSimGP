# spatialSimGP

## Introduction

`spatialSimGP` is a simulation tool that generates spatial
transcriptomics data. The purpose of this package is to use a Gaussian
Process for each gene to simulate data with spatial variation. We use
the Poisson distribution to simulate the values on the raw counts scale.
The mean and variance are tied together in the Poisson distribution, so
we simulate the mean-variance relationship with our function. The
mean-variance relationship is a bias in real spatial transcriptomics
data, so we must make sure it is a feature of in silico data as well.
`spatialSimGP` provides the option to simulate data with a fixed or
unique length scale for each gene. The simulated data can be used to
evaluate the performance of spatial transcriptomics analysis methods.

Bioconductor houses the infrastructure to store and analyze spatially
resolved transcriptomics data for R users, including many SVG detection
methods. This simulation framework can be used to benchmark SVG
detection methods and to develop new methods for spatially resolved
transcriptomics data. Additionally, this package interfaces with the
widely used `SpatialExperiment` class from Bioconductor.

## Installation

The following code will install the latest release version of the
`spatialSimGP` package from Bioconductor. Additional details are shown
on the [Bioconductor](https://bioconductor.org/packages/spatialSimGP)
page.

```{r, eval=FALSE}
install.packages("BiocManager")
BiocManager::install("spatialSimGP")
```

The latest development version can also be installed from the `devel`
version of Bioconductor or from
[GitHub](https://github.com/kinnaryshah/spatialSimGP).

## Simulation Framework

The simulation framework is as follows:

$$\boldsymbol{c(s)}|\lambda(\boldsymbol{s}) \sim Poisson (\lambda(\boldsymbol{s})); \lambda(\boldsymbol{s})= exp(\boldsymbol{\beta} + \boldsymbol{C}(\sigma^2))$$

-   $\boldsymbol{s}$: spatial locations
-   $\boldsymbol{\beta}$: vector of means per gene
-   $\sigma^2$: spatial component of variance
-   $\boldsymbol{C}$: covariance function using a Matern kernel with
    squared exponential distance

The exponential covariance function is as follows:

$$(C_{ij}(\boldsymbol{\theta})) = \sigma^2\exp(\frac{-||\boldsymbol{s_i}-\boldsymbol{s_j}||}{l})$$

-   $\boldsymbol{\theta} = (\sigma^2, l)$
-   $l$: length scale parameter
    -   sets how quickly spatial correlation decays with distance
-   $||\boldsymbol{s_i}-\boldsymbol{s_j}||$: Euclidean distance between
    spatial locations

We calculate the covariance matrix using the exponential covariance
function. Using mean $\boldsymbol{0}$ and covariance
$C(\boldsymbol{\theta})$ in the multivariate Normal distribution, we
simulate a Gaussian Process per gene. We use the Gaussian process and
$\beta$ to calculate $\lambda$ and then use the Poisson distribution to
simulate the gene expression levels for each spot.

## Tutorial

**Load packages and data**

```{r}
library(MASS)
library(SpatialExperiment)
library(STexampleData)
library(ggplot2)
```

**Simulating Data with Prior Coordinates Matrix**

One way to simulate data is to provide a matrix of coordinates. In this
example, we use a subset of spots from
`STexampleData::Visium_humanDLPFC()`, which is available from
Bioconductor.

```{r}
spe_demo <- Visium_humanDLPFC()

colData(spe_demo)$subset <- ifelse(colData(spe_demo)$array_row > 20 & colData(spe_demo)$array_row < 65 & colData(spe_demo)$array_col > 20 & colData(spe_demo)$array_col < 65, 1, 0)
spe_demo <- spe_demo[, colData(spe_demo)$subset == 1]

coords <- spatialCoords(spe_demo)
```

We also have to define our remaining parameters before simulating the
data.

-   `n_genes` is the total number of genes to simulate. In this example,
    we simulate 10 genes.
-   `proportion` is the proportion of genes that will have no spatially
    varying patterns. In other words, these genes will just have random
    noise. In this example, 50% of the genes will have no spatial
    patterns.
-   `range_sigma.sq` is the range of the spatial variance parameter. In
    this example, the spatial variance parameter will range from 0.2 to
    3.  
-   `range_beta` is the range of the mean expression value. In this
    example, the mean parameter will range from 0.5 to 9.

```{r}
n_genes <- 5
proportion <- 0.4
range_sigma.sq <- c(1.5, 3)
range_beta <- c(3, 7)
```

**(A) Simulating Data with Fixed Length Scale**

We first simulate 5 genes with a fixed length scale parameter. The
length scale parameter determines how quickly the correlation decays
with distance. Larger length scale parameters simulate larger spatial
patterns. The `simulate` function returns a `SpatialExperiment` object
with the simulated data. Remember to set the seed for reproducibility.

```{r}
length_scale <- 60

set.seed(16)
spe <- spatial_simulate(n_genes, proportion, coords, range_sigma.sq, range_beta, length_scale, length_scale_option = "fixed")
```


**(B) Simulating Data with Unique Length Scale**

We can also simulate data with a unique length scale for each gene. This
process is slower than simulating data with a fixed length scale, but it
allows for more flexibility in the spatial patterns of each gene. Each
gene has a unique length scale parameter, so the Gaussian Process kernel
must be calculated for each gene, slowing down the simulation process.

```{r}
length_scale <- c(60, 40, 20, 80, 100)

set.seed(1)
spe <- spatial_simulate(n_genes, proportion, coords, range_sigma.sq, range_beta, length_scale, length_scale_option = "unique")
```

**Simulating Data with User-Created Coordinates Matrix**

If you have your own coordinates matrix, you can use that to simulate
data. We have included an example below.

```{r}
# 10 spots per side
n_side <- 20

# x and y coordinates for the grid
x_coords <- rep(1:n_side, each = n_side)
y_coords <- rep(1:n_side, times = n_side)

# combine into a matrix
coords <- cbind(x_coords, y_coords)
colnames(coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

# run the simulation
set.seed(1)
length_scale <- 60
spe <- spatial_simulate(n_genes, proportion, coords, range_sigma.sq, range_beta, length_scale, length_scale_option = "fixed")
```

Note: If you want to have complete control over each simulated gene, you
can set `n_genes` = 1, `proportion` = 0, `range_sigma.sq` = c(a,a), and
`range_beta` = c(b,b). This will allow you to simulate one gene at a
time at the exact spatial variance and mean expression level desired.
You could loop through this process to simulate multiple genes with
different parameters.

```{r}
set.seed(123) 
n_genes <- 1 
proportion <- 0 
range_sigma.sq <- c(1, 1)
range_beta <- c(3, 3)
length_scale <- 60

spe <- spatial_simulate(n_genes, proportion, coords, range_sigma.sq, range_beta, length_scale, length_scale_option = "fixed")
```

## Session Info

```{r}
sessionInfo()
```
