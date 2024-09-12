#' Simulate SpatialExperiment object
#'
#' Simulate a SpatialExperiment object with spatially varying genes
#'
#' @details This function simulates a SpatialExperiment object with spatially varying
#'  genes. The function takes in the number of genes to simulate, the proportion of
#'  genes that will have no spatially varying patterns, a matrix of coordinates, the
#'  range of the spatial variance parameter, the range of the mean expression value,
#'  the length scale parameter, and the length scale option.

#' @param n_genes an integer specifying the number of genes to simulate.
#' @param proportion a numeric value specifying the proportion of genes that will have
#'  no spatially varying patterns.
#' @param coords a matrix of coordinates.
#' @param range_sigma.sq a numeric vector of length 2 specifying the range of the spatial
#'  variance parameter.
#' @param range_beta a numeric vector of length 2 specifying the range of the mean
#'  expression value.
#' @param length_scale if length_scale_option is "fixed", a numeric value specifying the
#'  length scale parameter. If length_scale_option is "unique", a numeric vector of length
#'  n_genes specifying the length scale parameter for each gene.
#' @param length_scale_option a character string specifying the length scale option.
#'  Options are "fixed" for a single length scale for all genes or "unique" for a unique
#'  length scale for each gene.
#'
#' @return A SpatialExperiment object with the simulated data.
#'
#' @import MASS
#' @import SpatialExperiment
#' @export
#'
#' @examples
#' library(STexampleData)
#'
#' set.seed(1)
#' n_genes <- 1
#' proportion <- 0.5
#' range_sigma.sq <- c(0.2, 3)
#' range_beta <- c(0.5, 9)
#' length_scale <- 60
#'
#' spe_demo <- Visium_mouseCoronal()
#' colData(spe_demo)$subset <- ifelse(
#'   colData(spe_demo)$array_row > 20 &
#'   colData(spe_demo)$array_row < 65 &
#'   colData(spe_demo)$array_col > 30 &
#'   colData(spe_demo)$array_col < 65,
#'   TRUE, FALSE
#' )
#' spe_demo <- spe_demo[, colData(spe_demo)$subset]
#' coords <- spatialCoords(spe_demo)
#'
#' spe <- spatial_simulate(n_genes, proportion, coords, range_sigma.sq, range_beta,
#'    length_scale, length_scale_option = "fixed")
#'
spatial_simulate <- function(n_genes, proportion, coords,
                             range_sigma.sq, range_beta,
                             length_scale, length_scale_option = "fixed") {

  # check the proportion is between 0 to 1
  if (proportion < 0 || proportion > 1) {
    stop("proportion must be between 0 and 1")
  }

  # check if coords is a matrix
  if (!is.matrix(coords)) {
    stop("coords must be a matrix")
  }

  # check range_sigma.sq is a numeric vector of length 2
  if (!is.numeric(range_sigma.sq) || length(range_sigma.sq) != 2) {
    stop("range_sigma.sq must be a numeric vector of length 2")
  }

  # check range_beta is a numeric vector of length 2
  if (!is.numeric(range_beta) || length(range_beta) != 2) {
    stop("range_beta must be a numeric vector of length 2")
  }


  # check if length_scale_option is "fixed", then length_scale must be a single number
  if (length_scale_option == "fixed" && length(length_scale) != 1) {
    stop("length_scale must be a single number when using fixed length scale")
  }

  if (length_scale_option == "fixed") {
    length_scale_single <- length_scale
  } else if (length_scale_option == "unique") {
    if (length(length_scale) != n_genes) {
      stop("length_scale vector must match the number of genes when using unique length scales")
    }
  } else {
    stop("invalid length_scale_option. use 'fixed' or 'unique'")
  }

  #some genes have some nonzero sigma.sq
  sigma.sq <- runif(n_genes, range_sigma.sq[1], range_sigma.sq[2])
  #some genes have zero sigma.sq
  sigma.sq[sample(seq_len(n_genes), proportion*n_genes)] <- 0
  ground_truth_rank <- rank(-sigma.sq)

  #all genes have nonzero beta values
  beta <- runif(n_genes, log(range_beta[1]), log(range_beta[2]))

  params <- data.frame(sigma.sq, beta)

  points_coord <- coords
  n_points <- nrow(points_coord)

  rownames(points_coord) <- NULL
  points_coord <- as.matrix(points_coord)

  # step 1: create pairs of points to calculate Euclidean distances

  pair.points <- cbind(
    matrix( rep(points_coord, each = n_points), ncol = 2, byrow = FALSE),
    rep(1, times = n_points) %x% points_coord # Creating the combinations using kronecker product
  ) |> data.frame()

  colnames(pair.points) <- c("si.x", "si.y", "sj.x", "sj.y")

  #step 2: calculate Gaussian process/kernel

  kernel.fun <- function(si.x, si.y, sj.x, sj.y,  l = 0.2){
    exp(-1*sqrt(((si.x-sj.x)^2+(si.y-sj.y)^2))/l)
  }

  if (length_scale_option == "fixed") {
    C_theta <- with(pair.points, kernel.fun(si.x, si.y, sj.x, sj.y, l = length_scale_single)) |>
      matrix(nrow = n_points, ncol = n_points)
  }

  counts <- matrix(NA, nrow = n_genes, ncol = n_points)

  for (i in seq_len(n_genes)) {

    message(sprintf("Simulating gene %d", i))
    sigma.sq_i <- sigma.sq[i]
    beta_i <- beta[i]

    if (length_scale_option == "unique") {
      C_theta <- with(pair.points, kernel.fun(si.x, si.y, sj.x, sj.y, l = length_scale[i])) |>
        matrix(nrow = n_points, ncol = n_points)
    }

    #step 3: simulate Gaussian process per gene

    gp_dat <- mvrnorm(n = 1, rep(0,n_points), sigma.sq_i* C_theta)

    #step 4: calculate lambda = exp(beta + Gaussian process) per gene

    lambda_i <- exp(gp_dat + beta_i)

    #step 5: simulate values per gene using Poisson distribution

    counts_i <- rpois(n = n_points, lambda_i)

    #put all counts in matrix
    #orientation: genes x spots

    counts[i,] <- counts_i
  }

  #create SpatialExperiment object using counts and coords matrices

  spe <- SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = points_coord)

  rowData(spe)$ground_truth_rank <- ground_truth_rank
  rowData(spe)$ground_truth_sigma.sq <- sigma.sq
  rowData(spe)$ground_truth_beta <- beta

  return(spe)

}
