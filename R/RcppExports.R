# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

weightedGaussianLocalDensity <- function(weights, distance, nrow, dc) {
    .Call('_weightedClustSuite_weightedGaussianLocalDensity', PACKAGE = 'weightedClustSuite', weights, distance, nrow, dc)
}

weightedNonGaussianLocalDensity <- function(weights, truesize, distance, nrow, dc) {
    .Call('_weightedClustSuite_weightedNonGaussianLocalDensity', PACKAGE = 'weightedClustSuite', weights, truesize, distance, nrow, dc)
}

SumCutOff <- function(weights, distance, nrow, dc) {
    .Call('_weightedClustSuite_SumCutOff', PACKAGE = 'weightedClustSuite', weights, distance, nrow, dc)
}

extractRowColDistfromDistMatrix <- function(distance, num_row, row_inds, col_inds) {
    .Call('_weightedClustSuite_extractRowColDistfromDistMatrix', PACKAGE = 'weightedClustSuite', distance, num_row, row_inds, col_inds)
}

DensityPeak_smallest_dist_rho_order_coords <- function(ordered_rho, ordered_coords) {
    .Call('_weightedClustSuite_DensityPeak_smallest_dist_rho_order_coords', PACKAGE = 'weightedClustSuite', ordered_rho, ordered_coords)
}

DensityPeak_distancesToPeakC <- function(distance, rho) {
    .Call('_weightedClustSuite_DensityPeak_distancesToPeakC', PACKAGE = 'weightedClustSuite', distance, rho)
}

createWeightedDist <- function(distance, num_row, weights) {
    .Call('_weightedClustSuite_createWeightedDist', PACKAGE = 'weightedClustSuite', distance, num_row, weights)
}

