#' Simulation Functions for emphasis Package
#'
#' This file contains functions for simulating phylogenetic trees and evolutionary processes.
#'
#' @section Functions:
#' - simulate_tree: Simulate phylogenetic trees from scratch under a chosen model.
#' - augment_tree: Augment an existing tree using a user-defined importance sampler.
#' - generateNonHomogeneousExp: Simulate events from a non-homogeneous Poisson process.
#' - nhExpRand: Sample event times from a non-homogeneous exponential process.

#' Simulate a Phylogenetic Tree from Scratch
#'
#' @description
#' Simulate a new phylogenetic tree under a specified diversification model (e.g., birth-death, diversity-dependent).
#' The user must specify the model and its parameters.
#'
#' @param model A list or object specifying the diversification model and its parameters.
#' @param ... Additional arguments (e.g., tree size, time, etc.).
#' @return A simulated phylo object.
#' @examples
#' # Example (pseudo-code):
#' model <- list(type = "birth-death", lambda = 0.1, mu = 0.05)
#' tree <- simulate_tree(model, n_taxa = 10)
#' @export
simulate_tree <- function(model, ...) {
  stop("Not yet implemented. Please see vignettes for simulation examples.")
}

#' Augment an Existing Phylogenetic Tree
#'
#' @description
#' Augment a given phylogenetic tree by simulating missing or unobserved lineages using a user-defined importance sampler.
#' This is useful for likelihood-based inference with incomplete trees or missing taxa.
#'
#' @param tree A phylo object (from ape) representing the observed tree.
#' @param importance_sampler A function or object specifying the importance sampling process for augmentation.
#' @param ... Additional arguments (e.g., number of augmentations).
#' @return An augmented phylo object or a list of augmented trees.
#' @examples
#' # Example (pseudo-code):
#' importance_sampler <- function(tree, ...) { tree } # placeholder
#' augmented <- augment_tree(tree, importance_sampler, n_augments = 10)
#' @export
augment_tree <- function(tree, importance_sampler, ...) {
  stop("Not yet implemented. Please see vignettes for augmentation examples.")
}

#' Generate Non-Homogeneous Exponential Events
#'
#' @description Simulate waiting times for a non-homogeneous Poisson process.
#' @param rate_fun A function specifying the rate as a function of time.
#' @param n Number of samples.
#' @return Numeric vector of event times.
#' @examples
#' # Example usage (pseudo-code):
#' # times <- generateNonHomogeneousExp(rate_fun, n)
#' @export
generateNonHomogeneousExp <- function(rate_fun, n) {
  stop("Not yet implemented.")
}

#' Non-homogeneous Exponential Random Sampling
#'
#' @description Sample event times from a non-homogeneous exponential process.
#' @param rate_fun A function specifying the rate as a function of time.
#' @param n Number of samples.
#' @return Numeric vector of event times.
#' @examples
#' # Example usage (pseudo-code):
#' # samples <- nhExpRand(rate_fun, n)
#' @export
nhExpRand <- function(rate_fun, n) {
  stop("Not yet implemented.")
}

#' Generate Non-Homogeneous Exponential Events
#'
#' @description Simulate waiting times for a non-homogeneous Poisson process.
#' @param rate_fun A function specifying the rate as a function of time.
#' @param n Number of samples.
#' @return Numeric vector of event times.
#' @examples
#' # Example usage (pseudo-code):
#' # times <- generateNonHomogeneousExp(rate_fun, n)
#' @export
generateNonHomogeneousExp <- function(rate_fun, n) {
  stop("Not yet implemented.")
}

#' Non-homogeneous Exponential Random Sampling
#'
#' @description Sample event times from a non-homogeneous exponential process.
#' @param rate_fun A function specifying the rate as a function of time.
#' @param n Number of samples.
#' @return Numeric vector of event times.
#' @examples
#' # Example usage (pseudo-code):
#' # samples <- nhExpRand(rate_fun, n)
#' @export
nhExpRand <- function(rate_fun, n) {
  stop("Not yet implemented.")
}
