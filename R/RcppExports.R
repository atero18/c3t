# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @description Recherche de la position d'un élément d'une matrice symétrique
#' stockée sous forme d'un vecteur.
#'
#' @param i,j position en ligne et colonne (démarre à 0L) de l'élément.
#' (entiers compris entre 0 et `dim - 1`)
#' @param dim dimension de la matrice. (entier positif)
#' @param aDefautDiag `true` si la matrice a une valeur par défaut
#' sur sa diagonale (ex : 0.0 pour une matrice de distance), `false` sinon.
#' (booléen)
#' @noRd
NULL

#' @description Effectue le calcul d'une distance inter-clusters avec
#' une matrice de distances `SymVMat` (R).
#' @name compSMSymVMat
#' @noRd
NULL

#' @description Effectue le calcul d'une distance inter-clusters avec
#' une matrice de distances `SymMMat` (R).
#' @name compSMSymMMat
#' @noRd
NULL

#' @description Calcule toutes les distances inter-clusters demandées
#' à partir de distances inter-éléments stockées dans un objet
#' `SymVMat` (R).
#' @name distanceInterSymVMat
#' @noRd
distanceInterSymVMat <- function(values, partition, indexs, dim, aDefautDiag, defaultDiag, comp) {
    .Call(`_c3t_distanceInterSymVMat`, values, partition, indexs, dim, aDefautDiag, defaultDiag, comp)
}

#' @description Calcule toutes les distances inter-clusters demandées
#' à partir de distances inter-éléments stockées dans un objet
#' `SymVMat` (R).
#' @name distanceInterSymMMat
#' @noRd
distanceInterSymMMat <- function(values, partition, indexs, comp) {
    .Call(`_c3t_distanceInterSymMMat`, values, partition, indexs, comp)
}

#' @title Find the nearest neighbor of some points
#' @description For a set of points (`subsetPoints`), find the nearest
#' neighbor in another set (`subsetNeighbors`). Distances are given in a
#' `p' x n'` matrix (`distances`), with `p' >= p` (p the number of points)
#' and `n' >= n` (n the number of neighbors). If points and neighbors come from
#' the same set, put inner to `TRUE` for assuring the nearest neighbor of a
#' point will not be itself.
#' @param distances distances between points and neighbors. A row is for a
#' point and a column a neighbor.
#' @param subsetPoints A vector of strictly positive integers giving the
#' position of the points the nearest neighbor must be found. Positions
#' starts from 1 like in R. Must be at least of length one and it would be
#' better if values are unique.
#' @param subsetNeighbors A vector of strictly positive integers giving the
#' position of the elements of the set that can be consider has a neighbor
#' for each point of `subsetPoints`. Positions starts from 1 like in R.
#' Must be at least of length one and it would be
#' better if values are unique.
#' @param contiguity In the case there are contiguity constraints
#' between points. If an element of `subsetPoints` and one
#' of `subsetNeighbors` are not contiguous then the second element cannot
#' be a neighbor of the first one. If no contiguity constraints must be
#' considered a matrix with 0 rows and 0 columns must be sent. Otherwise
#' must be a logical matrix with same number of rows and columns than
#' `distante`.
#' @returns a vector of length `p` with for each point the position of its
#' nearest neighbor (0 if it doesn't exist). Positions starts from 1 like in R.
#' @keywords internal
#' @name cpp_nearest_neighbor
NULL

#' @describeIn cpp_nearest_neighbor case when `distances` is a matrix
#' @keywords internal
nearest_neighbor_matrix <- function(distances, subsetPoints, subsetNeighbors, contiguity, inner) {
    .Call(`_c3t_nearest_neighbor_matrix`, distances, subsetPoints, subsetNeighbors, contiguity, inner)
}

