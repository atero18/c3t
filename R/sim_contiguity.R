# Generates a contiguity matrix for an arbitrary grid. Different contiguity
# modes are accepted. The returned matrix is a sparse (Matrix package),
# symmetric, boolean matrix.
# All elements on the diagonal are set to True. The default value of the
# sparse matrix is False.

#' @param contiguityType type of contiguity to apply for the grid.
#' Two options:
#' * "Rook": two cells are contiguous if they share a side.
#' * "Queen": two cells are contiguous if they share a side or a vertex.
#' @param x_int,y_int dimensions of the grid (positive integers)
#' @name grid_contiguity_matrix
NULL


#' @importFrom checkmate assertCount assertString assertChoice
#' @importFrom Matrix sparseMatrix
grid_contiguity_matrix <- function(x_int, y_int, # nolint: cyclocomp_linter
                                   contiguityType = "Queen")
{

  # Checking arguments
  assertCount(x_int, positive = TRUE)
  assertCount(y_int, positive = TRUE)
  assertString(contiguityType)
  contiguityType <- tolower(contiguityType)
  assertChoice(contiguityType, c("queen", "rook"))

  noms_vec <- apply(expand.grid(seq_len(x_int), seq_len(y_int)), 1L,
                    function(r) paste0("(", r[1L], ",", r[2L], ")"))

  # Création d'une matrice carrée creuse symétrique booléenne (package Matix)
  M_mat <- sparseMatrix(1L, 1L,
                        x = FALSE,
                        dims = c(x_int * y_int, x_int * y_int),
                        dimnames = list(noms_vec, noms_vec))

  # Pour chaque carré de la grille, recherche de ses voisins.
  for (i_int in seq_len(x_int))
  {
    for (j_int in seq_len(y_int))
    {
      pos_elem <- xy_to_rank_grid(i_int, j_int, x_int, y_int)

      # Gestion des voisins type "Rook" (haut bas gauche droite)
      if (i_int > 1L)
      {
        pos_elem_gauche <- xy_to_rank_grid(i_int - 1L, j_int,
                                           x_int, y_int)

        M_mat[pos_elem, pos_elem_gauche] <- TRUE
      }
      if (i_int < x_int)
      {
        pos_elem_droit <- xy_to_rank_grid(i_int + 1L, j_int,
                                          x_int, y_int)

        M_mat[pos_elem, pos_elem_droit] <- TRUE
      }
      if (j_int > 1L)
      {
        pos_elem_haut <- xy_to_rank_grid(i_int, j_int - 1L,
                                         x_int, y_int)

        M_mat[pos_elem, pos_elem_haut] <- TRUE
      }
      if (j_int < y_int)
      {
        pos_elem_bas <- xy_to_rank_grid(i_int, j_int + 1L,
                                        x_int, y_int)

        M_mat[pos_elem, pos_elem_bas] <- TRUE
      }

      # Gestion des voisins supplémentaires "Queen" (diagonale)
      if (contiguityType == "queen")
      {
        if (i_int > 1L && j_int > 1L)
        {
          pos_elem_haut_gauche <-
            xy_to_rank_grid(i_int - 1L, j_int - 1L,
                            x_int, y_int)

          M_mat[pos_elem, pos_elem_haut_gauche] <- TRUE
        }
        if (i_int > 1L && j_int < y_int)
        {
          pos_elem_bas_gauche <- xy_to_rank_grid(i_int - 1L, j_int + 1L,
                                                 x_int, y_int)

          M_mat[pos_elem, pos_elem_bas_gauche] <- TRUE
        }
        if (i_int < x_int && j_int > 1L)
        {
          pos_elem_haut_droit <- xy_to_rank_grid(i_int + 1L, j_int - 1L,
                                                 x_int, y_int)

          M_mat[pos_elem, pos_elem_haut_droit] <- TRUE
        }
        if (i_int < x_int && j_int < y_int)
        {
          pos_elem_bas_droit <- xy_to_rank_grid(i_int + 1L, j_int + 1L,
                                                x_int, y_int)

          M_mat[pos_elem, pos_elem_bas_droit] <- TRUE
        }
      }

    }
  }

  # Tout élément est "voisin" de lui-même
  diag(M_mat) <- TRUE

  return(M_mat)

}
