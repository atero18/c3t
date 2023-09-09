#' @include abstractSymMat.R
#' @include distMat.R


#' @keywords internal
#' @importFrom igraph gorder
setRefClass(
  "pbCon",
  fields = list(sizes = "numeric", N = "numeric",
                # should be `igraph` but have some issues
                # using [methods::setOldClass()] # nolint
                contiguity = "ANY",
                contiguityMatrix = "matriceOrNULL",
                m = "numeric", M = "numeric",
                dataContiguity = "list"),
  contains = "DistMat",
  methods = list(
    initialize = function(distances = NULL, d = NULL, data = NULL,
                          contiguity = NULL, contiguityMatrix = NULL,
                          sizes = seq_len(nrow(distances)),
                          m = 0.0, M = Inf, dataContiguity = list(), ...)
    {
      if (is.null(distances))
        return()

      distances <<- distances # nolint: undesirable_operator_linter
      d <<- d # nolint: undesirable_operator_linter
      data <<- data # nolint: undesirable_operator_linter

      setContiguity(contiguity, contiguityMatrix, dataContiguity)

      setSizes(sizes)

      setSizeConstraints(m, M)
    },
    setSizes = function(sizes)
    {
      assertSizes(sizes)
      sizes <<- sizes # nolint: undesirable_operator_linter
      N <<- sum(sizes) # nolint: undesirable_operator_linter

    },
    #' @importFrom checkmate assertNumber
    setSizeConstraints = function(tempm = 0.0, tempM = Inf)
    {
      if (is.null(tempm) || (length(tempm) == 1L && is.na(tempm)))
         tempm <- 0.0
      else
      {
        if (!m_valid(tempm))
          stop("Incorrect value for m")

        tempm <- simplify_minSizeConst(tempm, sizes)
      }

      m <<- tempm # nolint: undesirable_operator_linter


      if (is.null(tempM) || (length(tempM) == 1L && is.na(tempM)))
        tempM <- Inf
      else
      {
        if (!M_valid(tempM))
          stop("Incrorrect value for M")

        tempM <- simplify_maxSizeConst(tempM, sizes)
      }

      M <<- tempM # nolint: undesirable_operator_linter
    },

    setContiguity = function(newContGraph,
                             newContMat = NULL,
                             nouvDataContiguite = list())
    {
      if (is.null(newContGraph))
      {
        newContGraph <- complete_contiguity_graph(n())
        newContMat <- complete_contiguity_matrix(n())
      }

      # nolint start: undesirable_operator_linter
      contiguity <<- newContGraph

      contiguityMatrix <<- newContMat

      dataContiguity <<- nouvDataContiguite
      # nolint end
    },
    setCompleteContiguity = function()
    {
      setContiguity(complete_contiguity_graph(n()), NULL,
                    list(completementConnecte = TRUE))
    },
    setOldContiguity = function(dataContiguity_list)
    {
      setContiguity(dataContiguity_list$contiguity,
                    dataContiguity_list$contiguityMatrix,
                    dataContiguity_list$dataContiguity)
    }
  )
) -> pbCon # nolint : assignment_linter

pbCon$methods(
  #' @importFrom igraph is_igraph
  isValid = function()
  {
    if (!m_valid())
      return("invalid minimum size constraint")

    if (!M_valid())
      return("invalid maximum size constraint")

    if (!is_igraph(contiguity))
      return("`contiguity` must be an `igraph` object")


    if (n() != gorder(contiguity))
      return(paste0("Number of elements differs between distance matrix and ",
                    "contiguity graph"))

    if (!valid_sizes())
      return("Incorrect size vector")

    TRUE
  }
)

verif_pbCon <- function(object)
{
  if (length(nrow(object)) == 0L)
    return(TRUE)

  object$isValid()
}

setValidity("pbCon", verif_pbCon)

#' @importFrom checkmate assertFlag
is_pbCon <- function(x, checkValidity = FALSE)
{
  assertFlag(checkValidity)
  if (!inherits(x, "pbCon"))
    return(FALSE)

  if (checkValidity)
    return(validObject(x))

  TRUE
}

pbCon$methods(
  hasMinConstraint = function() m > 0.0,
  hasMaxConstraint = function() is.finite(M),
  hasSizeConstraint = function() hasMinConstraint()  || hasMaxConstraint()
)

# Gestion des tailles
# -- Vérification
pbCon$methods(
  valid_sizes = function(sizes = .self$sizes) testSizes(sizes, len = n())
)

# Gestion des tailles
# -- Vérification
pbCon$methods(
  m_valid = function(m = .self$m)
  {
    testMinSizeConstraints(m, sizes, connected_components(), len = 1L)
  },
  M_valid = function(M = .self$M)
  {
    testMaxSizeConstraints(M, sizes, len = 1L)
  }
)

# -- Simplification

pbCon$methods(
  simplify_min_constraint = function(m_vec)
  {
    simplify_minSizeConst(m_vec, sizes)
  },
  simplify_max_constraint = function(M_vec)
  {
    simplify_maxSizeConst(M_vec, sizes, connected_components())
  }
)

# -- Nombre de clusters ne vérifiant pas les contraintes
pbCon$methods(
  tooSmallClusters = function(partition)
  {
    if (!hasMinConstraint())
      return(integer(0L))

    clusters_trop_petits(m, partition, sizes)
  },
  nbTooSmallClusters = function(partition)
  {
    if (!hasMinConstraint())
      return(0L)

    nb_clusters_trop_petits(m, partition, sizes)
  },
  tooBigClusters = function(partition)
  {
    if (!hasMaxConstraint())
      return(integer(0L))

    clusters_trop_gros(M, partition, sizes)
  },
  nbTooBigClusters = function(partition)
  {
    if (!hasMaxConstraint())
      return(0L)

    nb_clusters_trop_gros(M, partition, sizes)
  }
)

# -- Récupération des données
pbCon$methods(
  getContiguity = function()
  {
    list(contiguity        = .self$contiguity,
         contiguityMatrix = .self$contiguityMatrix,
         dataContiguity    = .self$dataContiguity)
  }
)

# -- On regarde si la contrainte est active


pbCon$methods(
  #' @importFrom igraph degree
  isFullyConnected = function(calcul = TRUE)
  {
    if (isContiguityPropertyCalculated("fullyConnected"))
      return(getContiguityProperty("fullyConnected"))

    if (calcul)
    {
      degrees <- degree(contiguity,
                        mode = "all",
                        loops = FALSE,
                        normalized = FALSE)
      fullyConnected <- all(degrees == n() - 1L)
      setContiguityProperty("fullyConnected", fullyConnected)
      return(fullyConnected)
    }

    else
      return(FALSE)

  },
  aContrainteContiguite = function()
  {
    !isFullyConnected(calcul = TRUE)
  }
)


# -- Réception et ajout des propriétés de contiguïté
pbCon$methods(
  isContiguityPropertyCalculated = function(propriete)
  {
    propriete %in% names(dataContiguity)
  },

  getContiguityProperty = function(propriete)
  {
    dataContiguity[[propriete]]
  },

  setContiguityProperty = function(propriete, value)
  {
    .self$dataContiguity[[propriete]] <- value
  }
)

# - Recherche du nombre de composantes connexes
pbCon$methods(
  connected_components = function()
  {
    if (isContiguityPropertyCalculated("connectedComponents"))
      return(getContiguityProperty("connectedComponents"))

    if (isFullyConnected(calcul = FALSE) || is.null(contiguity))
      connectedComp <- rep(1L, n())

    else
      connectedComp <- components(contiguity, "weak")$membership

    setContiguityProperty("connectedComponents", connectedComp)

    return(connectedComp)
  }
)


#' @importFrom igraph components
setGeneric("components", igraph::components)

#' @importFrom igraph components
setMethod("components", signature(graph = "pbCon", mode = "ANY"),
          function(graph, mode) graph$connected_components()
)


pbCon$methods(
  nbComposantesConnexes = function()
  {
    if (isFullyConnected(calcul = FALSE))
        return(1L)

      else
        return(length(unique(connected_components())))
  }
)

#' @importFrom igraph count_components
setGeneric("count_components", igraph::count_components)

#' @importFrom igraph count_components
setMethod(
  "count_components",
  signature(graph = "pbCon", mode = "ANY"),
  function(graph, mode) graph$nbComposantesConnexes()
)

pbCon$methods(isConnected = function(calcul = TRUE)
{
  if (isFullyConnected(calcul = FALSE))
    return(TRUE)

  if ((calcul || isContiguityPropertyCalculated("connectedComponents"))
      && nbComposantesConnexes() == 1L)
    return(TRUE)

  return(FALSE)
})

# -- Recherche du nombre de voisins

# Vérification des régionalisations
pbCon$methods(
  isPartition = function(partition) testPartition(partition, n()),
  isRegionalisation = function(partition)
  {
    if (isFullyConnected(calcul = FALSE))
      return(isPartition(partition))

    else
      return(is_regionalisation(partition, contiguity))
  },
  isRegionalisationInitValide = function(partition)
  {
    isRegionalisation(partition) &&
      nbTooBigClusters(partition) == 0L
  },
  isFeasibleSolution = function(partition)
  {
    isRegionalisation(partition) &&
      nbTooSmallClusters(partition) == 0L &&
      nbTooBigClusters(partition)   == 0L
  }
)

# Récupération de la liste de fusions inter-classes
pbCon$methods(
  listeContiguiteClasses = function(partition,
                                    withFusionSize = TRUE)
  {
    contiguiteClasses <-
      as.matrix(clusters_contiguity_list(partition, contiguity))

    if (hasMaxConstraint() || withFusionSize)
    {
      tailleFusion <-
        size_after_fusion(contiguiteClasses, clusters_sizes(partition,
                                                            sizes))
    }

    if (hasMaxConstraint())
    {
      masqueFusionOk <- tailleFusion <= M
      contiguiteClasses <- contiguiteClasses[masqueFusionOk, , drop = FALSE]
      if (withFusionSize)
        tailleFusion <- tailleFusion[masqueFusionOk]
    }

    if (withFusionSize)
    {
      noms <- colnames(contiguiteClasses)
      contiguiteClasses <- cbind(contiguiteClasses,
                                 tailleFusion)
      if (!is.null(noms))
        colnames(contiguiteClasses) <- c(noms, "fusionSize")

    }

    contiguiteClasses
  }
)

#' @inheritParams arguments_distMat
#' @param contiguity A contiguity matrix or an `igraph` contiguity graph. If not
#' provided, the problem is considered completely contiguous (all elements are
#' neighbors of each other).
#' @param sizes Represents the size of each element. By default, it is set
#' to `1` for each element (the size of a cluster becomes its cardinal).
#' All data must be positive or zero. (positive real numeric vector)
#' @param m Minimum size constraint. Must be positive or zero and small enough
#' for the problem to be feasible. Default is `0` (no constraint).
#' (positive number)
#' @param M Maximum size constraint. Must be positive, superior or equal to `m`
#' and large enough for the problem to be feasible.
#' Default is `Inf` (no constraint). (positive number)
#' @name arguments_problem
#' @keywords internal
#' @returns An object of class `pbCon`.
constructor_pbCon <- function(distances = NULL, contiguity = NULL,
                              sizes = NULL,
                              d = NULL, data = NULL, m = 0.0, M = Inf,
                              standardQuant = FALSE, binarQual = FALSE,
                              storageMode = "matrix")
{
  mat <- constructor_DistMat(distances = distances,
                             d = d,
                             data = data,
                             standardQuant = standardQuant,
                             binarQual = binarQual,
                             storageMode = storageMode)

  n <- nrow(mat)
  contiguityMatrix <- NULL
  dataContiguity <- list()
  if (is.matrix(contiguity) || inherits(contiguity, "Matrix"))
    contiguity <- contiguity_matrix_to_graph(contiguity)

  else if (is.null(contiguity))
  {
    contiguity <- complete_contiguity_graph(n)
    dataContiguity[["fullyConnected"]] <- TRUE
  }
  else
    assertContiguityGraph(contiguity)

  if (is.null(sizes))
    sizes <- rep(1.0, n)

  pbCon$new(distances = mat$distances, d = mat$d, data = mat$data,
            contiguity = contiguity,
            contiguityMatrix = contiguityMatrix,
            sizes = sizes, m = m, M = M,
            dataContiguity = dataContiguity)
}

#' Get the contiguity matrix of a `pbCon` object.
#'
#' @param pb A `pbCon` object representing a constraint problem.
#' @returns The contiguity matrix of the `pbCon` object.
#' @noRd
getmatContPbCon <- function(pb)
{
  if (is.null(pb$contiguityMatrix))
  {
    if (pb$isFullyConnected(calcul = FALSE))
      return(complete_contiguity_matrix(pb))
    else
      graph_to_contiguity_matrix(pb)
  }

  return(pb$contiguityMatrix)
}

pbCon$methods(
  #' Divide a `pbCon` object into subproblems for connected components.
  #'
  #' This function divides the `pbCon` object into multiple subproblems based on
  #' its connected components. Each subproblem will correspond to one connected
  #' component in the `pbCon` object's contiguity graph.
  #' @returns A list of `pbCon` objects, each corresponding to a connected
  #' component in the `pbCon` object's contiguity graph.
  #' @importFrom igraph induced_subgraph
  #' @noRd
  split_connected_spb = function()
  {
    # Retrieve connected components
    connectedComponents <- connected_components()

    # Create a problem for each connected component
    tapply(seq_len(n()), connectedComponents,
           function(listeElements) {
      grapheCont <- induced_subgraph(contiguity, listeElements)
      distM <- sousMatrice_distMat(.self, listeElements)
      pbCon$new(distances = distM$distances, d = distM$d, data = distM$data,
                contiguity = grapheCont, sizes = sizes[listeElements],
                m = m, M = M,
                dataContiguity =
                  list(connectedComp = rep(1L, length(listeElements))))
           })
  }
)

#' @rdname abstractSymMat_replace
#' @keywords internal
setReplaceMethod(
  "[",
  signature(x = "pbCon", i = "matrix", j = "missing", value = "numeric"),
  function(x, i, j, value) callNextMethod()
)

pbCon$methods(
  #' @importFrom checkmate assertFlag
  quality_partition = function(partition, criterion,
                               linkage = NULL,
                               connected = FALSE,
                               valueOnly = TRUE)
  {
    assertFlag(connected)
    connected <- connected && !isConnected()

    clustering_criterion(partition, criterion,
                         .self, d, data,
                         FALSE, FALSE,
                         linkage, contiguity, connected,
                         valueOnly = valueOnly)
  },
  update_quality_partition = function(partitionBefore, criterion,
                                      dataCriterion,
                                      donor, receiver, givenElements,
                                      distances = distances,
                                      d = d,
                                      dataElements = data,
                                      linkage = NULL)
  {
    update_criterion(partitionBefore, criterion,
                     dataCriterion,
                     donor, receiver, givenElements,
                     .self, d, data,
                     linkage)
  }
)
