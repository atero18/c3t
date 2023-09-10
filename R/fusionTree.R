fusionNode <-
  setRefClass("FusionNode",
              fields = list(pb = "pbCon",
                            childrenNodes = "list",
                            parentNode = "ANY",
                            regionalisation = "vector",
                            clustersSizes = "vector",
                            checkMinSizeConst = "logical",
                            contiguities = "matrix",
                            nbClusters = "numeric"))


setAs(
  "FusionNode", "Partition",
  function(from)
  {
    constructor_Partition(pb = from$pb, partition = from$regionalisation,
                          method = "fusionTree", contiguity = TRUE,
                          minConstraint = from$checkMinSizeConst,
                          maxConstraint = TRUE)
  }
)

setAs(
  "FusionNode", "AHCTree",
  function(from)
  {
    partitions <- list()
    nbClusters <- from$nbClusters

    partitions[[as.character(nbClusters)]] <-
      as(from$regionalisation, "Partition")

    childNode <- from
    while (!is.null(childNode$parentNode))
    {
      childNode <- childNode$parentNode
      nbClusters  <- nbClusters + 1L

      partitions[[as.character(nbClusters)]] <-
        as(childNode$regionalisation, "Partition")
    }

    AHCTree$new(pb = from$pb, partitions = partitions)
  }
)


fusionTree <-
  setRefClass("FusionTree",
              fields = list(pb = "pbCon",
                            rootNode = "FusionNode",
                            feasibleNodes = "list",
                            unfeasibleLeaves = "list",
                            nbCrossedBranches = "numeric",
                            search = "list"))


#' @importFrom tibble as_tibble
setMethod(
  "as_tibble",
  signature(x = "FusionTree"),
  function(x)
  {
    if (length(x$feasibleNodes) > 0L)
      partitions <- lapply(x$feasibleNodes, as, "Partition")

    else
      partitions <- list()

    Partition_list_to_tibble(partitions)
  }
)


CHILDCHOICES <- c("random", "first_valid_node", "distance")

create_child_node <- function(node,
                              k,
                              childChoice = "first_valid_node",
                              tree = NULL)
{

  childNode <- node$copy(shallow = TRUE)
  childNode$parentNode <- node

  # Define the pair of clusters used in the fusion
  cluster1 <- node$contiguities[k, 1L]
  nomCluster1 <- as.character(cluster1)
  cluster2 <- node$contiguities[k, 2L]
  nomCluster2 <- as.character(cluster2)

  # Update regionalisation
  regionalisation <- childNode$regionalisation
  regionalisation[regionalisation == cluster2] <- cluster1
  childNode$regionalisation <- regionalisation

  # Update contiguities
  contiguities <- update_contiguity_list(node$contiguities,
                                         cluster1,
                                         cluster2,
                                         childNode$clustersSizes)

  # Update clusters sizes
  childNode$clustersSizes[nomCluster1] <-
  childNode$clustersSizes[nomCluster1] +
  childNode$clustersSizes[nomCluster2]

  childNode$clustersSizes <-
    childNode$clustersSizes[names(childNode$clustersSizes) != nomCluster2]

  childNode$nbClusters <- node$nbClusters - 1L

  # Remove of the fusions that would not check the maximum size constraint
  contiguities <- size_respect_contiguities(node$pb,
                                            contiguities,
                                            childNode$clustersSizes)

  reduction <- priorisation <- FALSE
  orderByDistance <- childChoice == "distance"
  alea <- childChoice == "random"
  # Calculation of missing linkage distances (if necessary)
  if (orderByDistance)
  {
    masqueDistancesManquantes <- is.na(contiguities[, "distance"])
    if (any(masqueDistancesManquantes))
    {
      contiguities[masqueDistancesManquantes, "distance"] <-
        .calcul_distances_inter(childNode$pb, regionalisation,
                                contiguities[masqueDistancesManquantes,
                                             1L:2L, drop = FALSE],
                                tree$search$linkage)
    }
  }

  # Sorting fusions depending on `childChoice`
  if (alea || reduction || priorisation || orderByDistance)
  {
    contiguities <-
      transform_contiguity_list(node$pb, contiguities,
                                childNode$clustersSizes,
                                table(regionalisation),
                                priorisation = priorisation,
                                alea = alea,
                                reduction = reduction,
                                orderByDistance = orderByDistance)
  }

  childNode$contiguities <- contiguities

  # We check if the child node verifies the minimum size constraint
  # (which is automatically check if its parent verified it)
  childNode$checkMinSizeConst <-
    node$checkMinSizeConst || all(childNode$clustersSizes >= node$pb$m)

  if (childNode$checkMinSizeConst && !is.null(tree))
    tree$feasibleNodes <- c(tree$feasibleNodes, childNode)

  childNode
}

SEARCHS <- c("random", "DFS", "BFS")

#' @importFrom fastmap faststack fastqueue
.create_node_list <- function(tree, search)
{
  if (search == "DFS")
  {
    nodesToEvaluate <- faststack()
    nodesToEvaluate$push(tree$rootNode)
  }

  else if (search == "BFS")
  {
    nodesToEvaluate <- fastqueue()
    nodesToEvaluate$add(tree$rootNode)
  }

  else if (search == "random")
    nodesToEvaluate <- list(tree$rootNode)
  else
  {
    stop("Unrecognized `search` mode")
  }

  nodesToEvaluate
}

.update_node_list <- function(nodesToEvaluate, childrenNodes,
                              search)
{
  if (search == "DFS")
    nodesToEvaluate$mpush(.list = childrenNodes)

  else if (search == "BFS")
    nodesToEvaluate$madd(.list = childrenNodes)

  else
    nodesToEvaluate <- c(nodesToEvaluate, childrenNodes)

  nodesToEvaluate
}

search_fusionTree <- function(tree,
                              search = "random",
                              childChoice = "first_valid_node")
{

  nodesToEvaluate <- .create_node_list(tree, search)

  it <- 0L
  while (length(nodesToEvaluate) > 0L && !is_stopCriterion_checked(tree))
  {
    it <- it + 1L
    if (search == "DFS")
      node <- nodesToEvaluate$pop()

    else if (search == "BFS")
      node <- nodesToEvaluate$remove()

    else if (search == "random")
    {
      posNextNode <- sample(seq_along(nodesToEvaluate), size = 1L)
      node <- nodesToEvaluate[[posNextNode]]
      nodesToEvaluate <- nodesToEvaluate[-posNextNode]
    }

    node <- create_child_node(node$parentNode,
                              node$idFusion,
                              childChoice,
                              tree)

    nbChildrenNodes <- nrow(node$contiguities)
    if (node$checkMinSizeConst)
      tree$feasibleNodes[[length(tree$feasibleNodes) + 1L]] <- node

    if (nbChildrenNodes == 0L)
    {
      tree$nbCrossedBranches <- tree$nbCrossedBranches + 1L
      if (!node$checkMinSizeConst)
        tree$unfeasibleLeaves <- c(tree$unfeasibleLeaves, node)

      next
    }

    order <- switch(search,
                    DFS = seq(nbChildrenNodes, 1L),
                    seq_len(nbChildrenNodes))

    childrenNodes <- lapply(order, function(k)
      list(parentNode = node,
           idFusion = k))

    nodesToEvaluate <- .update_node_list(nodesToEvaluate,
                                         childrenNodes,
                                         search)
  }

  return(tree)
}


recursive_DFS_fusionTree <- function(tree,
                                     node,
                                     childChoice = "first_valid_node")
{
  if (node$checkMinSizeConst)
    tree$feasibleNodes[[length(tree$feasibleNodes) + 1L]] <- node

  nbChildrenNodes <- nrow(node$contiguities)

  if (nbChildrenNodes == 0L)
  {
    tree$nbCrossedBranches <- tree$nbCrossedBranches + 1L
    tree$unfeasibleLeaves[[length(tree$unfeasibleLeaves) + 1L]] <- node
    return(FALSE)
  }

  if (is_stopCriterion_checked(tree))
    return(TRUE)

  stop <- FALSE
  k <- 1L

  while (!stop && k <= nbChildrenNodes)
  {
    childNode <- create_child_node(node, k, childChoice, tree)
    node$childrenNodes[[k]] <- childNode

    recursive_DFS_fusionTree(tree, node, childChoice)

   if (is_stopCriterion_checked(tree))
     stop <- TRUE

    k <- k + 1L
  }

  FALSE
}


STOPCRITERIA <- c("full",
                  "first_valid_node",
                  "first_valid_branch")

is_stopCriterion_checked <- function(tree)
{
  criterion <- tree$search$stopCriterion
  if (criterion == "full")
    return(FALSE)

  else if (criterion == "first_valid_node" && length(tree$feasibleNodes) > 0L)
    return(TRUE)

  else if (criterion == "first_valid_branch" &&
          length(tree$feasibleNodes) > 0L &&
          tree$nbCrossedBranches > 0L)
    return(TRUE)

  FALSE
}

#' @importFrom checkmate assertChoice assertFlag
search_fusion_tree_pb <- function(pb, search = "DFS",
                                  stopCriterion = "first_valid_branch",
                                  childChoice = "distance",
                                  regionalisation = seq_len(pb$n()),
                                  linkage = "complete",
                                  useRecursivity = FALSE)
{
  # Checking arguments
  validObject(pb)
  assertChoice(search, SEARCHS)
  assertChoice(stopCriterion, STOPCRITERIA)
  assertChoice(childChoice, CHILDCHOICES)

  if (search == "DFS")
    assertFlag(useRecursivity)

  contiguities <- clusters_contiguity_list(regionalisation, pb)

  if (childChoice == "distance")
  {
    assertLinkage(linkage)
    contiguities <- cbind(contiguities, NA)
    colnames(contiguities)[3L] <- "distance"
  }

  if (pb$hasMaxConstraint())
  {
    contiguities <- cbind(contiguities, NA)
    colnames(contiguities)[3L] <- "fusionSize"

    clustersSizes  <- .clusters_sizes(regionalisation, pb$sizes)

    if (any(clustersSizes > pb$M))
    {
      stop("Maximum size constraint is not verified")
    }

    contiguities <- size_after_fusion(contiguities,
                                      clustersSizes)
    contiguities <- size_respect_contiguities(pb,
                                              contiguities,
                                              clustersSizes)
  }

  if (nrow(contiguities) == 0L)
  {
    stop("No possible fusion")
  }

  if (childChoice == "distance")
  {
    contiguities[, "distance"] <-
      .calcul_distances_inter(pb,
                              regionalisation,
                              contiguities,
                              linkage)
    contiguities <-
      contiguities[order(contiguities[, "distance"]), ]
  }

  else if (childChoice == "random")
  {
    contiguities <-
      contiguities[sample(seq_len(nrow(contiguities))), ]
  }


  rootNode <- fusionNode$new(pb = pb, childrenNodes = list(),
                             regionalisation = regionalisation,
                             clustersSizes = clustersSizes,
                             contiguities = contiguities,
                             checkMinSizeConst =
                               !pb$hasMinConstraint() ||
                               all(clustersSizes >= pb$m),
                             nbClusters = length(clustersSizes),
                             parentNode = NULL)

  if (rootNode$checkMinSizeConst)
    feasibleNodes <- list(rootNode)
  else
    feasibleNodes <- list()

  tree <- fusionTree$new(pb = pb, rootNode = rootNode,
                         feasibleNodes = feasibleNodes,
                         search =
                           list(linkage = linkage,
                                stopCriterion = stopCriterion),
                         nbCrossedBranches = 0L)

  if (search == "DFS" && useRecursivity)
  {
    recursive_DFS_fusionTree(tree = tree,
                             node = rootNode,
                             childChoice = childChoice)
  }

  else
  {
    search_fusionTree(tree = tree,
                      search = search,
                      childChoice  = childChoice)

  }

  return(tree)
}


#' Search of a fusion tree under connectivity and/or size constraints.
#'
#' @inheritParams arguments_problem
#' @param regionalisation regionalisation to start the tree traversal.
#' @param search Type of traversal to apply on the fusion tree.
#' Three possibilities:
#' * "`BFS`" (Depth First Search)
#' * "`DFS`" (Breadth First Search)
#' * "`random`" (random)
#' @param childChoice In the case of breadth-first or depth-first search,
#' the way in which child nodes are ordered (and accessed). Multiple choices:
#' * "`random`": Children are ordered randomly.
#' * "`first_valid_node`": Initial order is preserved.
#' * "`distance`": Nodes are ordered by increasing inter-cluster distance.
#' @param stopCriterion Indicates when the traversal should stop.
#' The possibilities are as follows:
#' * "`full`": The algorithm stops when all possible fusions from the initial
#' partition have been evaluated.
#' * "`first_valid_node`": Stop when a first feasible solution is found.
#' Useful when a minimum size constraint is present.
#' * "`first_valid_branch`": Stop as soon as a node without fusion has been
#' selected and a feasible solution has been found.
#' @details The triplet
#' (`search` = "`DFS`", `childChoice` = "`distance`",
#' `stopCriterion` = "`first_valid_branch`")
#' corresponds to an Agglomerative Hierarchical Clustering (AHC) with
#' contiguity and size constraints, not necessarily stopping at the fusion of
#' two clusters but when a feasible partition is found.
#' The two methods are equivalent when there is no minimum size constraint.
#' `r badge('experimental')`
#' @name search_fusion_tree
#' @rdname search_fusion_tree
#' @export
search_fusion_tree <- function(distances = NULL, contiguity = NULL,
                               sizes = NULL,
                               d = NULL, data = NULL, m = 0L, M = Inf,
                               standardQuant = FALSE, binarQual = FALSE,
                               search = "DFS",
                               stopCriterion = "first_valid_branch",
                               childChoice = "distance",
                               regionalisation = NULL,
                               linkage = "complete")
{
  pb <- constructor_pbCon(distances = distances,
                          contiguity = contiguity,
                          sizes = sizes,
                          d = d, data = data,
                          m = m, M = M,
                          standardQuant = standardQuant,
                          binarQual = binarQual)

  if (is.null(regionalisation))
    regionalisation <- seq_len(pb$n())

  tree <- search_fusion_tree_pb(pb, search = search,
                                stopCriterion = stopCriterion,
                                childChoice = childChoice,
                                regionalisation = regionalisation,
                                linkage = linkage)

  as_tibble(tree)

}

#' @describeIn search_fusion_tree Shortcut for `search_fusion_tree`
#' in the case of an AHC / AHR search.
#' @param linkage used when `chilChoice = ` `"distance"`. Can
#' be a function or a string. This determines the linkage criterion and
#' so the selection of the children nodes.
#' If using custom  distance function, it must take the pairwise distance
#' matrix as an argument and return the linkage distance. Default implemented
#' linkages can be seen with [available_linkages()].
#' @seealso [AHR()]
#' @export
AHR_fusion_tree <- function(distances = NULL, contiguity = NULL,
                            sizes = NULL, d = NULL, data = NULL,
                            m = 0L, M = Inf,
                            standardQuant = FALSE, binarQual = FALSE,
                            regionalisation = NULL,
                            linkage = "complete")
{

  search <- "DFS" # nolint: object_usage_linter
  stopCriterion <- "first_valid_branch" # nolint: object_usage_linter
  childChoice <- "distance" # nolint: object_usage_linter

  eval(body(search_fusion_tree), envir = environment())
}
