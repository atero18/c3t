#' @include partition_class.R
#' @include fusion_constraint.R
#' @include AHCTree.R
#' @include pbCon.R


#' @importFrom igraph V E make_graph get.edge.ids neighbors ends gsize
#' @importFrom igraph V<- E<- set_edge_attr vertex_attr<-
#' @importFrom igraph delete_edges add_edges get.edge.attribute
#' @importFrom igraph get.vertex.attribute set_vertex_attr
#' @importFrom checkmate testString
AHR_single <- function(pb, linkage = "complete", partitionInit = NULL, # nolint: cyclocomp_linter
                       minNbClusters = 1L,
                       arretMin = FALSE, modeEvaluationDistInter = 2L,
                       calculDistComplet = FALSE,
                       fusionConstraint = FALSE,
                       fusionConstraintMode = "deterministic")
{
  assertLinkage(linkage)
  linkage <- corresponding_linkage(linkage)
  assertFusionConstraint(fusionConstraint)
  fusionConstraint <- corresponding_fusion_const(fusionConstraint)

  n <- pb$n()

  if (is.null(partitionInit))
    partitionInit <- seq_len(n)

  # Gestion de la contrainte et du mode de fusion
  penaliteTaille <- !is.na(fusionConstraint) &&
    fusionConstraint == "sizePenalty"

  if (fusionConstraint %in% c(NA_character_, "sizePenalty"))
    fusionConstraint <- FALSE

  else
  {
    typeContrainteFusion <- codeGroupesLiaison(fusionConstraint)

    if (typeContrainteFusion == codeGroupesLiaison("default"))
      fusionConstraint <- FALSE

    else
      fusionConstraint <- TRUE
  }

  if (fusionConstraint)
    fusionConstraintDeter <- fusionConstraintMode == "deterministic"

  # Initial partition data
  partition <- standardize_partition(partitionInit)
  nbClusters <- max(partition)
  nomClusters <- as.character(seq_len(nbClusters))
  listeClusters <- lapply(seq_len(nbClusters),
                          function(c) which(partition == c))
  names(listeClusters) <- nomClusters
  nbElementsClusters <- lengths(listeClusters)
  nbSingletons <- sum(nbElementsClusters == 1L)
  clustersRestants <- rep(TRUE, nbClusters)

  clustersSizes <- .clusters_sizes(partition, pb$sizes)

  # Check of the minimum size constraint
  nbClustersMin <- pb$nbTooSmallClusters(partition)
  contrainteMinVerifiee <- nbClustersMin == 0L

  # Definition of a "correct size" for a region
  tailleCorrecte <- 0.0
  if (pb$hasMinConstraint())
    tailleCorrecte <- pb$m

  if (pb$hasMaxConstraint())
    tailleCorrecte <- mean(c(tailleCorrecte, pb$M))

  # Look for available fusions
  contiguiteClasses <-
    pb$listeContiguiteClasses(partition, withFusionSize = TRUE)

  fusionSize <- contiguiteClasses[, "fusionSize"]
  contiguiteClasses <-
    contiguiteClasses[, colnames(contiguiteClasses) != "fusionSize",
                      drop = FALSE]

  # Creation of the contiguity graph with storage of linkage distances
  # and post-fusion sizes
  grapheClasses <-
    make_graph(edges = as.vector(t(contiguiteClasses)),
               n = nbClusters, directed = FALSE)

  V(grapheClasses)[as.numeric(names(clustersSizes))]$size <-
    unname(clustersSizes)

  E(grapheClasses)$distance <- NA_real_
  E(grapheClasses)$fusionSize <- fusionSize
  rm(clustersSizes, fusionSize)

  # Gestion de la distance de liaison
  distanceMinMax <- FALSE
  distancemedoid <- FALSE
  if (testString(linkage))
  {
    if (is_single_linkage(linkage))
    {
      distanceMinMax <- TRUE
      comp <- pmin
    }

    else if (is_complete_linkage(linkage))
    {
      distanceMinMax <- TRUE
      comp <- pmax
    }

    else if (is_centroid_linkage(linkage))
    {
      if (is.null(pb$data) || is.null(pb$d))
      {
        stop("With centroid linkage inter-elements distance `d` and context `data` must be given") # nolint: line_length_linter
      }

      if (!all(unlist(lapply(pb$data, is.numeric))))
      {
        stop("With centroid linkage all variables must be quantitative")
      }
    }

    else if (is_medoid_linkage(linkage))
    {
      distancemedoid <- TRUE
      V(grapheClasses)$medoid <- medoids_partition(pb, partition)
    }
  }

  # Function to obtain linkage distances with potential fusion constraint
  recup_distances_classes <- function()
  {
    # Distances order modification if a fusion constraint is activated
    fusionConstraintMade <- FALSE
    if (fusionConstraint)
    {
      coeffAlea <- 1.0 - nbClusters / n
      if (fusionConstraintDeter ||
          sample(c(FALSE, TRUE), size = 1L,
                 prob = c(1.0 - coeffAlea, coeffAlea)))
      {
        fusionConstraintMade <- TRUE
        masqueContraintesFusion <-
          transform_contiguity_list(pb,
                                    ends(grapheClasses, E(grapheClasses)),
                                    V(grapheClasses)$size,
                                    nbElementsClusters = nbElementsClusters,
                                    alea = FALSE,
                                    priorisation = FALSE,
                                    reduction = typeContrainteFusion,
                                    retour = "indices",
                                    byName = FALSE,
                                    orderByDistance = FALSE)


        if (!any(masqueContraintesFusion))
          return(double(0L))

        distances <- E(grapheClasses)[masqueContraintesFusion]$distance
        names(distances) <- masqueContraintesFusion
      }
    }

    if (!fusionConstraintMade)
    {
      distances <- E(grapheClasses)$distance

      if (length(distances) == 0L)
        return(double(0L))

      names(distances) <- masqueContraintesFusion <- seq_along(distances)
    }

    masque <- is.na(distances)
    # If some distances must be calculated
    if (any(masque))
    {
      aretesDistancesACalculer <- masqueContraintesFusion[masque]
      distancesACalculer <-
        ends(grapheClasses, aretesDistancesACalculer, names = FALSE)

      if (distanceMinMax && calculDistComplet)
      {
        distances[masque] <- distancesInter_mat[distancesACalculer]
      }

      else
      {
        if (distancemedoid)
        {
          distances[masque] <-
            .calcul_distances_inter(pb, partition, distancesACalculer,
                                    linkage, modeEval = modeEvaluationDistInter,
                                    nbElementsClusters = nbElementsClusters,
                                    listeClusters = listeClusters,
                                    medoids = V(grapheClasses)$medoid)
        }
        else
        {
          distances[masque] <-
            .calcul_distances_inter(pb, partition,
                                    distancesACalculer, linkage,
                                    modeEval = modeEvaluationDistInter,
                                    nbElementsClusters = nbElementsClusters,
                                    listeClusters = listeClusters)
        }
      }

      # Update of the distances in the contiguity graph
      assign("grapheClasses",
             set_edge_attr(grapheClasses, "distance",
                           index = aretesDistancesACalculer,
                           value = distances[masque]),
             envir = parent.frame())
    }

    # Add of a penality if fusionConstraint = "sizePenalty"
    if (penaliteTaille)
      distances <- distances *
      E(grapheClasses)$fusionSize[as.numeric(names(distances))]


    return(distances)
  }

  # Special cases of linkage initialisation
  # -- When linkage is single or complete and we allow us to calculate
  # all distances directly, even those with non-contiguous clusters
  if (distanceMinMax && calculDistComplet)
  {
    if (nbClusters == n)
    {
      ordre <- order(partition)
      distancesInter_mat <- pb[][ordre, ordre]
    }
    else
    {
      distancesInter_mat <- matrix(0.0, nrow = nbClusters, ncol = nbClusters)
      rownames(distancesInter_mat) <-
        colnames(distancesInter_mat) <-
        nomClusters

      grille <- cbind(rep(seq_len(nbClusters), each = nbClusters),
                      rep(seq_len(nbClusters), times = nbClusters))
      grille <- grille[grille[, 1L] < grille[, 2L], ]
      distancesInter_mat[grille] <- distancesInter_mat[grille[, c(2L, 1L)]] <-
        calculDistancesInter_Cpp(pb, partition, grille, linkage)
    }
  }


  # Definition of a partition model for recording in `AHCTree`
  contraintes <-         c(TRUE, FALSE,  TRUE)
  names(contraintes) <- c("connectivity", "min", "max")

  if (contrainteMinVerifiee)
    scoreMinSizeConst <- 0.0
  else
  {
    scoreMinSizeConst <-
      score_constraints_min(pb$m, clustersSizes = V(grapheClasses)$size)
  }




  partitionRef <-
    constructor_Partition(pb, partition = partition, method = "RAH",
                          contiguity = TRUE,
                          minConstraint = contrainteMinVerifiee,
                          maxConstraint = TRUE,
                          scoreMinSizeConst = scoreMinSizeConst,
                          scoreMaxSizeConst = 0.0)


  # Classe ArbreCAH stockant les partitions
  arbreCAH <- constructor_AHCTree(pb, partitionRef)
  rm(partitionRef)

  feasible <- TRUE
  nbClusters <- nbClusters - 1L

  while (!(arretMin && contrainteMinVerifiee) &&
        feasible &&
        nbClusters >= minNbClusters)
  {
    distancesClasses_vec <- recup_distances_classes()

    if (length(distancesClasses_vec) == 0L)
    {
      feasible <- FALSE
      break
    }

    nomsDistances <- names(distancesClasses_vec)
    posDistanceMin <- which.min(distancesClasses_vec)[1L]

    areteAssociee_int <- as.numeric(nomsDistances[posDistanceMin])
    clusters <- ends(grapheClasses, areteAssociee_int, names = TRUE)

    cluster1 <- clusters[1L, 1L]
    tailleCluster1 <- V(grapheClasses)[cluster1]$size

    cluster2 <- clusters[1L, 2L]
    tailleCluster2 <- V(grapheClasses)[cluster2]$size

    if (nbSingletons > 0L)
    {
      cluster1EstSingleton <- nbElementsClusters[cluster1] == 1L
      cluster2EstSingleton <- nbElementsClusters[cluster2] == 1L

      nbSingletons <- nbSingletons - cluster1EstSingleton - cluster2EstSingleton
      if (fusionConstraint && nbSingletons == 0L)
        typeContrainteFusion <- max(codeGroupesLiaison("minMin"),
                                    typeContrainteFusion)
    }
    if (nbClustersMin > 0L)
    {
      cluster1EstMin <- tailleCluster1 < pb$m
      cluster2EstMin <- tailleCluster2 < pb$m

      # If at least one of the clusters are a "min", we lose at least
      # one min cluster
      if (cluster1EstMin || cluster2EstMin)
        nbClustersMin <- nbClustersMin - 1L

      # If both are "min" and fusion implies a big enough cluster we loose
      # two "min" clusters
      if (cluster1EstMin &&
          cluster2EstMin &&
          tailleCluster1 + tailleCluster2 >= pb$m)
        nbClustersMin <- nbClustersMin - 1L

      if (nbClustersMin == 0L)
      {
        fusionConstraint <- FALSE
        contrainteMinVerifiee <- TRUE
      }
    }

    # Partition update
    partition[partition == cluster2] <- cluster1
    listeClusters[[cluster1]] <- c(listeClusters[[cluster1]],
                                   listeClusters[[cluster2]])
    listeClusters[[cluster2]] <- numeric(0L)
    clustersRestants[cluster2] <- FALSE

    # Medoid update (if actual linkage distance is medoid)
    if (distancemedoid)
    {
      elementsCluster1 <- listeClusters[[cluster1]]
      grapheClasses <-
        set_vertex_attr(grapheClasses, "medoid", cluster1,
                        elementsCluster1[medoid(pb[elementsCluster1,
                                                   elementsCluster1,
                                                   drop = FALSE])])
    }

    partitionArbre <- arbreCAH[[nbClusters]]
    partitionArbre$partition <- standardize_partition(partition)
    partitionArbre$contraintes["min"] <- contrainteMinVerifiee

    if (contrainteMinVerifiee)
      partitionArbre$scoreSizeConsts[c("min", "total")] <- 0.0
    else
    {
      partitionArbre$scoreSizeConsts[c("min", "total")] <-
        score_constraints_min(pb$m, clustersSizes = V(grapheClasses)$size)
    }

    # Add of the new contiguities to the merged cluster
    voisinsCluster1 <- as.numeric(neighbors(grapheClasses, cluster1))
    voisinsCluster1 <- voisinsCluster1[voisinsCluster1 != cluster2]
    voisinsCluster2 <- as.numeric(neighbors(grapheClasses, cluster2))
    voisinsCluster2 <- voisinsCluster2[voisinsCluster2 != cluster1]
    nouvellesConnexites <- setdiff(voisinsCluster2, voisinsCluster1)
    if (length(nouvellesConnexites) > 0L)
    {
      fusionSizeNouvellesConnexites <- tailleCluster1 + tailleCluster2 +
        get.vertex.attribute(grapheClasses, "size", nouvellesConnexites)

      # We don't add contiguities that would imply a disrespect of the
      # maximum size constraint
      if (pb$hasMaxConstraint())
      {
        masqueAjout <- fusionSizeNouvellesConnexites <= pb$M
        nouvellesConnexites <- nouvellesConnexites[masqueAjout]
        fusionSizeNouvellesConnexites <-
          fusionSizeNouvellesConnexites[masqueAjout]
      }

      if (length(nouvellesConnexites) > 0L)
      {
        nouvellesAretes <- rep(nouvellesConnexites, each = 2L)
        nouvellesAretes[2L * seq_along(nouvellesConnexites)] <- cluster1

        grapheClasses <-
          add_edges(grapheClasses, nouvellesAretes,
                    attr = list(distance = NA_real_,
                                fusionSize = fusionSizeNouvellesConnexites))
      }
    }

    # Sizes update
    if (length(voisinsCluster1) > 0L)
    {
      listeVoisins1_vec <- rep(voisinsCluster1, each = 2L)
      listeVoisins1_vec[2L * seq_along(voisinsCluster1)] <- cluster1
      idsAretesCluster1 <- get.edge.ids(grapheClasses, listeVoisins1_vec)
      nouvellesTaillesFusion <-
        get.edge.attribute(grapheClasses, "fusionSize", idsAretesCluster1) +
        tailleCluster2

      # If some fusions would imply a disrespect of the maximum size constraint
      # they are removed
      if (pb$hasMaxConstraint())
        masqueModif <- nouvellesTaillesFusion <= pb$M

      else
        masqueModif <- rep(TRUE, length(nouvellesTaillesFusion))

      masqueSuppression <- !masqueModif
      if (any(masqueModif))
      {
        grapheClasses <-
          set_edge_attr(grapheClasses, "fusionSize",
                        idsAretesCluster1[masqueModif],
                        nouvellesTaillesFusion[masqueModif])
      }
      if (any(masqueSuppression))
      {
        grapheClasses <- delete_edges(grapheClasses,
                                      idsAretesCluster1[masqueSuppression])
        voisinsCluster1 <- voisinsCluster1[masqueModif]
      }
    }


    V(grapheClasses)[c(cluster1, cluster2)]$size <-
      c(tailleCluster1 + tailleCluster2, tailleCorrecte)
    nbElementsClusters[cluster1] <- nbElementsClusters[cluster1] +
                                    nbElementsClusters[cluster2]
    nbElementsClusters[cluster2] <- 0L

    liaisonsAEffacer <- voisinsCluster1


    # Special cases distance updates
    if (distanceMinMax)
    {
      if (calculDistComplet)
      {
        clustersRestants[cluster1] <- FALSE
        distancesInter_mat[cluster1, clustersRestants] <-
          distancesInter_mat[clustersRestants, cluster1] <-
          comp(distancesInter_mat[cluster1, clustersRestants, drop = TRUE],
               distancesInter_mat[cluster2, clustersRestants, drop = TRUE])
        clustersRestants[cluster1] <- TRUE
        distancesInter_mat[cluster2, ] <- distancesInter_mat[, cluster2] <- Inf

        ensembleVoisins1 <- c(voisinsCluster1, nouvellesConnexites)
        if (length(ensembleVoisins1) > 0L)
        {
          listeNoeuds1_vec <- rep(ensembleVoisins1, each = 2L)
          listeNoeuds1_vec[2L * seq_along(ensembleVoisins1)] <- cluster1
          idsAretesCluster1 <- get.edge.ids(grapheClasses, listeNoeuds1_vec)

          nouvellesDistances <-
            unname(distancesInter_mat[cluster1, ensembleVoisins1, drop = TRUE])

          grapheClasses <- set_edge_attr(grapheClasses, "distance",
                                         index = idsAretesCluster1,
                                         value = nouvellesDistances)
        }

        liaisonsAEffacer <- NULL
      }

      else
      {
        voisinageCommun <- intersect(voisinsCluster1, voisinsCluster2)
        if (length(voisinageCommun) > 0L)
        {

          listeNoeuds1_vec <- rep(voisinageCommun, each = 2L)
          listeNoeuds1_vec[2L * seq_along(voisinageCommun)] <- cluster1
          idsAretesCluster1 <- get.edge.ids(grapheClasses, listeNoeuds1_vec)


          distanceACluster1 <-
            get.edge.attribute(grapheClasses, "distance", idsAretesCluster1)

          listeNoeuds2_vec <- rep(voisinageCommun, each = 2L)
          listeNoeuds2_vec[2L * seq_along(voisinageCommun)] <- cluster2
          idsAretesCluster2 <- get.edge.ids(grapheClasses, listeNoeuds2_vec)


          distancesACluster2 <-
            get.edge.attribute(grapheClasses, "distance", idsAretesCluster2)

          grapheClasses <-
            set_edge_attr(grapheClasses, "distance", idsAretesCluster1,
                          comp(distanceACluster1, distancesACluster2))

          # Reinitialisation of the cluster1-* distances determined before
          # fusion that have not been updated
          liaisonsAEffacer <-
            setdiff(voisinsCluster1, c(voisinageCommun, cluster2))
        }
      }
    }


    # Reinitialisation of some distances
    if (length(liaisonsAEffacer) > 0L)
    {
      listeNoeuds_vec <- rep(liaisonsAEffacer, each = 2L)
      listeNoeuds_vec[2L * seq_along(liaisonsAEffacer)] <- cluster1
      grapheClasses <-
        set_edge_attr(grapheClasses, "distance",
                      index = get.edge.ids(grapheClasses, listeNoeuds_vec),
                      value = NA_real_)
    }

    # second cluster is removed from the graph
    voisinsCluster2 <- c(voisinsCluster2, cluster1)
    listeAreteSupprimer <- rep(voisinsCluster2, each = 2L)
    listeAreteSupprimer[2L * seq_along(voisinsCluster2)] <- cluster2
    grapheClasses <-
      delete_edges(grapheClasses, get.edge.ids(grapheClasses,
                                               listeAreteSupprimer))

    nbClusters <- nbClusters - 1L

    # If there are no fusion available loop is stopped
    if (gsize(grapheClasses) == 0L)
      feasible <- FALSE
  }

  return(arbreCAH[!(names(arbreCAH) %in% 1L:nbClusters)])
}
