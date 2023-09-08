set.seed(123L)
grid  <- simple_grid(4L, 5L)
data <- grid$context
sizes <- grid$repartition$nbIndividuals
contiguity <- grid$contiguity

m <- 200.0
M <- 800.0

d <- "euclidean"

criterion <- 'CHI'

AHR(contiguity = contiguity,
    d = d, data = data,
    sizes = sizes,
    m = m, M = M,
    criteria = criterion,
    nbTries = 1L,
    fusionConstraints = available_fusion_constraints(),
    fusionConstraintModes = available_fusion_modes(),
    parallele = FALSE)
