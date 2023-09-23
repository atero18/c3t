set.seed(123L)
grid  <- simple_grid(4L, 5L)
data <- grid$context
sizes <- grid$repartition$nbIndividuals
contiguity <- grid$contiguity

m <- 200.0
M <- 800.0

d <- "euclidean"

unfeasibleSolution <-
  c(1L, 2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 4L,
    4L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)

resolve_unfeasible(contiguity = contiguity,
                   sizes = sizes,
                   data = data,
                   d = d, m = m, M = M,
                   regionalisation = unfeasibleSolution,
                   verbose = FALSE)
