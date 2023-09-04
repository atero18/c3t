source("tools/profiling/setup_profiling.R") # nolint
# nolint start: undesirable_function_linter
library(tibble)
library(dplyr)
library(ggplot2)
# nolint end
loadNamespace("fpc")

# Profiling sur un petit nombre de données
# -- Cas avec indice uniquement
data <- c3t_grid_simulation(30L, 50L)$data
partition <- sample(seq_len(nrow(data)), replace = TRUE, size = nrow(data))

profil <- profvis(calinski_harabasz(data, partition, valueOnly = TRUE))
profil

# -- Cas avec retour global
profil2 <- profvis(calinski_harabasz(data, partition, valueOnly = FALSE))
profil2

# Comparaison des performances entre fpc::calinhara et c3t



data <- c3t_grid_simulation(30L, 50L)$data

performances_ch <- function(calculCH)
{

  kMax <- ceiling(nrow(data) / 2.0)


  set.seed(123L)
  nbTries <- 40L

  res <- double(nbTries)
  while (nbTries > 0L)
  {
    n <- sample(seq_len(nrow(data)), size = 1L)
    partition <- sample(seq_len(min(n, kMax)), size = n, replace = TRUE)
    partition <- standardize_partition(partition)
    calculCH(data[seq_len(n), ], partition)
    nbTries <- nbTries - 1L
  }

  res
}


m <- mark(fpc = performances_ch(fpc::calinhara),
          c3t = performances_ch(c3t:::calinski_harabasz)) # nolint



data("grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0")
data <- grid_queen_vide0_metropole0_x50_y70_indivMoy100_quant3_qual0$context

performances_ch2 <- function(calculICH)
{
  nValues <- seq(10L, nrow(data), 10L)
  nbTriesPerN <- 10L

  minK <- 10L

  res <- tibble(n = rep(nValues, each = nbTriesPerN),
                k = NA_integer_,
                time = NA_real_)

  set.seed(123L)
  pos <- 1L
  for (n in nValues)
  {
    actualData <- data[seq_len(n), ]

    kValues <- ceiling(seq(from = ifelse(minK < n - 1L, minK, 1L), to = n - 1L,
                           length.out = nbTriesPerN))
    for (k in kValues)
    {

      partition <- c(sample(seq_len(k)),
                     sample(seq_len(k), replace = TRUE, size = n - k))


      t0 <- Sys.time()
      calculICH(actualData, partition)
      res[pos, "time"] <- as.numeric(Sys.time() - t0, units = "secs")
      res[pos, "k"] <- k

      pos <- pos + 1L

    }
  }

  res
}

resFPC <- performances_ch2(fpc::calinhara)
resFPC$package <- "fpc"
resC3T <- performances_ch2(c3t:::calinski_harabasz) # nolint
resC3T$package <- "c3t"

comp <- rbind(resFPC,
              resC3T)

rm(resFPC, resC3T)

comp %>%
  group_by(package, n) %>%
  summarize(medTimeN = median(time)) %>%
  ggplot(aes(x = n, y = medTimeN)) +
  geom_point(aes(colour = package)) +
  ggtitle("Time complexity (s) depending of the number of elements n")

comp %>%
  group_by(n, k) %>%
  mutate(time = time / min(time)) %>%
  ungroup() %>%
  group_by(package, n) %>%
  summarize(medTimeN = median(time)) %>%
  ggplot(aes(x = n, y = medTimeN)) +
  geom_point(aes(colour = package)) +
  ggtitle("Time complexity (s) depending of the number of elements n")


performances_ch3 <- function()
{
  nValues <- seq(10L, nrow(data), 10L)

  res <- tibble(n = rep(nValues, each = 2L),
                k = ceiling(rep(nValues, each = 2L) / 2L),
                package = rep(c("fpc", "c3t"), times = length(nValues)))

  memory <- NULL

  set.seed(123L)
  pos <- 1L
  for (n in nValues)
  {
    actualData <- data[seq_len(n), ]

    k <- ceiling(n / 2L)

    partition <- c(sample(seq_len(k)),
                   sample(seq_len(k), replace = TRUE, size = n - k))


    m <- mark(fpc = fpc::calinhara(actualData, partition),
              c3t = c3t:::calinski_harabasz(actualData, partition), # nolint
              iterations = 4L)

    memory <- c(memory, m$mem_alloc / m$n_itr)

    pos <- pos + 2L

  }

  cbind(res, memory)
}

res <- performances_ch3()

res %>%
  ggplot(aes(x = n, y = memory)) +
  geom_point(aes(colour = package)) +
  ggtitle("Memory complexity (bytes) depending of the number of elements n")

ratio <- tibble(n = res$n[seq(1L, nrow(res), 2L)],
                memory = res[res$package == "fpc", "memory"] /
                  res[res$package == "c3t", "memory"])

ggplot(ratio, aes(x = n, y = memory)) + geom_line()
