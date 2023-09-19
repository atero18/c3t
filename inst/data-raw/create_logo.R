library(igraph)
library(hexSticker)
library(ggplot2)
library(cowplot)
library(magick)

emptyLetter <- matrix(FALSE, ncol = 3L, nrow = 3L)


contiguityGraphC3T <- matrix(NA_character_)



# "c"

edgesC <- c(1, 2, 2, 3,
            1,
            4,
            4,
            5,
            5, 6, 6, 7)

nbElementsC <- 7L

graphC <- make_graph(edgesC,
                    directed = FALSE)

layoutGraphC <- matrix(c(0L, 0L, 1L, 0L, 2L, 0L,
                         0L, 1L,
                         0L, 2L, 1L, 2L, 2L, 2L), ncol = 2L, byrow = TRUE)

plot.igraph(graphC, layout = layoutGraphC)

V(graphC)$name <- seq_len(gorder(graphC))

# "3"

edges3 <- c(1, 2, 2, 3,
            3,
            4,
            5, 4,
            4,
            6,
            8, 7, 7, 6)

nbElements3 <- 8L
graph3 <- make_graph(edges3, directed = FALSE)

layoutGraph3 <- matrix(c(0L, 0L,
                         0L, 1L,
                         0L, 2L,
                         1L, 2L,
                         1L, 1L,
                         2L, 2L,
                         2L, 1L,
                         2L, 0L), byrow = TRUE, ncol = 2L)
layoutGraph3 <- layoutGraph3[, 2L:1L]


plot.igraph(graph3, layout = layoutGraph3)

V(graph3)$name <- gorder(graphC) + seq_len(gorder(graph3))

# "t"

edgesT <- c(1,
            2,
            2,3,
            2,
            4)

nbElementsT <- 4L
graphT <- make_graph(edgesT, directed = FALSE)

layoutGraphT <- matrix(c(0.0,  0,
                         0.0, -0.8,
                         3.0, -0.8,
                         0.0, -2),
                         byrow = TRUE, ncol = 2L)


edges <- c(edgesC, edges3 + max(edgesC), edgesT + max(edgesC) + max(edges3))
c3tGraph <- make_graph(edges)

layout <- rbind(layoutGraphC,
                layoutGraph3,
                layoutGraphT)
layout[seq(nbElementsC, length.out = nbElements3), 1L] <-
  layout[seq(nbElementsC, length.out = nbElements3), 1L] + max(layoutGraphT[,1L])

layout[seq(nbElementsC + nbElements3, length.out = nbElementsT)] <-
  layout[seq(nbElementsC + nbElements3, length.out = nbElementsT), 1L] +
  max(layoutGraphC[, 1L]) + max(layoutGraph3[, 1L])

plot.igraph(c3tGraph, layout = layout)

p <- ggdraw() +
  draw_image("C:/Users/tangi/OneDrive - GENES/Année 2 (2022-2023)/Stage/Logos/2023-09-19 - Logo avec c3t et clustering, logomaker.ai.png")


# Inspired from
# https://shixiangwang.github.io/home/en/post/2019-06-20-how-i-create-ucscxenatools-logo/
sticker(p, package = "c3t",
        url = "atero18.github.io/c3t/",
        filename = "man/figures/logo.png",
        p_size = 0L, s_x = 1.0, s_y = 1.0, s_width = 2.5, s_height = 2.5,
        p_x = 1.1, p_y = 0.9, h_fill = "white", h_color = "black",
        u_x = 1.25, u_y = 0.25, u_size = 3, u_color = "blue")
