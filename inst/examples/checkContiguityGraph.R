require(igraph)

graph <- make_graph(c(1, 2, 2, 3, 3, 1), directed = FALSE)

checkContiguityGraph(graph) # TRUE
