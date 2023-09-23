p1 <- c(2, 4, 3)
p2 <- c("yes", "no", "maybe")
equivalent_partitions(p1, p2) # TRUE

p3 <- c(2, 4, 3)
p4 <- c("yes", "no", "no")
equivalent_partitions(p3, p4) # FALSE
