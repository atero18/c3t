contiguity1 <- matrix(c(FALSE, TRUE, FALSE,
                        TRUE, FALSE, TRUE,
                        FALSE, TRUE, FALSE), nrow = 3L)

is_connected(contiguity1) # TRUE

contiguity2 <- matrix(c(FALSE, FALSE,
                        FALSE, FALSE), nrow = 2L)

is_connected(contiguity2) # FALSE
