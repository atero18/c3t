contiguity <- simple_grid(3L, 2L)$contiguity

regionalisation <- c(1L, 1L, 1L,
                     2L, 2L, 2L)

checkRegionalisation(regionalisation, contiguity) # TRUE

notRegionalisation <- c(1L, 2L, 1L,
                        1L, 2L, 1L)

checkRegionalisation(notRegionalisation, contiguity) # message

