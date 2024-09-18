## code to prepare `capeverde_pops` dataset
capeverde_pops <- c(rep("Mandenka",22), rep("Gambian",109),
                    rep("Cape Verdean",44), rep("Iberian",107),
                    rep("French",28), rep("British",89))
usethis::use_data(capeverde_pops, overwrite = TRUE)
