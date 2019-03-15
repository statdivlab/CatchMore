setwd("P:/git/CatchAll/misc")
library(CatchAll)

poisson_sim_once <- function(CC = 500, lambda = 3) {
  counts <- rpois(CC, lambda)

  ## generate frequency count table
  fc_tab <- as.data.frame(table(counts))
  names(fc_tab) <- c("index", "frequency")
  fc_tab$index <- as.numeric(fc_tab$index)
  fc_tab$frequency <- as.numeric(fc_tab$frequency)
  fc_tab <- fc_tab[fc_tab$index > 0, ]
  bmm <- best_model(fc_tab)
  return(list(Model = bmm$model,
              obs = sum(counts > 0),
              est = bmm$est))
}


geom_sim_once <- function(CC = 500, prob = 0.2) {
  counts <- rgeom(CC, prob)

  ## generate frequency count table
  fc_tab <- as.data.frame(table(counts))
  names(fc_tab) <- c("index", "frequency")
  fc_tab$index <- as.numeric(fc_tab$index)
  fc_tab$frequency <- as.numeric(fc_tab$frequency)
  fc_tab <- fc_tab[fc_tab$index > 0, ]
  bmm <- best_model(fc_tab)
  return(list(Model = bmm$model,
              obs = sum(counts > 0),
              est = bmm$est))
}

for(CC in seq(100, 1000, by = 100)) {
  pois_results <- replicate(100, poisson_sim_once(CC = CC, lambda = 3))
  geom_results <- replicate(100, geom_sim_once(CC = CC, prob = 0.2))
  save(pois_results, file = paste("./output/pois_C_", CC, ".RData", sep = ""))
  save(geom_results, file = paste("./output/geom_C_", CC, ".RData", sep = ""))
}
