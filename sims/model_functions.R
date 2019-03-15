## @knitr models

## Count data follow a Poisson distribution
poisson_counts <- function(CC, lambda) {
  new_model(name = "Poisson-Model",
            label = sprintf("Poisson Model (C = %s, lambda = %s)", CC, lambda),
            params = list(CC = CC, lambda = lambda),
            simulate = function(CC, lambda, nsim) {
              counts_list <-vector(mode = "list", length = nsim)
              for (ii in 1:nsim) {
                counts <- rpois(CC, lambda)
                
                ## generate frequency count table
                fc_tab <- as.data.frame(table(counts))
                names(fc_tab) <- c("index", "frequency")
                fc_tab$index <- as.numeric(fc_tab$index)
                fc_tab$frequency <- as.numeric(fc_tab$frequency)
                fc_tab <- fc_tab[fc_tab$index > 0, ]
                counts_list[[ii]] <- fc_tab
              } 
              return(counts_list) # make each col its own list element
            })
}

## Count data follow a geometric distribution
geom_counts <- function(CC, lambda) {
  new_model(name = "Geometric-Model",
            label = sprintf("Geometric Model (C = %s, prob = %s)", CC, prob),
            params = list(CC = CC, prob = prob),
            simulate = function(CC, prob, nsim) {
              counts_list <-vector(mode = "list", length = nsim)
              for (ii in 1:nsim) {
                counts <- rgeom(CC, prob)
                
                ## generate frequency count table
                fc_tab <- as.data.frame(table(counts))
                names(fc_tab) <- c("index", "frequency")
                fc_tab$index <- as.numeric(fc_tab$index)
                fc_tab$frequency <- as.numeric(fc_tab$frequency)
                fc_tab <- fc_tab[fc_tab$index > 0, ]
                counts_list[[ii]] <- fc_tab
              } 
              return(counts_list) # make each col its own list element
            })
}