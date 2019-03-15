setwd("P:/git/CatchAll/misc")
library(dplyr)

geom_grid <- data.frame(C = rep(1:10 * 100, each = 100),
                        C_hat = NA,
                        OBS = NA,
                        correct_model = NA)

for(CC in (1:10 * 100)) {
  load(paste("./output/geom_C_", CC, ".RData", sep = ""))

  models <- sapply(1:100, function(ii) geom_results[1,][[ii]])
  obs <- sapply(1:100, function(ii) geom_results[2,][[ii]])
  ests <- sapply(1:100, function(ii) geom_results[3,][[ii]])

  geom_grid$C_hat[geom_grid$C == CC] <- ests
  geom_grid$OBS[geom_grid$C == CC] <- obs
  geom_grid$correct_model[geom_grid$C == CC] <- models == "Geometric"
}

pdf("geom_boxplots.pdf")
boxplot(C_hat ~ C, data = geom_grid, main = "Geometric distributed counts")
points(x = 1:10 * 100, y = 1:10 * 100, col = "red", pch = 4)
dev.off()

mean_C_hat <- geom_grid %>%
  group_by(C) %>%
  summarise(C_AVE = mean(C_hat))

save(mean_C_hat, file = "mean_C_geom.RData")

pois_grid <- data.frame(C = rep(1:10 * 100, each = 100),
                        C_hat = NA,
                        OBS = NA,
                        correct_model = NA)

for(CC in (1:10 * 100)) {
  load(paste("./output/pois_C_", CC, ".RData", sep = ""))

  models <- sapply(1:100, function(ii) pois_results[1,][[ii]])
  obs <- sapply(1:100, function(ii) pois_results[2,][[ii]])
  ests <- sapply(1:100, function(ii) pois_results[3,][[ii]])

  pois_grid$C_hat[pois_grid$C == CC] <- ests
  pois_grid$OBS[pois_grid$C == CC] <- obs
  pois_grid$correct_model[pois_grid$C == CC] <- models == "poisetric"
}

pdf("poisson_counts.pdf")
boxplot(C_hat ~ C, data = pois_grid, main = "Poisson distributed counts")
points(x = 1:10 * 100, y = 1:10 * 100, col = "red")
dev.off()

mean_C_hat <- pois_grid %>%
  group_by(C) %>%
  summarise(C_AVE = mean(C_hat))

save(mean_C_hat, file = "mean_C_pois.RData")
