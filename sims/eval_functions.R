## @knitr metrics

bias_richness <- new_metric("bias", "Bias",
                        metric = function(model, out) {
                          return(out$est - model$CC)
})

se_richness <- new_metric("SE", "S.E.",
                        metric = function(model, out) {
                          return(sd(out$est))
                        })
