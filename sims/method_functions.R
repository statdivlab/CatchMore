## @knitr methods


catchall_best <- new_method("catchall", "Best model with CatchAll",
                        method = function(model, draw) {
                          bmm <- CatchAll::best_model(draw)
                          list(Model = bmm$model,
                               obs = sum(counts > 0),
                               est = bmm$est)
                        })

