compare.tables <- function(full.res = NULL,
                           true.labs = NULL){
  if(is.null(full.res)) stop("argument full.res missing with no default: input list generated from function compare.clusters")
  if(is.null(true.labs)) stop("argument true.labs missing with no default: input vector of true class labels")
  cat("Kmeans")
  print(table(true.labs, full.res$km$cluster))
  cat("\n")
  cat("BootKmeans MAP-oob")
  print(table(true.labs, full.res$bkm$clusters))
  if(full.res$what=="all"){
    cat("\n")
    cat("FuzzKmeans MAP")
    print(table(true.labs, full.res$fkm$clus[,1]))
  }
}
