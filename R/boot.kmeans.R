boot.kmeans <- function(data = NULL,
                  groups = NULL,
                  iterations=500,
                  nstart=1,
                  export=FALSE,
                  display=FALSE,
                  pval=0.05,
                  itermax=10,
                  maxsamp=1000,
                  verbose=FALSE,
                  returnall=FALSE){

  if(is.null(data)) stop("argument 'data' is missing with no default")
  if(is.null(groups)) stop("argument 'groups' is missing with no default: specify number of clusters")

  if(is.matrix(groups)){
    centgiven = TRUE
    centsgiven <- groups
    groups <- nrow(centsgiven)
  }
  else{ centgiven <- FALSE }
  centers <- list()
  bootmat <- cbind()
  ooblist <- kmlist <- obsdistlist <- list()
  oobmat <- obsdist <- matrix(0, nrow=nrow(data), ncol=groups)
  soslist <- c()
  p.value <- -1
  inititer <- iterations

  if(verbose) cat("Running bootstrapped samples\n")
  i = 0
  straightData <- FALSE
  while(!straightData) {
    i=i+1
    if(verbose) cat("\r", "Iterations:", i, ifelse(p.value != -1, paste("p.value:",p.value), ""))
    if(i==1) {
      sampl <- sort(sample(nrow(data),replace=TRUE))
      newData <- data[sampl,]
      if(centgiven){
        kmlist[[i]] <- km <- kmeans(newData,centers=centsgiven, nstart=nstart, iter.max=itermax)
      }
      else{
        kmlist[[i]] <- km <- kmeans(newData,centers=groups, nstart=nstart, iter.max=itermax)
      }
    }
    else {
      sampl <- sort(sample(nrow(data),replace=TRUE))
      newData <- data[sampl,]
      kmlist[[i]] <- km <- kmeans(newData,centers=km$centers, iter.max = itermax)
    }
    ooblist[[i]] <- setdiff(1:nrow(data),sampl)
    centers[[i]] <- km$centers

    for(j in 1:groups) obsdist[,j] <- mahalanobis(data, km$centers[j,], diag(ncol(data)), inverted=TRUE)

    obsdistlist[[i]] <- obsdist
    allocations <- apply(obsdist, 1, which.min)
    soslist[i] <- sum(apply(obsdist, 1, min))
    bootmat <- cbind(bootmat, allocations)

    if(display == TRUE) plot(tail(soslist, inititer))
    if(i>=inititer) {
      p.value = bgtest(tail(soslist, inititer) ~ c(1:inititer))$p.value
      if(p.value<pval & iterations < maxsamp) {
        iterations = iterations +1
      }
      else {
        straightData = TRUE
        if(verbose) cat("\r", "Iterations:", i, ifelse(p.value != -1, paste("p.value:",p.value), ""))
      }
    }
    if(export) {
      jpeg(paste("plot",i,".jpg", sep=''))
      plot(tail(soslist, inititer))
      dev.off()
    }
  }

  convcenters <- apply(abind(centers[(iterations-inititer+1):iterations], along=3), c(1,2), mean)
  oobsind <- as.data.frame(lapply(ooblist, function(v) replace(rep(0,nrow(data)),v,1)), col.names=1:iterations)
  oobcounts <- oobsind[,(iterations-inititer+1):iterations] * bootmat[,(iterations-inititer+1):iterations]
  bigoobz <- matrix(0, nrow=nrow(data), ncol=groups)

  for(j in 1:groups) bigoobz[,j] <-  apply(oobcounts, 1, function(v) sum(v==j))

  fuzzoob <- bigoobz/rowSums(bigoobz)
  hardoob <- apply(fuzzoob, 1, which.max)
  if(returnall == TRUE){
    value <- list(U = fuzzoob,
                clusters = hardoob,
                centers = convcenters,
                centerlist = centers,
                ooblist = ooblist,
                p.value = p.value,
                iterations= iterations,
                occurences = bootmat,
                size = groups,
                soslist = soslist,
                kmlist = kmlist)
    attr(value, "class") <- "BSKMeans"
  } else {
    value <- list(U = fuzzoob,
                  clusters = hardoob,
                  centers = convcenters,
                  centerlist = NULL,
                  ooblist = NULL,
                  p.value = p.value,
                  iterations= iterations,
                  occurences = bootmat,
                  size = groups,
                  soslist = soslist,
                  kmlist = NULL)
    attr(value, "class") <- "BSKMeans"
  }
  return(value)
}



