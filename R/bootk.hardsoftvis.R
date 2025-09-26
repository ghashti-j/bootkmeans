bootk.hardsoftvis <- function(data = NULL, res, plotallvars = FALSE, var1 = NULL, var2 = NULL){
  if(is.null(data)) stop("Argument data empty: input dataset")
  if(is.null(res)) stop("Input a boot.kmeans result list (BSKMeans)")
  col <- numeric()
  for(i in 1:nrow(data)) ifelse(any(res$U[i,] == 1), col[i] <- 1, col[i] <- 2)
  if(plotallvars){
    op <- par(no.readonly = TRUE)
    pairs(data, col = col + 2, oma=c(3,3,3,10))
    par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend("right", legend = c("Hard", "Soft"), col = c(3, 4), pch = 19, xpd=TRUE)
    par(op)
  } else{
    if(is.null(var1) || is.null(var2)) stop("Specify two numeric values for var1 and var2 to plot")
    plot(data[,var1], data[,var2], col = col+2, xlab = names(data)[var1], ylab = names(data)[var2])
    legend("topright", legend = c("Hard", "Soft"), col = c(3, 4), pch = 19)
  }
}
