## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE,warning = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(stats)
library(MASS)
library(bootkmeans)
library(abind)
library(fclust)
library(lmtest)
library(ggplot2)
library(patchwork)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# library(devtools)
# install_github("https://github.com/ghashti-j/bootkmeans")
# library(kdml)

## -----------------------------------------------------------------------------
set.seed(1)
x <- as.matrix(iris[, -5])

fit <- boot.kmeans(
  data = x, 
  groups = 3,
  iterations = 500,   
  itermax   = 1000,
  nstart    = 5,
  verbose   = FALSE,
  maxsamp = 1000
)

## -----------------------------------------------------------------------------
cat("bootkmeans centres:\n")
fit$centers

cat("true cluster centers")
aggregate(. ~ Species, data = iris, FUN = mean)

## -----------------------------------------------------------------------------
fit$U[1:10, ]

## -----------------------------------------------------------------------------
fit$U[51:60,]

## -----------------------------------------------------------------------------
cat("\n p-value for BG-test: \n")
fit$p.value

## ----fig.align='center', fig.width=8,fig.height=8-----------------------------
bootk.hardsoftvis(x, fit, plotallvars = TRUE)

## -----------------------------------------------------------------------------
res <- compare.clusters(
  data = x, 
  groups = 3, 
  nstart = 10,
  what =  "all"
)

compare.tables(res, true.labs = iris$Species)

## -----------------------------------------------------------------------------
set.seed(1)
numClusters <- 10
numObsPerCluster <- 30
radius <- 6
dim <- 2
centreNumObs <- numObsPerCluster
angles <- seq(0, 2 * pi, length.out = numClusters)

outerMean <- t(sapply(1:(numClusters - 1),
                      function(i) c(radius * cos(angles[i]), 
                                    radius * sin(angles[i]))))
centreMean <- c(0, 0)
allMeans <- rbind(outerMean, centreMean)
covMat <- diag(1, dim)

data <- do.call(rbind, 
                c(lapply(1:(numClusters - 1), 
                         function(i) MASS::mvrnorm(numObsPerCluster, 
                                                   outerMean[i, ], covMat)),
                  list(MASS::mvrnorm(centreNumObs, centreMean, covMat))
))

# within cluster density calculation
densities <- sapply(1:numClusters, 
                    function(i) mvtnorm::dmvnorm(data, mean = allMeans[i, ], 
                                                 sigma = covMat))

totDensity <- rowSums(densities) # overall densities
U <- densities / totDensity # probabilistic (fuzzy) cluster assignments
hardU <- apply(U, 1, which.max) # assigning hard cluster labels 

## ----fig.align='center', fig.width=8,fig.height=8-----------------------------
plotDF <- data.frame(X = data[,1],
                     Y = data[,2],
                     U = 1- apply(U, 1, max),
                     Z = factor(hardU))

ggplot(plotDF, aes(x = X, y = Y, col = Z, size = U)) +
  geom_point() +
  theme_classic() +
  coord_equal() +
  labs(x = expression(X[1]), 
         y = expression(X[2]), 
         col = "Hard Cluster",
         size = "Uncertainty") +
  theme(panel.border = element_rect(NA, "black", 1))

## ----fig.align='center', fig.width=8,fig.height=8-----------------------------
set.seed(1) 

compareFUZZ <- compare.clusters(data, groups = 10)
BKMres <- compareFUZZ$bkm

BKMplotDF <- data.frame(X = data[,1], Y = data[,2],
                        U = 1 - apply(BKMres$U, 1, max),
                        Z = factor(BKMres$clusters))

FKMres <- compareFUZZ$fkm
FKMplotDF <- data.frame(X = data[,1], Y = data[,2], 
                        U = 1 - apply(FKMres$U, 1, max),
                        Z = factor(FKMres$clus[,1]))

cpal <- scales::hue_pal()(10)

mplots <- function(df, title) {
  ggplot(df, aes(x = X, y = Y, col = Z, size = U)) +
    geom_point(alpha = 0.85) +
    coord_equal() +
    scale_color_manual(values = cpal) +
    scale_size_continuous(limits = c(0,1), range = c(1,5)) +
    labs(x = expression(X[1]),
         y = expression(X[2]),
         col = "Cluster",
         size = "Uncertainty",
         title = title) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, color = "black"),
          legend = "bottom")
}

p1 <- mplots(BKMplotDF, "bootkmeans")
p2 <- mplots(FKMplotDF, "fuzzy cmeans")
(p1 + p2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "center")


## -----------------------------------------------------------------------------
cat("Bootkmeans clustering accuracy: \n")
sum(diag(Thresher::matchLabels(table(BKMres$clusters, hardU))))/nrow(data)
cat("\nFuzzy cmeans clustering accuracy: \n")
sum(diag(Thresher::matchLabels(table(FKMres$clus[,1], hardU))))/nrow(data)

## -----------------------------------------------------------------------------
cat("Bootkmeans FARI: \n")
fari(U, BKMres$U)
cat("\nFuzzy cmeans FARI: \n")
fari(U, FKMres$U)

