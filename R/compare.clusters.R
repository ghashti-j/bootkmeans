compare.clusters <- function(data = NULL,
                             groups = NULL,
                             seed = 13462,
                             nstart = 50,
                             what = "all") {
  if (is.null(data)) stop("argument 'data' is missing with no default")
  if (is.null(groups)) stop("argument 'groups' is missing with no default: specify number of clusters")
  if (!is.null(seed)) set.seed(seed)

  km  <- stats::kmeans(data, groups, nstart = nstart)
  bkm <- boot.kmeans(data, groups, nstart = nstart, returnall = FALSE)

  if (identical(what, "all")) {
    fkm <- fclust::FKM(data, groups, RS = nstart)
    list(km = km, bkm = bkm, fkm = fkm, what = what)
  } else {
    list(km = km, bkm = bkm, what = what)
  }
}
