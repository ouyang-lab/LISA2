# library(Rtsne)
# library(igraph)
# library(plotly)
# library(ggplot2)
# library(pdfCluster)
# library(RANN)
# library(flashClust)
# library(dimRed)
# library(RColorBrewer)
# library(pheatmap)
# library(colorRamps)
# library(RSpectra)
# library("netbiov")
#remotes::install_github("jlmelville/rnndescent")

#' Do log2 and filtering genes and cells
#' @param data matrix, the input data with row (genes) and column (cells)
#' @param dolog logical, whether to do log2, default FAlSE, not do log2
#' @param mingenev number, if the maximum expression value for the gene lower than the threshold,
#'                the gene is filttered.
#' @param cellPor number, if the proportion of non-zero genes of the cell is lower than the threshold,
#'                 the cell is filttered
#' @param cvth number, if the variance of the gene is lower than the threshold,
#'                 the gene is filttered
#' @return
#'    the log2 and filtered data
#' @export
#'
#' @examples
log2Filter <- function(data, dolog=FALSE, mingenev = 1, cellPor = 0.01, cvth = 1){
  if(dolog){
    data <- log2(data+1)
  }
  data <- data[apply(data, 1, max)>=mingenev,]  # Filter low express genes
  cellv <- apply(data,2, function(x){ length(x[x>0])/length(x) }) # Filter low express cells
  data <- data[,cellv>cellPor]
  cv <-  apply(data, 1, sd)/rowMeans(data) # Filter low variance genes
  data <- data[cv > cvth,] #
  data
}


#' Find nearest peak for each valley
#' @param valley
#' @param peaks
#' @param dst vector, the density vector of each point
#' @param nid matrix, the knn matrix
#' @return
#'    the vector of nearest peak for each valley
#' @export
#'
#' @examples
climbMount <- function(valley, peaks, dst, nid){
  vlen = length(valley)
  vmemb = valley
  for(i in 1:vlen){
    va = valley[i]
    neighbors = nid[va,2:ncol(nid)]
    neidst = dst[neighbors]
    nmax = neighbors[which.max(neidst)]
    while(!(nmax %in% peaks)){
      neighbors = nid[nmax,2:ncol(nid)]
      neidst = dst[neighbors]
      nmax = neighbors[which.max(neidst)]
    }
    vmemb[i] = nmax
  }
  vmemb
}

#' Get the minimum spanning tree
#' @param Y The data matrix to build the MST
#' @return
#'    the igraph object
#' @export
#'
#' @examples
getCenMST <- function(Y){
  B<-as.matrix(dist(Y))
  gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
  dp_mst0 = mst(gp)
  return(dp_mst0)
}


#' Calculate the distance form point A to edge OB
#' @param A point
#' @param O a point of edge OB
#' @param B a point of edge OB
#' @return
#'    a vector of distance:
#'    length of edge AP
#'    length of edge OP
#'    length of edge OB
#'    the cosin value of the degree of edge OA and OB
#' @export
#'
#' @examples
#' A=c(2,1) or A=c(1,1) or A=c(1,0)
#' B=c(2,0)
#' O=c(0,0)
#' ds_p2OB(A,O,B)
ds_p2OB <- function(A,O,B){##
  OA = A - O
  OB = B - O
  dOB = sqrt(sum(OB^2))
  dOA = dist(rbind(A,O))[1]
  dAB = dist(rbind(A,B))[1]
  if(dOA>0&&dAB>0){
    cosin = sum(OA*OB) / ( dOB * dOA)
    OP = OB*(sqrt(sum(OA^2)) * cosin /sqrt(sum(OB^2)))
    AP = OP - OA
    dAP = sqrt(sum(AP^2))
    dOP = sqrt(sum(OA^2)) * cosin
  }else if(dOA==0){
    OP = OA
    dAP = 0
    cosin = 0
    dOP = 0
  }else if(dAB==0){
    OP = OB
    dAP = 0
    cosin = 0
    dOP = sqrt(sum(OB^2))
  }
  sign = getLR(O,B,A)
  a=c(dAP, dOP, dOB,cosin, sign)
  names(a) = c('dAP',"dOP",'dOB',"cosin", "sign")
  a
}



mean.of.logs <- function(x, base = 2)
{
  log(mean((base^x) - 1) + 1, base = base)
}


prop.exp <- function (x)
{
  return(length(which(x > 0))/length(x))
}


getMapId <- function(rawMap, vecs){
  allid = sapply(vecs, function(x){which(rawMap==x)})
  x <- as.integer(allid)
  na.omit(x)
}


defaultURDContinuousColors <- function (with.grey = F)
{
  if (with.grey) {
    return(c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0",
             "#307DF0", "#0065FF", "#0078FF", "#008DFF", "#00A1FF",
             "#00B5FF", "#00CAFF", "#00DEFF", "#00F2FF", "#27FFD7",
             "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", "#FFC900",
             "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00",
             "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200",
             "#FF0000"))
  }
  else {
    return(c("#0000FF", "#0013FF", "#0028FF", "#003CFF",
             "#0050FF", "#0065FF", "#0078FF", "#008DFF", "#00A1FF",
             "#00B5FF", "#00CAFF", "#00DEFF", "#00F2FF", "#27FFD7",
             "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", "#FFC900",
             "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00",
             "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200",
             "#FF0000"))
  }
}

getLR <- function(O,B,x){# 0 on the line, and +1 on one side, -1 on the other side.
  sign((B[1] - O[1]) * (x[2] - O[2]) - (B[2] - O[2]) * (x[1] - O[1]))
}



getCenters <- function(dat, cluid){
  allM = as.integer(names(table(cluid)))
  allM = allM[order(allM)]
  print(allM)
  Y = t(sapply(allM,
               function(i,dat) {
                 ind = 1:length(cluid)
                 ind = ind[cluid==i]
                 if(length(ind)>1){
                   rowMeans(dat[,ind])
                 }else if(length(ind)==1){
                   dat[,ind]
                 }
               },t(dat) ))

  rownames(Y) <- paste("Y_", 1:nrow(Y), sep = "")
  return(Y)
}



getLayout2 <- function(X, cluid, root=1){
  Y <- getCenters(X, cluid) #   cluidL
  mst <-  getCenMST(Y)
  #plot(mst, layout=layout_as_tree(mst, root=root))
  ####1. convert undirect to direct based root
  dmst = as.directed(mst, mode = "mutual" )

  ##delete edges direct to root
  nb = neighbors(dmst, v=root, mode =  "all")
  as.vector(nb)
  dnb = bfs(dmst, root=root, neimode="out")
  od = as.vector(dnb$order)
  for(i in 2:length(od)){
    ch = od[i]
    for(j in 1:(i-1)){
      fa = od[j]
      nb = as.vector(neighbors(dmst, v=ch, mode =  "all"))
      if(fa %in% nb){
        dmst = delete_edges(dmst, edges= paste("Y_",ch,"|","Y_",fa,sep="") )
      }
    }
  }
  plot(dmst)
  ########2.get layout
  gp = as_data_frame(dmst)
  ggraph(dmst, 'igraph', algorithm = 'tree') +
    geom_edge_link() +
    ggforce::theme_no_axes()

  igraph.layout <- create_layout(dmst, layout = "dendrogram")
  igraph.layout
}


getCluCenterID <- function(tZ, cluid){
  ## estimate density of t-SNE
  cs = as.integer(names(table(cluid)))
  peaks = 1:length(cs)
  for(i in  1:length(cs)){
    id = which(cluid==i)
    if(length(id)==1){
      peaks[i] = id
    }else{
      cenX = as.matrix(colMeans(tZ[id,]))
      ds = proxy::dist( tZ[id,], t(cenX))
      sid = which.min(ds)
      peaks[i] = id[sid]
    }
  }
  peaks
}

getCluCenters <- function(tZ, cluid){
  ## estimate density of t-SNE
  cs = as.integer(names(table(cluid)))
  mx = matrix(0, nrow = length(cs), ncol = ncol(tZ))
  rownames(mx) = cs
  peaks = 1:length(cs)
  for(i in  1:length(cs)){
    id = which(cluid==i)
    if(length(id)==1){
      mx[i,] = tZ[id,]
    }else{
      mx[i,]  = as.matrix(colMeans(tZ[id,]))
    }
  }
  mx
}

### Build L-ISOMAP with automatically L1 and L2 size to build KNN graph
autoLISOMAP <- function(Z, peaks=NULL, ndim = 3, L1=NULL, L2=NULL,
                        N=2)
{## connect the cells to its nearest landmark points
  ####### Landmakr points
  #peaks = getCluCenterID(Z, cluid)
  if(is.null(peaks)){
    message(Sys.time(), "Please set the peaks!!!")
    return()
  }

  indata = Z[peaks, ]
  if(is.null(L1) && is.null(L2)){
    L1 = length(peaks) -1
    L2 = length(peaks) -1
  }
  if(is.null(L1)){
    L1 = 2
  }
  if(is.null(L2)){
    L2 = 2
  }

  ## geodesic distances
  message(Sys.time(), ": constructing knn graph")
  knng <- makeKNNgraph(x = indata, k = L1, eps = 0)
  n = components(knng$g, mode =  "strong" )
  while(n$no>1){
    print("L1 is too small to build a connected graph!")
    L1 = L1 + 1
    knng <- makeKNNgraph(x = indata, k = L1, eps = 0)
    n = components(knng$g, mode =  "strong" )
    print(L1)
  }
  # plot(knng$g)
  message(Sys.time(), ": calculating geodesic distances")
  geodist <- igraph::distances(knng$g, algorithm = "dijkstra")

  message(Sys.time(), ": Classical Scaling")
  ## TODO: add regularization
  k <- geodist ^ 2
  k <- .Call(stats:::C_DoubleCentre, k)
  k <- - k / 2
  ## TODO: explicit symmetrizing
  ## TODO: return eigenvectors?
  e <- RSpectra::eigs_sym(k, ndim, which = "LA", opts = list(retvec = TRUE))
  e_values <- e$values
  e_vectors <- e$vectors
  neig <- sum(e_values > 0)
  if (neig < ndim) {
    warning("Isomap: eigenvalues < 0, returning less dimensions! Please set a larger L1!!!")
    e_values <- e_values[seq_len(neig)]
    e_vectors <- e_vectors[, seq_len(neig), drop = FALSE]
    return()
  }
  e_vectors <- e_vectors * rep(sqrt(e_values), each = nrow(e_vectors))
  colnames(e_vectors) <- paste("iso", seq_len(neig))

  orgdata = indata
  ################ Apply to left points
  indata <- Z
  nindata <- nrow(indata)
  norg <- nrow(orgdata)

  message(Sys.time(), ": constructing knn graph")

  lknng <- makeKNNgraph(rbind(indata, orgdata), k = L2, eps = 0)
  n = components(lknng$g, mode =  "strong" )
  while(n$no>1){
    print("L2 is too small to build a connected graph!")
    L2 = L2 + 1
    lknng <- makeKNNgraph(rbind(indata, orgdata), k = L2, eps = 0)
    n = components(lknng$g, mode =  "strong" )
    print(L2)
  }

  message(Sys.time(), ": calculating geodesic distances")
  #lgeodist <- igraph::distances(lknng$g,
  #                              seq_len(nindata),
  #                              nindata + seq_len(norg))

  #t1 = proc.time()
  #lgeodist1 <- t(igraph::distances(lknng$g,
  #                                nindata + seq_len(norg),
  #                                seq_len(nindata)))
  #proc.time() -t1
  if(N>1){
    x = nindata + seq_len(norg)
    chunk = split(x,cut(x,quantile(x,(0:N)/N), include.lowest=TRUE, labels=FALSE))
    registerDoParallel(N)
    system.time({
      r <- foreach(i=1:N, .combine=cbind) %dopar% {
        t(igraph::distances(lknng$g,
                            chunk[[i]],
                            seq_len(nindata)))
      }
    })
  }else{
    lgeodist <- t(igraph::distances(lknng$g,
                                    nindata + seq_len(norg),
                                    seq_len(nindata)))
  }


  message(Sys.time(), ": embedding")

  dammu <- sweep(lgeodist ^ 2, 2, colMeans(geodist ^ 2), "-")
  Lsharp <- sweep(e_vectors, 2, e_values, "/")
  out <- -0.5 * (dammu %*% Lsharp)
  rownames(out) = rownames(Z)
  return(out)
}



makeKNNgraph <- function (x, k, eps = 0, diag = FALSE){
  ## requireNamespace("RANN")
  ## requireNamespace("igraph")

  ## consts
  INF_VAL <- 1.340781e+15
  NA_IDX  <- 0
  BDKD_LIM <- 1000000                 #todo: figure out a good value here

  ## select parameters
  M <- nrow(x)
  treetype <- "kd"                # if (M < BDKD_LIM) "kd" else "bd"
  # see:
  # https://github.com/jefferis/RANN/issues/19
  searchtype <- if (eps == 0) "standard" else "priority"

  ## RANN::nn2 returns the points in data with respect to query
  ## e.g. the rows in the output are the points in query and the
  ## columns the points in data.
  nn2res <- RANN::nn2(data = x, query = x, k = k + 1, treetype = treetype,
                      searchtype = searchtype, eps = eps)

  ## create graph: the first ny nodes will be y, the last nx nodes
  ## will be x, if x != y
  ## it is not really pretty to create a
  ## directed graph first and then make it undirected.
  g <- igraph::make_empty_graph(M, directed = TRUE)
  g[from = if (diag) rep(seq_len(M), times = k + 1)
    else      rep(seq_len(M), times = k),
    to   = if (diag) as.vector(nn2res$nn.idx)
    else      as.vector(nn2res$nn.idx[, -1]),
    attr = "weight"] <-
    if (diag)  as.vector(nn2res$nn.dists)
  else as.vector(nn2res$nn.dists[, -1])

  g = igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "first")

  return(list(g=g, nn2res=nn2res ))
}

optDim <- function(sdev, pvalue){
  len = length(sdev)
  sa = sdev[2:len] - sdev[1:(len-1)]
  sm = summary(sa)

  id = which(sa < sm[5])

  df = data.frame(y = sa, x = 1:length(sa))
  rf = lm(y~ x, df[1:nrow(df),])
  mr = mean(rf$residuals)


  x=1:nrow(df)
  y=predict(rf,newdata=list(x=x))

  sv = (sa -y )^2
  sv = sv[1:100]
  mv = summary(sv)
  #fit <- fitdistr(sv, "normal")
  #class(fit)
  th = qnorm(1-pvalue, mean = mv[3], sd = mv[3], lower.tail = TRUE, log.p = FALSE)
  pcdim = which(sv > th)
  length(pcdim)

  if(1==0){
    plot(sa)
    abline(h=sm[5], col="blue")
    plot(rf,which = 1)
    plot(x, sa)
    points(x,y, col = "red")
    plot(sv[10:30])
  }
  pcdim
}


isLeafNode <- function(hid, leaves){
  for(i in 1:length(leaves)){
    if(hid %in% leaves[[i]]) return(i)
  }
  return(0)
}

