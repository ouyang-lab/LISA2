
#' laodExample: load example dataset
#'
#' @param
#' @return
#'   data matrix
#' @export
#'
#' @examples
loadExample <- function(){
  load('SLS3279.RData')
  rownames(SLS3279$data) = paste("gene", 1:nrow(SLS3279$data), sep="")
  colnames(SLS3279$data)= paste("cell", 1:ncol(SLS3279$data), sep="")
  list(data=SLS3279$data, t=SLS3279$t)
}

#' reduceByPT: Do PCA and t-SNE for given data matrix
#'
#' @param data
#'        a matrix with genes * cells
#' @param threshold
#'     numeric, (default 3), select the PCs based on proportion of variance of PCs
#' @param nDim
#'     numeric, (default NULL), Set the top principle components(PCs) number. If it is not NULL, the threshold will not be used.
#' @param center
#'    logical, (default FALSE), it is used in PCA to do center for data
#' @param scale
#'    logical, (default FALSE), it is used in PCA to do scale for data
#' @return list()
#'    a list with three elements : pca_out, tsne_out, and Z
#' @export
#'
#' @examples  load(data);
#'            reducedData <- reduceByPT(data)
#'
#'
reduceByPT <- function(data, threshold=0.001, nDim = NULL, center=FALSE, scale=FALSE){
  ## run PCA
  pca_out <-  prcomp(t(data), center = center,scale. = scale)
  eigs <- pca_out$sdev^2
  psum <- rbind(SD = sqrt(eigs), Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))
  Proportion = eigs/sum(eigs)
  ta=table(Proportion>threshold)
  if(length(ta)>1){
    pcdim =  1:ta[2]
  }else{
    pcdim =  1:ta[1]
  }
  if(!is.null(nDim) ){
    if(length(nDim)==1){
      pcdim = 1:nDim
    }else{
      pcdim = nDim
    }
  }
  Z = pca_out$x[,pcdim]
  print(dim(Z))
  ## run tSNE
  set.seed(1)
  tsne_out <- Rtsne(Z, dims = 2, pca=FALSE, initial_dims = ncol(Z), check_duplicates = FALSE)
  tZ =  tsne_out$Y
  print(dim(tZ))
  return(list(pca_out=pca_out, tsne_out=tsne_out, Z=Z))
}


#' pcaL: Do PCA and plot the PCs importance
#'
#' @param data
#'        a matrix with genes * cells
#' @param threshold
#'     numeric, (default 3), select the PCs based on proportion of variance of PCs
#' @param nDim
#'     numeric, (default 1:20), Set the top principle components(PCs) number. If it is not NULL, the threshold will not be used.
#' @param center
#'    logical, (default TRUE), it is used in PCA to do center for data
#' @param scale
#'    logical, (default FALSE), it is used in PCA to do scale for data
#' @return list()
#'    a list with three elements : pca_out and Z
#' @export
#'
#' @examples  load(data);
#'            reducedData <- pcaL(data)
#'
#'
pcaL <- function(data, threshold=0.001, nDim = NULL, center=TRUE, scale=FALSE){
  message(Sys.time(), ": Start for PCA")
  ## run PCA
  pca_out <-  prcomp(t(data), center = center,scale. = scale)
  eigs <- pca_out$sdev^2
  psum <- rbind(SD = sqrt(eigs), Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))
  Proportion = eigs/sum(eigs)

  ta=table(Proportion>threshold)
  if(length(ta)>1){
    pcdim =  1:ta[2]
  }else{
    pcdim =  1:ta[1]
  }
  plot(Proportion[pcdim], main =
         paste("The variance importance of top PCs (>threhold:", threshold, ")" ))

  if(!is.null(nDim) ){
    if(length(nDim)==1){
      pcdim = 1:nDim
    }else{
      pcdim = nDim
    }
  }
  Z = pca_out$x[,pcdim]
  print(dim(Z))
  message(Sys.time(), ": End of PCA")
  return(list(pca_out=pca_out, Z=Z))
}

#' pcaL1: Do PCA and select the number of PCs
#'
#' @param data
#'        a matrix with genes * cells
#' @param pvalue
#'     numeric, (default 0.05), select the PCs based on significant value for normal distribution
#' @param pcaSel
#'     numeric, (default 1 or 2), if pcaSel=1, matrix size < 10^8  109.78s;
#'     if sel=2, For matrix size < 10^10
#' @param sel
#'     numeric, (default 1 or 2), if sel=1, use normal distibution to select top PCs;
#'     if sel=2, do linear model to select top PCs from the top 100 PCs; it is faster than sel=1.
#'     And the number of PCs selected are often larger than sel=1.
#' @param center
#'    logical, (default TRUE), it is used in PCA to do center for data
#' @param scale
#'    logical, (default FALSE), it is used in PCA to do scale for data
#' @return list()
#'    a list with three elements : pca_out and the number of PCs selected
#' @export
#'
#' @examples  load(data);
#'            reducedData <- pcaL1(data) # For matrix size < 10^8
#'            reducedData <- pcaL1(data,sel=2) # For matrix size < 10^10
#'
#'
pcaL1 <- function(data, pvalue=0.05, pcaSel = 1, sel=1, center=TRUE, scale=FALSE, ith=1e-2){
  message(Sys.time(), ": Start for PCA")
  ## run PCA
  t1 = proc.time()
  if(pcaSel==1){#For matrix size< 10^8  109.78s
    set.seed(1)
    pca_out <-  prcomp(t(data), center = center,scale. = scale)
  }else if(pcaSel==2){#500PC:598s, 100PC:21s  #For matrix size < 10^10
    set.seed(1)
    pca_out <- prcomp_irlba(t(data), n = 100, center = center, scale. = scale)
    #summary(p1)
  }
  print(proc.time() - t1)
  message(Sys.time(), ": PCA finished")

  eigs <- pca_out$sdev^2 #   standard deviations
  psum <- rbind(SD = sqrt(eigs), Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))
  Proportion = eigs/sum(eigs)

  if(sel==1){
    fit <- fitdistr(Proportion, "normal")
    class(fit)
    para <- fit$estimate
    m2 = max(Proportion)
    m1 = min(Proportion)
    #x <- seq(m1, m2, by = (m2-m1)/100)
    #y <- dnorm(x, mean = para[1], sd = para[2])
    th = qnorm(1-pvalue, mean = para[1], sd = para[2], lower.tail = TRUE, log.p = FALSE)
    #pnorm(th, mean = para[1], sd = para[2], lower.tail = TRUE, log.p = FALSE)
    pcdim = which(Proportion > th)
  }else if(sel==2){
    pcdim = optDim(sdev = pca_out$sdev, pvalue)
  }

  threshold = Proportion[pcdim[length(pcdim)]]

  if(threshold > ith) threshold = ith
  plot(Proportion[pcdim], main =
         paste("The variance importance of top PCs" ))

  message(Sys.time(), ": End of PCA")
  return(list(pca_out=pca_out, n=length(pcdim)))
}


getPCAN <- function(pca_out, pvalue=0.05, sel=1, ith=1e-2){
  eigs <- pca_out$sdev^2 #   standard deviations
  psum <- rbind(SD = sqrt(eigs), Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))
  Proportion = eigs/sum(eigs)

  if(sel==1){
    fit <- fitdistr(Proportion, "normal")
    class(fit)
    para <- fit$estimate
    m2 = max(Proportion)
    m1 = min(Proportion)
    #x <- seq(m1, m2, by = (m2-m1)/100)
    #y <- dnorm(x, mean = para[1], sd = para[2])
    th = qnorm(1-pvalue, mean = para[1], sd = para[2], lower.tail = TRUE, log.p = FALSE)
    #pnorm(th, mean = para[1], sd = para[2], lower.tail = TRUE, log.p = FALSE)
    pcdim = which(Proportion > th)
  }else if(sel==2){
    pcdim = optDim(pca_out$sdev, pvalue)
  }

  threshold = Proportion[pcdim[length(pcdim)]]

  if(threshold > ith) threshold = ith
  plot(Proportion[pcdim], main =
         paste("The variance importance of top PCs" ))

  return(length(pcdim))
}



#' tsneL: Do tSNE
#'
#' @param Z
#'        PCA
#' @param tdim
#'     numeric, (default 2), how many dimensions for tSNE
#' @param center
#'    logical, (default FALSE), it is used in PCA to do center for data
#' @param scale
#'    logical, (default FALSE), it is used in PCA to do scale for data
#' @return list()
#'    tsne_out
#' @export
#'
#' @examples
#'
#'
#'
tsneL <- function(Z,  tdim=2,  center=FALSE, scale=FALSE, seeds=1){
  ## run tSNE
  message(Sys.time(), " Start of tSNE")
  set.seed(seeds)
  tsne_out <- Rtsne(Z, dims = tdim, pca=FALSE, initial_dims = ncol(Z), check_duplicates = FALSE)
  tZ =  tsne_out$Y
  print(dim(tZ))
  message(Sys.time(), " End of tSNE")
  return(tsne_out)

}

#' umpL: Do UMAP
#'
#' @param Z
#'        PCA
#' @param tdim
#'     numeric, (default 2), how many dimensions for tSNE
#' @param center
#'    logical, (default FALSE), it is used in PCA to do center for data
#' @param scale
#'    logical, (default FALSE), it is used in PCA to do scale for data
#' @param seeds
#'    numeric, (default 1), it is used to reproduce same UMAP
#' @return list()
#'    tsne_out
#' @export
#'
#' @examples
#'
#'
#'
umpL <- function(Z, seeds=1){
  ## run tSNE
  message(Sys.time(), " Start of UMAP")
  set.seed(seeds)
  emb2 <-  umap::umap(Z)
  uY = emb2$layout
  message(Sys.time(), " End of UMAP")
  return(uY)
}


searchPeaks <- function(gg, cluid){
  # select cell which has most neighbors and minimum distance to neighbors
  deg = degree(gg,v=1:length(cluid))
  #plot(degree_distribution(gg))
  nc = as.integer(names(table(cluid)))

  peaks = NULL

  for(i in 1:length(nc)){
    c1 = which(cluid==nc[i])
    p1 = deg[c1]
    mxid = c1[which(p1==max(p1))] # 1. select node which has most neighbors
    if(length(mxid)>1){
      mds = NULL
      for (j in 1:length(mxid)) {# 2. Select node which has smallest mean distance to neighbors
        nn = as.vector(neighbors(gg, v=mxid[j]))
        mds = c(mds, mean(distances(gg, v=mxid[j], to=nn)))
      }
      kid = which.min(mds)
      peaks = c(peaks, mxid[kid])
    }else{
      peaks = c(peaks, mxid)
    }
  }
  return(peaks)
}

#' calDensity: Estimate the density of t-SNE
#'
#' @param tZ matrix, the t-SNE results
#' @param fixed logical, (default TRUE), Either 'fixed' or 'adaptive',
#'  corresponding to a kernel estimator with fixed or adaptive bandwidths respectively. See kepdf().
#' @param ksize integer, (default 50), it is used to build KNN
#' @return
#'  a list with three elements : pdf,peaks,valley
#'  pdf, the results of kepdf, users can plot the 3D density plot using pdf. See kepdf()
#'  peaks, the peaks nodes
#'  valley, all valley nodes
#' @export
#'
#' @examples
calDensity <- function(tZ, fixed= TRUE, merge=FALSE, ksize=50){
  ## estimate density of t-SNE
  set.seed(1)
  if(fixed){
    pdf <- kepdf(tZ)
  }else{
    pdf <- kepdf(tZ, bwtype = "adaptive")
  }
  ## Get peaks and valley by k-NN graph
  dst = pdf@estimate
  zo = rep(0, length(pdf@par$h))
  h = proxy::dist(t(zo), t(pdf@par$h))[1]

  kdt <- nn2(tZ, k= ksize, searchtype="standard")
  nid = kdt$nn.idx
  peaks = NULL
  valley = NULL
  for(i in 1:nrow(nid)){
    a = dst[nid[i,]]
    if(1 == which.max(a)){
      peaks=c(peaks, i)
    }
    if(1==which.min(a)){
      valley = c(valley, i)
    }
  }
  cluN = length(peaks)

  if(merge){
    mergePeaks = list()
    for(i in 1:length(peaks)){
      k = peaks[i]
      id = which(nid[k,] %in% peaks)
      np = nid[k, id]
      if(length(id)>1){
        mergePeaks[[i]]= np
        pdens = dst[np]
        mp = np[order(np)]
        rid = which(peaks %in% mp[-1])
        peaks = peaks[-rid]
      }else{
        mergePeaks[[i]]=NULL
      }
    }
  }

  print(length(peaks))
  list(pdf=pdf, peaks=peaks, valley=valley)
}



#' denknnCluT: Do hierarchical clustering based on density of t-SNE
#'
#' @param desRes list, the results of calDensity()
#' @param tZ matrix, the t-SNE results
#' @return
#'  memb, The clustering id
#' @export
#'
#' @examples
denknnCluT <- function(denRes, tZ, cmethod1 = "ward",cmethod2 = "ward" ){

  peaks <- denRes$peaks
  valley <- denRes$valley
  pdf <- denRes$pdf
  if(is.null(denRes$cluN)){
    cluN = length(peaks)
  }else{
    cluN = denRes$cluN
  }

  print(cluN)
  ### 1-level cluster
  hc <- flashClust::hclust(dist(tZ), method = cmethod1)
  hc$height <- hc$height^2
  memb <- cutree(hc, k = cluN)

  ### 2-level cluster
  if(cluN>=2){
    tp = table(memb[peaks])
    ncenter = length(table(memb))
    while (length(tp)<length(peaks)) { ## Divid the cluster until each cluster contains only one peak
      curM = 1
      rmb = memb
      aid = 1:nrow(tZ)
      for(i in 1:ncenter){
        lid = aid[memb==i]
        pf = peaks %in% lid
        pid = (1:length(peaks))[pf]
        len = length(pid)
        if(len>1){## do clustering
          hc <- flashClust::hclust(dist(tZ[lid,]), method = cmethod2)
          hc$height <- hc$height^2
          sm <- cutree(hc, k = len) + curM -1
          rmb[lid] = sm
          curM = curM + len
          # plot(tZ[lid,], col=sm)
        }else if(len>=0){ ### if one clusters contains valley from different cluster, divide them
          vf = valley %in% lid
          vid = (1:length(valley))[vf]
          tb = table(memb[vid])
          if(length(tb)>1){
            tid = 1:length(tb)
            vlen = length(tid[tb>1])
            if(vlen>1){
              hc <- flashClust::hclust(dist(tZ[lid,]), method = cmethod2)
              hc$height <- hc$height^2
              sm <- cutree(hc, k = vlen) + curM -1
              rmb[lid] = sm
              curM = curM + vlen
            }else{
              rmb[lid] = curM
              curM = curM + 1
            }
          }else{
            rmb[lid] = curM
            curM = curM + 1
          }
        }else{
          rmb[lid] = curM
          curM = curM + 1
        }
      }

      ncenter = curM-1
      print(table(memb))
      print(table(rmb))
      memb = rmb
      tp = table(memb[peaks])
    }

    ### Assign cluster without peaks to cluster with one peak
    pm = memb[peaks]
    allM = as.integer(names(table(memb)))
    ef = which(allM %in% unique(pm))
    noPeakMem = allM[-ef]
    if(length(noPeakMem)>0){
      df = length(noPeakMem)
      af = df
      while(df>0){
        pcluDist =matrix(0, nrow = length(noPeakMem), ncol = length(peaks))
        for(i in 1:length(noPeakMem)){
          ef = memb==noPeakMem[i]
          mz = tZ[ef,]
          if(is.null(dim(mz))){
            mz = t(as.matrix(mz))
          }
          a= proxy::dist(tZ[peaks,], mz)
          pcluDist[i,] = apply(a, 1, min)
        }

        minDist = apply(pcluDist, 1, min)
        mind = which.min(minDist)
        ## assign mind to its nearest peaks
        minP = pm[which.min(pcluDist[mind,])]
        svmNum = svmCluster(tZ,i=noPeakMem[mind],j=minP,memb)
        if(svmNum>0){
          af = memb==noPeakMem[mind]
          memb[af] = minP
        }
        ###
        #allM = as.integer(names(table(memb)))
        #ef = which(allM %in% unique(pm))
        #noPeakMem = allM[-ef]
        noPeakMem = noPeakMem[-mind]
        df = length(noPeakMem)
      }
    }
  }
  print(table(memb))
  ta = as.integer(names(table(memb)))

  rmb = memb
  for(i in 1:length(ta)){
    ef = memb==ta[i]
    memb[ef] = i
  }
  print(max(memb))
  memb
}


#' setLeignPath: set the path of python package which installed the leidenalg package
#'
#' @param pyPath string, python command path
#' @param condaPath string, conda install path
#' @param libPath string, python library install path
#' @return NULL
#'
#' @export
#'
#' @examples
setLeignPath <- function(pyPath="/Users/cheny48/miniconda3/bin/python",
                         condaPath="/Users/cheny48/miniconda3",
                         libPath='/Users/cheny48/miniconda3/lib/python3.8/site-packages'){
  library(reticulate)
  reticulate::use_condaenv(condaPath, conda = "auto", required = T)
  reticulate::use_python(pyPath, required = T)
  reticulate::import_from_path("leidenalg", path =libPath,
                               convert = TRUE, delay_load = FALSE)
  reticulate::py_config()
  reticulate::py_module_available("leidenalg")
}

#' getLeidenCluster: apply leiden community clustering
#'
#' @param graph0 graph object, input
#' @return partition, vector, clustering results
#'
#' @export
#'
#' @examples
getLeidenCluster <- function(graph0,seed=1){
  ##user need to set local python path and leidenalg install path
  library("leiden")
  partition <- leiden(graph0, seed=seed)
  partition
}


#' cmcluster: Build clustering based KNN graph and community clustering
#'
#' @param Z matrix, the PCA results
#' @param ksize integer, (default 2), neighbor size to build knn graph
#' @param cmmethod string, (default louvain), user can also try fast_greedy method (cluster_fast_greedy),
#'                 walktrap ( cluster_walktrap() ) and infomap ( cluster_infomap())
#' @param stepN integer, used in walktrap method cluster_walktrap()
#' @param clutrialN integer, (default 10), for infomap method cluster_infomap()
#' @return list
#'         knng, knn graph
#'         cluid, clustering results
#' @export
#'
#' @examples
cmcluster <- function(Z, ksize = 2, cmmethod="louvain", stepN=4, clutrialN = 10,seed=1){
  ## geodesic distances
  message(Sys.time(), ": constructing knn graph")
  knng <- makeKNNgraph(x = Z, k = ksize, eps = 0)
  n = components(knng$g, mode =  "strong" )

  # plot(knng$g)
  message(Sys.time(), ": calculating geodesic distances")
  geodist <- igraph::distances(knng$g, algorithm = "dijkstra")

  message(Sys.time(), ": do clustering")

  if(cmmethod=="louvain"){
    set.seed(seed)
    clures = cluster_louvain(knng$g)
  }else if(cmmethod=="fast_greedy"){
    set.seed(seed)
    clures = cluster_fast_greedy(knng$g, merges = TRUE, modularity = TRUE, membership = TRUE)
  }else if(cmmethod=="walktrap"){
    set.seed(seed)
    clures = cluster_walktrap(knng$g,  steps = stepN, merges = TRUE,
                              modularity = TRUE, membership = TRUE)
  }else if(cmmethod=="infomap"){
    set.seed(seed)
    clures = cluster_infomap(knng$g, e.weights = NULL, v.weights = NULL,
                             nb.trials = clutrialN, modularity = TRUE)
  }else if(cmmethod=="leiden"){
    cluid = getLeidenCluster(knng$g,seed)
  }
  if(cmmethod %in% c("louvain", "fast_greedy", "walktrap", "infomap")){
    cm = communities(clures)
    cluid = 1:nrow(Z)
    for(i in 1:length(cm)){ cluid[cm[[i]]] = i }
  }
  message(Sys.time(), ": Finished clustering")
  return(list(knng=knng, cluid=cluid))
}


plotCommNet <- function(Z, ksize, cluid, colname= rownames( brewer.pal.info)[12],
                        gtype = c("star", "Kamada-Kawai", "abstractStar", "abstractFruchtermanReingold"),
                        seedN=1){
  data("color_list")
  clus = unique(cluid)
  clus = clus[order(clus)]
  knng <- makeKNNgraph(x = Z, k = ksize, eps = 0)
  mc1 = list()
  for(i in 1:length(clus)){
    mc1[[i]] = which(cluid==clus[i])
  }
  names(mc1) = clus
  cs = colorRampPalette(brewer.pal(brewer.pal.info[colname,1],colname))(length(clus))
  nodecs = sapply(cluid, function(x) { cs[x ] })

  set.seed(seedN)
  if(gtype=="star"){
    id <- plot.modules(knng$g, layout.function = c(layout.fruchterman.reingold),mod.list = mc1,
                       layout.overall = layout.star,
                       sf=40, tkplot=FALSE, modules.color=cs, mod.lab=TRUE, bg="#ffffff")

  }else if(gtype=="Kamada-Kawai"){
    id <- plot.modules(knng$g, mod.list = mc1,
                       layout.function = layout.graphopt,layout.kamada.kawai,
                       sf=60, tkplot=FALSE, modules.color=cs, mod.lab=TRUE, bg="#ffffff")

  }else if(gtype=="abstractStar"){
    id <- plot.abstract.nodes(knng$g, nodes.color ="grey",layout.function=layout.star,
                              edge.colors= sample(color.list$bright), tkplot =FALSE,
                              v.sf=-40,
                              mod.list = mc1, modules.color=cs, mod.lab=TRUE, bg="#ffffff")

  }else if(gtype=="abstractFruchtermanReingold"){
    id<- plot.abstract.nodes(knng$g, layout.function=layout.fruchterman.reingold,
                             edge.colors= sample(color.list$bright),nodes.color="grey",
                             v.sf=-40,
                             mod.list = mc1, mod.lab=TRUE,bg="#ffffff")
  }
}

#' doLmISOMAP: Do landmark ISOMAP based on peaks and valleys
#'
#' @param Z matrix, the PCA results
#' @param ndim integer, (default 2), set the dimensions of ISOMAP
#' @param peaks vector, the selected density peaks nodes
#' @param valley vector, the selected density valley nodes
#' @param landMark logical, (default TRUE), whether use landmark, using landmark can save a lots of time
#' @return
#'    X, the ISOMAP data
#' @export
#'
#' @examples
doLmISOMAP <- function(Z, ndim=2, peaks, valley, landMark=TRUE){
  lid = c(peaks,valley)
  smpL = length(lid)

  if(landMark){
    dat = dimRedData(Z) #
    emb2 <- embed(dat[lid,], "Isomap",knn= smpL-1, ndim=ndim) # smpL-1
    emb3 <- emb2@apply(dat)
    X =  emb3@data
  }else{
    dat = dimRedData(Z) #
    emb2 <- embed(dat, "Isomap",knn= 50, ndim=ndim) # smpL-1
    emb3 <- emb2@apply(dat)
    X =  emb3@data
  }
  return(X)
}


#' autoLISOMAP: Build L-ISOMAP with automatically L1 and L2 size to build KNN graph
#'
#' @param Z
#'    matrix, the PCA results
#' @param peaks
#'    peaks vector, the selected density peaks nodes
#' @param ndim
#'    integer, (default 2), set the dimensions of ISOMAP
#' @param L1
#'     integer, (default NULL),the kNN size for peaks
#' @param L2
#'     integer, (default NULL),the kNN size for all cells
#' @param N
#'     integer, (default 1), if N > 1, do parallel
#' @return
#'    X, the ISOMAP matrix
#' @export
#'
#' @examples
autoLISOMAP <- function(Z, peaks=NULL, ndim = 3, L1=NULL, L2=NULL, N=1){
  ## connect the cells to its nearest landmark points
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
    L1 = length(peaks) -1
  }
  if(is.null(L2)){
    L2 = length(peaks) -1
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
    lgeodist <- foreach(i=1:N, .combine=cbind) %dopar% {
      t(igraph::distances(lknng$g,
                          chunk[[i]],
                          seq_len(nindata)))
    }

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
  message(Sys.time(), ": Done")
  return(out)
}


#' doLmISOMAPTemp: Do landmark ISOMAP using the cluster centers as peaks
#'
#' @param Z matrix, the PCA results
#' @param ndim integer, (default 2), set the dimensions of ISOMAP
#' @param peaks vector, the selected density peaks nodes
#' @param valley vector, the selected density valley nodes
#' @param landMark logical, (default TRUE), whether use landmark, using landmark can save a lots of time
#' @return
#'    X, the ISOMAP data
#' @export
#'
#' @examples
doLmISOMAPTemp <- function(Z, ndim=2, peaks, Lsize=50, landMark=TRUE, clustering=FALSE){
  lid = peaks
  smpL = length(lid)

  if(landMark){# do Landmark ISOMAP
    if(clustering==FALSE){
      dat = dimRedData(Z) #
      emb2 <- embed(dat[lid,], "Isomap",knn= smpL-1, ndim=ndim) # smpL-1
      emb3 <- emb2@apply(dat)
      X =  emb3@data
    }else{
      hc <- flashClust::hclust(dist(Z), method = "ward")
      hc$height <- hc$height^2
      clusters <- cutree(hc, k = Lsize)
      tb = as.numeric(names(table(clusters)))
      cid <-  sapply(tb, function(x) {
                                ind = 1:length(clusters)
                                ind = ind[clusters==x]
                                if(length(ind)>1){
                                  cen = colMeans(Z[ind,])
                                  ds = proxy::dist(t(as.matrix(cen)), Z[ind,])
                                  ind[which.min(ds)]
                                }else if(length(ind)==1){
                                  ind
                                }
                              })

      dat = dimRedData(Z) #
      emb2 <- embed(dat[cid,], "Isomap",knn= Lsize-1, ndim=ndim) # smpL-1
      emb3 <- emb2@apply(dat)
      X =  emb3@data
    }
  }else{# do ISOMAP
    dat = dimRedData(Z) #
    emb2 <- embed(dat, "Isomap",knn= 50, ndim=ndim) # smpL-1
    emb3 <- emb2@apply(dat)
    X =  emb3@data
  }
  return(X)
}



#' calSimilarity: calculate the similarity among the clusters in PCA or ISOMAP space.
#'
#' @param X matrix, the ISOMAP matrix
#' @param Z matrix, the PCA matrix
#' @param cluid vector, the cluster id
#' @param type integer, (default 1), 1: X (ISOMAP) is used for silhouette index  ,
#'        2: Z is used for silhouette index
#'        3: Z is used for Pearson correlation instead of silhouette index
#' @return list()
#'    list, rootDist: the similarity of clusters, cenDist: the distance of clusters
#' @export
#'
#' @examples
calSimilarity <- function(X,Z, cluid, type=1){

  D = Z
  cluN = max(cluid)
  rootDist = matrix(-1, nrow = cluN, ncol = cluN)
  cenDist = matrix(-1, nrow = cluN, ncol = cluN)

  t1 = proc.time()
  for(i in 2:cluN){
    for(j in 1:(i-1)){
      sid = which(cluid %in% c(i,j) )
      if(type==1){
        si2 <- silhouette(cluid[sid], dist(X[sid,], "euclidean"))
        rootDist[i,j] = mean(si2[,3])
        rootDist[j,i] =  rootDist[i,j]
      }else if(type==2){
        si2 <- silhouette(cluid[sid], dist(Z[sid,], "euclidean"))
        rootDist[i,j] = mean(si2[,3])
        rootDist[j,i] =  rootDist[i,j]
      }else if(type==3){
        cors = 1 - cor(t(Z[sid,]))
        rootDist[i,j] = mean(cors)
        rootDist[j,i] = rootDist[i,j]
      }

      id1 = which(cluid == i)
      id2 = which(cluid == j)
      if(length(id1)==1 & length(id2)>1){
        ds = proxy::dist(t(as.matrix(D[id1,])), D[id2,])
      }else if(length(id1)>1 & length(id2)==1){
        ds = proxy::dist(D[id1,], t(as.matrix(D[id2,])))
      }else if(length(id1)==1 &length(id2)==1){
        ds = proxy::dist( t(as.matrix(D[id1,])) , t(as.matrix(D[id2,])))
      }else{
        ds = proxy::dist(D[id1,], D[id2,])
      }
      cenDist[i,j] = mean(ds)
      cenDist[j,i] = mean(ds)
    }
  }


  rootDist = rootDist +1
  diag(rootDist) = 0

  print(proc.time() - t1)

  return(list(rootDist=rootDist, cenDist=cenDist))
}

#' kNNDist: create KNN graph
#'
#' @param rootDist
#'    similariy matrix of clusters
#' @param ks
#'    integer, default is 2
#' @return graph
#'    g1
#' @export
#'
#' @examples
#'
kNNDist <- function(rootDist, ks=2){
  g1 <- igraph::make_empty_graph(nrow(rootDist), directed = FALSE)
  egList = NULL
  ek = 1
  if(ks > nrow(rootDist)){
    ks = nrow(rootDist)
  }
  for(i in 1:nrow(rootDist)){
    rs = order(rootDist[i,])
    for(j in 2:ks){
      g1 = g1 + edge(i, rs[j])
      E(g1)$weight[ek] = rootDist[i,rs[j]]
      ek = ek + 1
    }
  }
  g1 = simplify(g1)
  #plot(g1,vertex.label.cex = 1.5,vertex.size=0.1)
  g1
}

#' getMSTLeafRootISO: create spanning tree with specified root and leaves
#' @param X
#'     ISOMAP matrix
#' @param rootDist
#'     similariy matrix of clusters
#' @param cenDist
#'     distance matrix of clusters
#' @param cluid
#'     vector, the cluster id
#' @param L2
#'    integer, default is 2, knn size
#' @param root
#'    integer, default is 1, the root cluster
#' @param leaves
#'    list, specify the leave groups
#' @return list
#'    list(g5=g5,rootDist=rootDist,w=sum(E(g5)$weight))
#' @export
#'
#' @examples
#'
getMSTLeafRootISO <- function(rootDist, cenDist, cluid, L2=2, root=1,
                              leaves = list(a=c(9,12,13), b=17,c=c(2,6),d=3,e=c(4,5),f=14))
{
  allleaf = NULL
  for(i in 1:length(leaves)){
    allleaf = c(allleaf, leaves[[i]])
  }
  nonleaf = (1:max(cluid))[-allleaf]
  ########1.Calcualte the distance matrxi among each clusters
  ###1.Compute the minimum distance from the leaf cluster to the root cluster
  M = length(unique(cluid))

  if(length(nonleaf)==1){
    ## convert into all nodes graph
    g1 <- igraph::make_empty_graph(nrow(rootDist), directed = FALSE)
  }else{

    ########2.Build kNN graph using the distance matrix of nonleaf nodes
    g0 = kNNDist(rootDist[nonleaf,nonleaf], ks=L2)
    dg = degree(g0, v=1:length(nonleaf))
    dg
    n = components(g0, mode =  "strong" )
    mark=TRUE
    nonleafroot = nonleaf[-which(nonleaf==root)]
    while(n$no>1){
      print("L2 is too small to build a connected graph!")
      L2 = L2 + 1
      g0 <- kNNDist(rootDist[nonleaf,nonleaf], ks=L2)
      n = components(g0, mode =  "strong" )
      print(L2)
    }
    #
    edList= cbind(get.edgelist(g0), E(g0)$weight)
    #plot(g0,layout=layout.reingold.tilford(g0, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

    ## convert into all nodes graph
    g1 <- igraph::make_empty_graph(nrow(rootDist), directed = FALSE)
    for(i in 1:nrow(edList)){
      v1 = nonleaf[edList[i,1]]
      v2 = nonleaf[edList[i,2]]
      g1 = g1 + edge(v1, v2)
      E(g1)$weight[i] = edList[i,3]
    }
  }

  #plot(g1,layout=layout.reingold.tilford(g1, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

  ########3. Add leaf node into the nearest nonleaf node

  len = length(leaves)
  ## Keep the degree of leaf node group to its parent is one
  for(i in 1:len){
    lv = leaves[[i]]
    ##### add edge among leave using MST
    if(length(lv)>1){
      B <- as.matrix(rootDist[lv,lv] )
      gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
      dp_mst = mst(gp)
      se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
      k = length(E(g1)$weight) + 1
      for(z in 1:nrow(se)){
        g1 = g1 + edge(lv[se[z,1]] ,lv[se[z,2]])
        E(g1)$weight[k] = se[z,3]
        k = k + 1
      }
    }
    ## Find nearest nonleaf node of the leaf group
    ls = lv
    ld = lv
    for(j in 1:length(lv)){
      ds = rootDist[lv[j], nonleaf]
      ls[j] = nonleaf[which.min(ds)]
      ld[j] = min(ds)
    }
    id = which.min(ld)
    k = length(E(g1)$weight) + 1
    g1 = g1 + edge(lv[id], ls[id])
    E(g1)$weight[k] = ld[id]
  }
  g1 = simplify(g1)
  edList= cbind(get.edgelist(g1), E(g1)$weight)
  print(nrow(edList))
  #plot(g1,layout=layout.reingold.tilford(g1, circular=F))
  # pieTreeRoot(mst=g1, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)


  ###### 2. Add known path into
  g4 = igraph::make_empty_graph(M, directed = FALSE)
  paths = list() #get.all.shortest.paths(g1,from = leaves, to=root, mode = "all")
  c=1
  for(i in 1:length(leaves)){
    lv = leaves[[i]]
    ds = distances(g1, v=lv, to= root)
    id = which.min(ds[,1])
    p1 = get.shortest.paths(g1, from=root, to=lv[id])
    if(length(p1$vpath)==0){
      message(Sys.time(), "Error: Please set larger L2 !!!")
      return()
    }
    ev = as.vector(p1$vpath[[1]])

    ##### add edge among leaves using MST
    if(length(ev)>1){
      h = length(ev)-1
      for(j in 1:h){
        ei <- get.edge.ids(g4, c(ev[j],ev[j+1]))
        if(ei==0){
          g4 = g4 + edge(ev[j],ev[j+1])
          E(g4)$weight[c] = distances(g1, v=ev[j], to= ev[j+1])[1,1]
          c = c + 1
        }
      }

      if(length(lv) > 1){
        B <- as.matrix(rootDist[lv,lv] )
        gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
        dp_mst = mst(gp)
        se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
        for(z in 1:nrow(se)){
          g4 = g4 + edge(lv[se[z,1]] ,lv[se[z,2]])
          E(g4)$weight[c] = se[z,3]
          c = c + 1
        }
      }

    }
  }
  g4 = simplify(g4)
  # plot(g4,layout=layout.reingold.tilford(g4, circular=F))
  # pieTreeRoot(mst=g4, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  ###### 4. Merge isolated nonleaf node into the nearest edges
  nonleafroot = nonleaf#[-which(nonleaf==root)]
  dg = degree(g4, v = nonleafroot)
  id = which(dg==0)

  g5 = g4
  while (length(id)>0) {
    isonode = nonleafroot[id]
    if(length(isonode)==1){

    }
    curnode = nonleaf[ !(nonleaf %in% isonode) ]
    #rds = cenDist[root,curnode] #distances(g0, v=root, to=curnode)
    #ids = cenDist[root,isonode] #distances(g0, v=root, to=isonode)
    ids =  cenDist[root,isonode]
    k = which.min(ids)

    v1 = isonode[k]
    ds = rootDist[v1,curnode] #distances(g0, v=v1, to=curnode)
    od = order(ds)
    curnode = curnode[od]
    j=0
    if(j==0){
      nb = curnode[1]
      ns = neighbors(g5,v=nb) # get its neighbors
      ns = as.vector(ns)
      ns = ns[!(ns %in% c(isonode,v1))]
      rs = rootDist[v1, ns]  #distances(g1, v=v1, to= ns )
      od = order(rs)
      ns = ns[od]
      j = 1
    }

    if(length(ns)>0){
      ## add edge between v1 and ns[j]
      g5 = g5 + edge(v1, nb)
      E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
      g5 = g5 + edge(v1, ns[j])
      E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,ns[j]] # distances(g0,v=v1,to= ns[j])[1,1]
      edges = paste(nb,"|",ns[j],sep = "")
      g5 = delete_edges(g5, edges)
      es = distances(g5,v=nb,to=ns[j])
    }else{
      print("Error")
      if(nb %in% allleaf){## add nonleaf node into nearest leaf
        #zs = distances(g5, nb, nonleaf)
        zs = cenDist[nb,nonleaf]

        nn = nonleaf[which.min(zs[1,])]
        pth =  get.shortest.paths(g5, from=nn, to=nb)
        pth = as.vector(pth$vpath[[1]])
        nb = pth[2]

        ## add edge between v1 and ns[j]
        g5 = g5 + edge(v1, nb)
        E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
        g5 = g5 + edge(v1, nn)
        E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nn] # distances(g0,v=v1,to= ns[j])[1,1]
        edges = paste(nb,"|",nn,sep = "")
        g5 = delete_edges(g5, edges)

      }else{
        print("Error")
        g5 = g5 + edge(v1, nb)
        E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
      }
    }

    dg = degree(g5, v = nonleafroot)
    id = which(dg==0)

  }
  g5 = simplify(g5)
  #  plot(g5,layout=layout.reingold.tilford(g5, circular=F), vertex.size=2,  vertex.label.cex=2)
  # pieTreeRoot(mst=g5, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  list(g5=g5,rootDist=rootDist,w=sum(E(g5)$weight))
  #gr=list(g5=g5,rootDist=rootDist)
}


#' simpPTime: Infer pseduo time based on minimum spanning tree (MST) and given root in ISOMAP space
#'
#' @param X, the ISOMAP data;
#'
#' @param cluid, the clustering id,
#'
#' @param cenTree, the tree for mapping the cells. If it is null, use the clusters to build MST
#' @param root, the root node of MST
#' @return
#'    Y, the nodes to build MST;
#'    cenTree, the MST;
#'    pt, pseduo time;
#'    cellDs, the list
#' @export
#'
#' @examples

simpPTime<-function ( X, cluid, cenTree=NULL, root = 1)
{

  Y <- getCenters(X, cluid)
  print(rownames(Y))
  if(is.null(cenTree)){
    cenTree <- getCenMST(Y)
  }

  allcors = matrix(0, nrow = nrow(Y), ncol = 2)
  allpt = matrix(0, nrow = nrow(Y), ncol = nrow(X))
  if (degree(cenTree, paste("Y_", root, sep = "")) >= 1) {
    allres = newOrder(X, Y, dp_mst=cenTree, cluid, root = root)
    pt = allres[, 7]
  }
  return(list(Y = Y, cenTree = cenTree, pt = pt, cellDs = allres))
}


#' newOrder: Infer pseduo time based on minimum spanning tree (MST) and given root in ISOMAP space
#'
#' @param X, the ISOMAP data;
#'
#' @param Y, the cluster centers
#'
#' @param dp_mst, the tree for mapping the cells. If it is null, use the clusters to build MST
#'
#' @param cluid, the cluster id
#'
#' @param root, the root node of MST
#'
#' @return
#'    ptime, matrix contians the pseudotime and distance from cell to the tree
#' @export
#'
#' @examples
newOrder <- function (X, Y, dp_mst = NULL, cluid, root = 1)
{
  if (is.null(dp_mst)) {
    dp <- as.matrix(dist(Y))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
  }
  if (root < 0 || root > nrow(Y)) {
    stop("Error: Not existing such root!!!")
  }
  root_state = paste("Y_", root, sep = "")
  edMat = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
  edist = rep(-1, nrow(edMat))
  for (i in 1:nrow(edMat)) {
    d1 = distances(dp_mst, root_state, edMat[i, 1])
    d2 = distances(dp_mst, root_state, edMat[i, 2])
    v1 = as.integer(strsplit(edMat[i, 1], split = "_")[[1]][2])
    v2 = as.integer(strsplit(edMat[i, 2], split = "_")[[1]][2])
    edist[v1] = d1
    edist[v2] = d2
    if (d1 > d2) {
      edMat[i, 1] = paste("Y_", v2, sep = "")
      edMat[i, 2] = paste("Y_", v1, sep = "")
    }
  }
  cellNum = nrow(X)
  eNum = nrow(edMat)
  ptime = matrix(0, nrow = cellNum, ncol = 8)
  distances_X_to_Y <- proxy::dist(X, Y)
  distances_Y <- proxy::dist(Y, Y)
  ##################################
  ##For each cluster, map the cells in the cluster to the its connected edges
  nearestEdge = 1:nrow(X)

  for(i in 1:nrow(Y)){
    nb = as.vector(neighbors(dp_mst, v=i))
    kid = which(cluid==i)
    curC = X[kid,]
    dtCellEdge = matrix(0, nrow = nrow(curC), ncol = length(nb))
    dtCellProjL = matrix(0, nrow = nrow(curC), ncol = length(nb))
    dtCellPLR = matrix(0, nrow = nrow(curC), ncol = length(nb))

    ### compute the
    for(j in 1:length(nb)){
      if(distances(dp_mst, root_state, i) < distances(dp_mst, root_state, nb[j])){
        O = Y[i, ]
        B = Y[nb[j], ]
      }else{
        O = Y[nb[j], ]
        B = Y[i, ]
      }

      OB = B - O
      md = apply(curC, 1, function(x) {
        ds_p2OB(x, O, B)
      })
      dtCellEdge[, j] = md[1, ]
      dtCellProjL[, j] = md[2, ]
      dtCellPLR[, j] = md[5, ]
    }
    stime = matrix(0, nrow = nrow(curC), ncol = 6)

    cellNE = apply(dtCellEdge, 1, which.min)

    for (k in 1:nrow(curC)) {
      v1 = nb[cellNE[k]]
      if(distances(dp_mst, root_state, i) < distances(dp_mst, root_state, v1)){
        stime[k, 1] =  i
        stime[k, 2] =  v1
        stime[k, 3] = distances(dp_mst, root_state, i) + dtCellProjL[k, cellNE[k]]
      }else{
        stime[k, 1] =  v1
        stime[k, 2] =  i
        stime[k, 3] = distances(dp_mst, root_state, v1) + dtCellProjL[k, cellNE[k]]
      }
      stime[k, 4] = dtCellEdge[k, cellNE[k]]
      stime[k, 5] = dtCellPLR[k, cellNE[k]]
      stime[k, 6] = dtCellProjL[k, cellNE[k]]
    }
    ptime[kid,1:ncol(stime)] = stime
  }

  ##########################################

  for (i in 1:eNum) {
    v1 = as.integer(strsplit(edMat[i, 1], split = "_")[[1]][2])
    v2 = as.integer(strsplit(edMat[i, 2], split = "_")[[1]][2])
    cid = which(ptime[,1]==v1 & ptime[,2]==v2)
    id1 = cid[which(cluid[cid]==v1)]
    id2 = cid[which(cluid[cid]==v2)]
    ds1 = ptime[id1, 6]
    ds2 = ptime[id2, 6]
    ds1 = ds1 - min(ds1)
    ds2 = ds2 - min(ds2)
    L = as.numeric(edMat[i,3])
    p0 = L/(max(ds1)+max(ds2))
    p1 = max(ds1)/(max(ds1)+max(ds2))
    rs = distances(dp_mst, root_state, v1)[1,1];
    ptime[id1,7] = ds1*p0 + rs
    ptime[id2,7] = max(ds1)*p0  + rs + ds2*p0
    ptime[id1,8] = ds1*p0
    ptime[id2,8] = max(ds1)*p0 + ds2*p0

    # summary(ds1*p0 + rs)
    # summary(max(ds1)*p0   + ds2*p0 + rs )
    # table(cluid[id1])
    # table(cluid[id2])
  }

  colnames(ptime) = c("c1", "c2", "dist2Root2", "dtCellEdge",
                      "dtCellPSign", "projLen", "dist2Root2", "dist2Parent")
  ptime
}

svmCluster <- function(X,i=1,j=2,cluid){
  id1 <- which(cluid==i)
  id2 <- which(cluid==j)
  id <- c(id1,id2)
  sZ = X[id, ]
  nc <- c(rep(1, length(id1)), rep(2, length(id2)))
  dat = data.frame(x=sZ,y= as.factor(nc))
  svmfit = svm(y ~ ., data = dat, kernel = "linear", cost = 10, scale = FALSE)
  pred <- predict(svmfit, dat)
  a = table(as.integer(pred),nc)
  if(nrow(a)>1){
    svmInd  <- a[1,2] + a[2,1]
  }else{
    svmInd  <- a[1,2]
  }
  return(svmInd)
}




#' getCluCellTh: get the main cell types in the cluster
#'
#' @param cluid
#'     cluster id
#' @param cell_labels
#'     cell annotations
#' @param th
#'     threshold for the proportion of the cell type in the cluster 0.1,
#' @param subN
#'     threshold of the number of character in the string
#' @return
#'     main cell types in each cluster
#' @export
#'
#' @examples
#'
getCluCellTh <- function(cluid, cell_labels, th=0.1, subN=4){
  xc <- c(cluid, rep("xx", max(cluid)))
  tb = table(cluid)
  tn = rep("", length(tb))
  for(i in 1:length(tb)){
    id = which(cluid==i)
    td = table(cell_labels[id])
    td = td/length(id)
    id = order(td, decreasing = TRUE)
    tn[i]=i
    for(j in 1:length(id)){
      if(td[id[j]] >= th){
        sbn = substr(names(td[id[j]]), 1,subN)
        tn[i] = paste(tn[i], sbn, sep="_")
      }
    }
  }
  tb = tn
}

#' getCellSM: the matrix of cell number from each cell type in each cluster
#'
#' @param cluL
#'     cluster id
#' @param cellft
#'     cell annotations
#' @return
#'     the matrix of cell number from each cell type in each cluster
#' @export
#'
#' @examples
#'
getCellSM <- function(cluL, cellft){
  ctyM = matrix(0, nrow = length(table(cluL)), ncol = length(table(cellft)))
  cluN = as.integer(names(table(cluL)))
  cellN = as.character(names(table(cellft)))
  for(i in 1:nrow(ctyM)){
    cid = which(cluL==cluN[i])
    subcell = cellft[cid]
    for(j in 1:ncol(ctyM)){
      ctyM[i,j] = length(which(subcell==cellN[j]))
    }
  }
  colnames(ctyM) =  cellN
  rownames(ctyM) =  cluN
  ctyM
}


#' distGraphClu: Compute the distance among the clusters in knn graph
#'
#' @param Z
#'     matrix, the selected PCs
#' @param cluid
#'     vector, cell clustering results
#' @param ksize
#'     integer, neighbor siez to build knn graph
#' @return
#'     matrix, the mean graph distance between each two clusters in graph
#' @export
#'
#' @examples
#'     dsM = distGraphClu(Z, cluid,ksize=4)
#'
distGraphClu <- function(Z, cluid,ksize=4){
  message(Sys.time(), ": Get neighbors ")

  knng <- makeKNNgraph(x = Z, k = ksize, eps = 0)
  nc = components(knng$g, mode =  "strong" )
  while(nc$no>1){
    print("ksize is too small to build a connected graph!")
    ksize = ksize + 1
    knng <- makeKNNgraph(x = Z, k = ksize, eps = 0)
    nc = components(knng$g, mode =  "strong" )
    print(ksize)
  }
  mG = knng$g

  graphDist = distances(mG)
  clus = unique(cluid)
  ## the sum of distance of neighbors in each two clusters
  dsM = matrix(0, nrow = length(clus), ncol = length(clus))

  for(i in 2:length(clus)){
    id1 = which(cluid == clus[i])
    for(j in 1:(i-1)){
      id2 = which(cluid == clus[j])
      cs = mean(graphDist[id1, id2])
      dsM[clus[i], clus[j]] = cs
      dsM[clus[j], clus[i]] = cs
    }
  }

  ##
  g0 = kNNDist(dsM, ks=ksize)
  g0 = simplify(g0)
  #plot(g0,vertex.label.cex = 1.5,vertex.size=0.1, bg="white")
  message(Sys.time(), ": End ")
  dsM
}


#' getCluNeighbors: Compute the neighbors number for each two clusters, if the number is larger than 0,
#'                  compute mean distance between two clusters.
#'
#' @param Z
#'     matrix, the selected PCs
#' @param cluid
#'     vector, cell clustering results
#' @param ksize
#'     integer, neighbor siez to build knn graph
#' @return
#'     matrix, the mean graph distance between each two clusters in graph
#' @export
#'
#' @examples
#'     dsM = getCluNeighbors(Z, cluid,ksize=4)
#'
getCluNeighbors <- function(Z,  cluid, ksize=4){
  #
  message(Sys.time(), ": Get neighbors ")

  clus = unique(cluid)
  mG = RANN::nn2(Z, k= ksize, searchtype="priority")

  #nnL = list()#store the neigbors for each two cluster
  #ndL = list()#store the distance vector for each two cluster

  ## storage the neighbors number for each two clusters
  ## If there is no neighbors, set it as -1
  nnM = matrix(0, nrow = length(clus), ncol = length(clus))
  ## the sum of distance of neighbors in each two clusters
  dsM = matrix(0, nrow = length(clus), ncol = length(clus))


  for(i in 1:length(clus)){
    fc = clus[i]
    scid = which(cluid==fc)
    nbid = as.vector(mG$nn.idx[scid,2:ncol(mG$nn.idx)])
    subns = as.vector(mG$nn.dists[scid,2:ncol(mG$nn.idx)])
    cc = cluid[nbid]
    nclu = unique(cc)

    for(j in 1:length(nclu)){
      tc = nclu[j]
      kid = which(cc==tc)
      nnM[fc,tc] = nnM[fc,tc] +  length(kid)
      dsM[fc,tc] = dsM[fc,tc] + summary(subns[kid])[4]
    }
  }


  for(i in 1:nrow(nnM)){
    id = which(nnM[i,]==0)
    nnM[i,id] = -1
  }

  msM = dsM/nnM

  for(i in 1:nrow(msM)){
    id = which(msM[i,]<=0)
    msM[i,id] = Inf
  }
  diag(msM) = 0
  for(i in 1:nrow(msM)){
    for(j in 1:nrow(msM)){
      mv = min(msM[i,j], msM[j,i])
      msM[i,j] = mv
      msM[j,i] = mv
    }
  }

  ##
  #g0 = kNNDist(msM, ks=k)
  #g0 = simplify(g0)
  #plot(g0,vertex.label.cex = 1.5,vertex.size=0.1, bg="white")
  message(Sys.time(), ": End ")
  msM
}


#' mergeCluster: merge clusters given specified cluster id
#'
#' @param mc
#'     vector the cluster id to merge
#' @param cluid
#'     vector, cell clustering results
#' @return
#'     vector, the new merged cluster id
#' @export
#'
#' @examples
#'     newCluid = getCluNeighbors(mc=c(2,4), cluid)
#'
mergeCluster <- function(mc=c(2,4), cluid){
  rc =NULL
  for(i in 1:length(mc)){
    id1 = which(cluid==mc[i])
    rc = c(rc, id1)
  }
  cluid[rc] = 0
  tc = unique(cluid)
  tc = tc[order(tc)]
  lc = length(tc)
  mxv = max(tc)
  for(i in 1:lc){
    cluid[which(cluid==tc[i])] = i + mxv
  }
  cluid = cluid -  mxv
  cluid
}



#' getMSTLeafRootISO: create spanning tree with specified root and leaves
#' @param uY
#'     UMAP/t-SNE matrix
#' @param Z
#'     pca matrix
#' @param dsM
#'     distance matrix of clusters
#' @param ck
#'    integer, used to build knn graph for umap
#' @param ksize
#'    integer, default is 2, to build knn graph for clusters
#' @param cluid
#'     vector, the cluster id
#' @param root
#'    integer, default is 1, the root cluster
#' @param leaves
#'    list, specify the leave groups
#' @return graph, g1, the trajectory.
#'
#' @export
#'
#' @examples
#'
buildG_UMAP_PCA <- function(uY, Z, dsM, ck =8, ksize = 4, cluid=cluid, root=rt,
                            leaves = list(a= 6 , b=8, c=15, d= 14,e=11, f=13,g=12 ) ){

  #
  message(Sys.time(), ": Build graph with k neighbors ")
  clus = unique(cluid)

  knng <- makeKNNgraph(x = uY, k = ck, eps = 0)
  n = components(knng$g, mode =  "strong" )
  while(n$no >= length(clus)){
    #print("ksize is too small to build a connected graph!")
    ck = ck + 2
    knng <- makeKNNgraph(x = uY, k = ck, eps = 0)
    n = components(knng$g, mode =  "strong" )
  }
  message(Sys.time(), paste0(": There are ", n$no, " components with ksize=", ck))
  mG = knng$g

  #plotAlls(X=uY,draw=1, t=NULL, cluid=as.character(n$membership), type=2,size =6,opacity=1,
  #         titlename=paste("Components"))

  ##1. Use neighbor distance of clusters in graph bulit by PCA, but not good
  #msM = getCluNeighbors(Z, k=ck, cluid=cluid)


  ####
  rootid = which(cluid == rt)
  tb = table(n$membership[rootid])
  rootC = as.integer(names(which.max(tb)))

  nc = unique(n$membership)

  leftGC = nc[-which(nc==rootC)]

  cluMeb = matrix(0, nrow = n$no, ncol= length(clus))

  for(i in 1:length(leftGC)){
    cid = which(n$membership==leftGC[i])
    ct = table(cluid[cid])
    cid = as.integer(names(ct))
    cluMeb[leftGC[i], cid] = as.integer(ct)
  }

  cluR = matrix(0, nrow = length(clus), ncol = length(clus))
  ncStat = matrix(0, nrow = n$no, ncol = length(clus))
  for(i in 1:n$no){
    rid = which(n$membership==i)
    ct = table(cluid[rid])
    cid = as.integer(names(ct))
    ncStat[i,cid] = ct
  }

  #### get each subgraph in each component
  for(i in 1:n$no){
    rid = which(n$membership==i)
    ct = table(cluid[rid])
    cid = as.integer(names(ct))
    ef = rep(FALSE, length(cid))
    for(b in 1:length(cid)){
      if(ncStat[i, cid[b]]==max(ncStat[, cid[b] ])){
        ef[b] = TRUE
      }
    }
    cid = cid[ef]

    if(length(cid)==1){
      rs = dsM[cid,]
      rs[cid] = Inf
      od = order(rs)
      for(w in 1:length(od)){
        st = od[w]
        if(dsM[rt, st] < dsM[rt, cid] &  isLeafNode(st, leaves)==0){
          break
        }
      }
      if(length(which(cluR[cid,]==1))==0){
        cluR[cid, st] = 1
        cluR[st, cid] = 1
      }
      next
    }

    ## the sum of distance of neighbors in each two clusters
    ldsM = matrix(0, nrow = length(cid), ncol = length(cid))
    rownames(ldsM) = cid
    colnames(ldsM) = cid
    for(s in 2:length(cid)){
      id1 = which(cluid == cid[s] &  n$membership==i)
      for(j in 1:(s-1)){
        id2 = which(cluid == cid[j]  &  n$membership==i)
        cs = mean(distances(mG,id1, id2))
        ldsM[s, j] = cs
        ldsM[j, s] = cs
      }
    }
    ### select the cluster which is nearest to the root based PCA
    mid = which.min(dsM[rt, cid])

    localRoot = mid
    pid = cluid[rid]
    nid = rep(0, length(pid))
    for(k in 1:length(cid)){
      nid[which(pid==cid[k])] = k
    }

    leafID = sapply(cid, function(x){ isLeafNode(x, leaves) })

    localLeaf = list()
    k = 1
    for(s in 1:length(leafID)){
      if(leafID[s] > 0){
        localLeaf[[k]] = s
        k = k +1
      }
    }

    if(length(localLeaf) == 0 | length(ct) < 3 ){
      rownames(ldsM) = 1:length(cid)
      colnames(ldsM) = 1:length(cid)
      gp <- graph.adjacency(ldsM, mode = "undirected", weighted = TRUE)
      g1 = mst(gp)

    }else if(length(localLeaf) < length(leaves)){
      gr = getLeafRootG(rootDist = ldsM, cenDist=ldsM, subcluid=nid, L2=min(ksize, nrow(ldsM)-1),
                        root=localRoot,  localLeaf= localLeaf)

      g1 = gr$g1
      #plot(g1,layout=layout.reingold.tilford(g1, circular=F), vertex.label.cex = 1.5,vertex.size=0.1)

    }else{
      g1 = getMSTLeafRootISO(rootDist = ldsM, cenDist=ldsM, cluid=nid, L2=min(ksize, nrow(ldsM)-1),
                             root=localRoot, leaves= localLeaf)
      g1 = g1$g5

    }
    #plot(g1,layout=layout.reingold.tilford(g1, circular=F), vertex.label.cex = 1.5,vertex.size=0.1)

    edgeL = get.edgelist(g1)
    for(e in 1:nrow(edgeL)){
      st = cid[as.integer(edgeL[e,1]) ]
      ed = cid[as.integer(edgeL[e,2])]
      if(ncStat[i, st] == max(ncStat[,st]) & ncStat[i, ed] == max(ncStat[,ed]) ){
        cluR[st, ed] = 1
        cluR[ed, st] = 1
      }
    }
  }

  g0 = graph_from_adjacency_matrix(cluR, mode = "undirected")
  g0 = simplify(g0)
  #plot(g0,vertex.label.cex = 2,vertex.size=0.01)


  sn = components(g0, mode =  "strong" )
  if(sn$no>1){# connnect each two components to its parent
    rc = sn$membership[root]
    ac = 1:sn$no
    ac = ac[-rc]
    for(h in ac){
      rcid = which(sn$membership==h)
      leafID = sapply(rcid, function(x){ isLeafNode(x, leaves) })
      if(length(which(leafID==0)) > 0) {
        rcid = rcid[leafID==0]
      }

      localRt = rcid[which.min(dsM[rt, rcid])]
      ##Find its parent node
      od = order(dsM[localRt,])
      for(m in od){
        if(m==localRt | isLeafNode(m, leaves) | sn$membership[m] == sn$membership[localRt]){

        }else{
          break
        }
      }

      cluR[localRt, m] = 1
      cluR[m, localRt] = 1

    }
    g0 = graph_from_adjacency_matrix(cluR, mode =   "undirected")
    g0 = simplify(g0)
    plot(g0,vertex.label.cex = 2,vertex.size=0.01)
  }

  ### Do mst for each leaves
  for(j in 1:length(leaves)){
    if(length(leaves[[j]]) > 1){
      lv = leaves[[j]]
      B <- as.matrix(dsM[lv,lv] )
      gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
      dp_mst = mst(gp)
      se = get.edgelist(dp_mst)
      for(z in 1:nrow(se)){
        cluR[lv[se[z,1]] ,lv[se[z,2]]] = 1
        cluR[lv[se[z,2]] ,lv[se[z,1]]] = 1
      }

      ps = 1:length(lv)
      dss = 1:length(lv)
      for(h in 1:length(lv)){
        u = lv[h]
        x = which(cluR[u,] == 1)
        x = x[!(x %in% lv)]

        if(length(x)>0){
          cluR[u,x] = 0
          cluR[x,u] = 0
          lid = which.min(dsM[u, x])
          x = x[lid]

          ps[h] = x
          dss[h] = dsM[u, x]
        }else{
          ps[h] = 0
          dss[h] = Inf
        }
      }

      mid = which.min(dss)
      cluR[ps[mid],  lv[mid]] = 1
      cluR[lv[mid],  ps[mid]] = 1
    }
  }

  g0 = graph_from_adjacency_matrix(cluR, mode = "undirected")
  g0 = simplify(g0)
  #plot(g0,vertex.label.cex = 2,vertex.size=0.01)

  ## Add distance
  se =  get.edgelist(g0)
  k = 1
  for(z in 1:nrow(se)){
    E(g0)$weight[k] = dsM[se[z,1] , se[z,2] ]
    k = k + 1
  }

  ########## remove redudant edges
  g1 = getMSTLR(g0, rootDist=dsM, cenDist=dsM, cluid, root,  leaves)
  #plot(g1,vertex.label.cex = 2,vertex.size=0.01)
  g1
}




#' getLeafRootG: create spanning tree with specified root and leaves
#' @param rootDist
#'     similariy matrix of clusters
#' @param cenDist
#'     distance matrix of clusters
#' @param subcluid
#'     vector, the cluster id
#' @param L2
#'    integer, default is 2, knn size
#' @param root
#'    integer, default is 1, the root cluster
#' @param localLeaf
#'    list, specify the leave groups
#' @return list
#'      list(g1=g1,edList=edList)
#' @export
#'
#' @examples
#'
getLeafRootG <- function(rootDist, cenDist, subcluid, L2=2, root=1,
                         localLeaf = list(a=c(9,12,13), b=17,c=c(2,6),d=3,e=c(4,5),f=14))
{
  allleaf = NULL
  for(i in 1:length(localLeaf)){
    allleaf = c(allleaf, localLeaf[[i]])
  }
  nonleaf = (1:max(subcluid))[-allleaf]
  ########1.Calcualte the distance matrxi among each clusters
  ###1.Compute the minimum distance from the leaf cluster to the root cluster
  M = length(unique(subcluid))

  if(length(nonleaf)==1){
    ## convert into all nodes graph
    g1 <- igraph::make_empty_graph(nrow(rootDist), directed = FALSE)
  }else{

    ########2.Build kNN graph using the distance matrix of nonleaf nodes
    g0 = kNNDist(rootDist[nonleaf,nonleaf], ks=L2)
    dg = degree(g0, v=1:length(nonleaf))
    dg
    n = components(g0, mode =  "strong" )
    mark=TRUE
    nonleafroot = nonleaf[-which(nonleaf==root)]
    while(n$no>1){
      print("L2 is too small to build a connected graph!")
      L2 = L2 + 1
      g0 <- kNNDist(rootDist[nonleaf,nonleaf], ks=L2)
      n = components(g0, mode =  "strong" )
      print(L2)
    }
    #
    edList= cbind(get.edgelist(g0), E(g0)$weight)
    #plot(g0,layout=layout.reingold.tilford(g0, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

    ## convert into all nodes graph
    g1 <- igraph::make_empty_graph(nrow(rootDist), directed = FALSE)
    for(i in 1:nrow(edList)){
      v1 = nonleaf[edList[i,1]]
      v2 = nonleaf[edList[i,2]]
      g1 = g1 + edge(v1, v2)
      E(g1)$weight[i] = edList[i,3]
    }
  }

  #plot(g1,layout=layout.reingold.tilford(g1, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

  ########3. Add leaf node into the nearest nonleaf node

  len = length(localLeaf)
  ## Keep the degree of leaf node group to its parent is one
  for(i in 1:len){
    lv = localLeaf[[i]]
    ##### add edge among leave using MST
    if(length(lv)>1){
      B <- as.matrix(rootDist[lv,lv] )
      gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
      dp_mst = mst(gp)
      se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
      k = length(E(g1)$weight) + 1
      for(z in 1:nrow(se)){
        g1 = g1 + edge(lv[se[z,1]] ,lv[se[z,2]])
        E(g1)$weight[k] = se[z,3]
        k = k + 1
      }
    }
    ## Find nearest nonleaf node of the leaf group
    ls = lv
    ld = lv
    for(j in 1:length(lv)){
      ds = rootDist[lv[j], nonleaf]
      ls[j] = nonleaf[which.min(ds)]
      ld[j] = min(ds)
    }
    id = which.min(ld)
    k = length(E(g1)$weight) + 1
    g1 = g1 + edge(lv[id], ls[id])
    E(g1)$weight[k] = ld[id]
  }
  g1 = simplify(g1)
  edList= cbind(get.edgelist(g1), E(g1)$weight)
  print(nrow(edList))
  #plot(g1,layout=layout.reingold.tilford(g1, circular=F))
  # pieTreeRoot(mst=g1, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)

  list(g1=g1,edList=edList)
}


#' getMSTLR: create spanning tree with specified root and leaves
#' @param g1
#'     input graph
#' @param rootDist
#'     similariy matrix of clusters
#' @param cenDist
#'     distance matrix of clusters
#' @param cluid
#'     vector, the cluster id
#' @param root
#'    integer, default is 1, the root cluster
#' @param leaves
#'    list, specify the leave groups
#' @return g1, graph
#'
#' @export
#'
#' @examples
#'
getMSTLR<- function(g1, rootDist, cenDist, cluid, root=1,
                    leaves = list(a=c(9,12,13), b=17,c=c(2,6),d=3,e=c(4,5),f=14))
{
  allleaf = NULL
  for(i in 1:length(leaves)){
    allleaf = c(allleaf, leaves[[i]])
  }
  nonleaf = (1:max(cluid))[-allleaf]
  M = length(unique(cluid))

  len = length(leaves)
  ###### 1. get shortest path from root to leaf
  g4 = igraph::make_empty_graph(M, directed = FALSE)
  paths = list() #get.all.shortest.paths(g1,from = leaves, to=root, mode = "all")
  c=1
  for(i in 1:length(leaves)){
    lv = leaves[[i]]
    ds = distances(g1, v=lv, to= root)
    id = which.min(ds[,1])
    p1 = get.shortest.paths(g1, from=root, to=lv[id])
    if(length(p1$vpath)==0){
      message(Sys.time(), "Error: Please set larger L2 !!!")
      return()
    }
    ev = as.vector(p1$vpath[[1]])

    ##### add edge among leaves using MST
    if(length(ev)>1){
      h = length(ev)-1
      for(j in 1:h){
        ei <- get.edge.ids(g4, c(ev[j],ev[j+1]))
        if(ei==0){
          g4 = g4 + edge(ev[j],ev[j+1])
          E(g4)$weight[c] = distances(g1, v=ev[j], to= ev[j+1])[1,1]
          c = c + 1
        }
      }

      if(length(lv) > 1){
        B <- as.matrix(rootDist[lv,lv] )
        gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
        dp_mst = mst(gp)
        se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
        for(z in 1:nrow(se)){
          g4 = g4 + edge(lv[se[z,1]] ,lv[se[z,2]])
          E(g4)$weight[c] = se[z,3]
          c = c + 1
        }
      }

    }
  }
  g4 = simplify(g4)
  # plot(g4,layout=layout.reingold.tilford(g4, circular=F))
  # pieTreeRoot(mst=g4, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  ###### 4. Merge isolated nonleaf node into the nearest edges
  nonleafroot = nonleaf#[-which(nonleaf==root)]
  dg = degree(g4, v = nonleafroot)
  id = which(dg==0)

  g5 = g4
  while (length(id)>0) {
    isonode = nonleafroot[id]
    if(length(isonode)==1){

    }
    curnode = nonleaf[ !(nonleaf %in% isonode) ]
    #rds = cenDist[root,curnode] #distances(g0, v=root, to=curnode)
    #ids = cenDist[root,isonode] #distances(g0, v=root, to=isonode)
    ids =  cenDist[root,isonode]
    k = which.min(ids)

    v1 = isonode[k]
    ds = rootDist[v1,curnode] #distances(g0, v=v1, to=curnode)
    od = order(ds)
    curnode = curnode[od]
    j=0
    if(j==0){
      nb = curnode[1]
      ns = neighbors(g5,v=nb) # get its neighbors
      ns = as.vector(ns)
      ns = ns[!(ns %in% c(isonode,v1))]
      rs = rootDist[v1, ns]  #distances(g1, v=v1, to= ns )
      od = order(rs)
      ns = ns[od]
      j = 1
    }

    if(length(ns)>0){
      ## add edge between v1 and ns[j]
      g5 = g5 + edge(v1, nb)
      E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
      g5 = g5 + edge(v1, ns[j])
      E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,ns[j]] # distances(g0,v=v1,to= ns[j])[1,1]
      edges = paste(nb,"|",ns[j],sep = "")
      g5 = delete_edges(g5, edges)
      es = distances(g5,v=nb,to=ns[j])
    }else{
      print("Error")
      if(nb %in% allleaf){## add nonleaf node into nearest leaf
        #zs = distances(g5, nb, nonleaf)
        zs = cenDist[nb,nonleaf]

        nn = nonleaf[which.min(zs[1,])]
        pth =  get.shortest.paths(g5, from=nn, to=nb)
        pth = as.vector(pth$vpath[[1]])
        nb = pth[2]

        ## add edge between v1 and ns[j]
        g5 = g5 + edge(v1, nb)
        E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
        g5 = g5 + edge(v1, nn)
        E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nn] # distances(g0,v=v1,to= ns[j])[1,1]
        edges = paste(nb,"|",nn,sep = "")
        g5 = delete_edges(g5, edges)

      }else{
        print("Error")
        g5 = g5 + edge(v1, nb)
        E(g5)$weight[length(E(g5)$weight)] = rootDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
      }
    }

    dg = degree(g5, v = nonleafroot)
    id = which(dg==0)

  }
  g5 = simplify(g5)
  #  plot(g5,layout=layout.reingold.tilford(g5, circular=F), vertex.size=2,  vertex.label.cex=2)
  # pieTreeRoot(mst=g5, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  #list(g5=g5,rootDist=rootDist,w=sum(E(g5)$weight))
  #gr=list(g5=g5,rootDist=rootDist)
  g5

}

#' getMSTLeafRootISO_3: create spanning tree with specified root and leaves
#' @param neigDist
#'     neighbor distance matrix of clusters
#' @param cenDist
#'     distance matrix of cluster centers
#' @param cluid
#'     vector, the cluster id
#' @param L2
#'    integer, default is 2, knn size
#' @param root
#'    integer, default is 1, the root cluster
#' @param leaves
#'    list, specify the leave groups, for each leave group, make it as a vector,
#'    For example,a=c(9,12,13,"linear"), b=c(2,5, "mst"), c=c(7,8,9,"parallel")
#' @return list
#'    list(g5=g5,neigDist=neigDist,w=sum(E(g5)$weight))
#' @export
#'
#' @examples
#'
getMSTLeafRootISO_3 <- function(neigDist, cenDist, cluid, L2=2, root=1,
                                leaves = list(a=c(9,12,13,"linear"), b=17,c=c(2,6),d=3,e=c(4,5),f=14))
{
  allleaf = NULL
  for(i in 1:length(leaves)){
    cl = leaves[[i]]
    if(length(leaves[[i]]) > 1) {
      allleaf = c(allleaf, cl[1:(length(cl)-1)])
    }else{
      allleaf = c(allleaf, cl)
    }
  }

  allleaf = as.integer(allleaf)
  nonleaf = (1:max(cluid))[-allleaf]
  nonleafroot = nonleaf[-which(nonleaf %in% root)]

  ###1.Build kNN graph using the distance matrix of all nodes
  M = length(unique(cluid))
  if(length(nonleaf)==1){
    ## convert into all nodes graph
    g1 <- igraph::make_empty_graph(nrow(cenDist), directed = FALSE)
  }else{
    ######## 1.1 Build kNN graph using the distance matrix of nonleaf nodes
    g0 = kNNDist(cenDist, ks=L2)
    dg = degree(g0, v=1:nrow(cenDist))
    dg
    n = components(g0, mode =  "strong" )
    mark=TRUE
    while(n$no>1){
      print("L2 is too small to build a connected graph!")
      L2 = L2 + 1
      g0 <- kNNDist(cenDist, ks=L2)
      n = components(g0, mode =  "strong" )
      print(L2)
    }
    #
    edList= cbind(get.edgelist(g0), E(g0)$weight)
    #plot(g0,layout=layout.reingold.tilford(g0, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

    for(i in 1:nrow(edList)){###1.2 Remove connection between leaf and nonleaf nodes
      if(edList[i,1] %in% allleaf  | edList[i,2] %in% allleaf){
        g0 = delete_edges(g0, edges=paste0(edList[i,1], "|", edList[i,2]))
      }
    }
    #plot(g0,layout=layout.reingold.tilford(g0, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)
    g1 <- g0
  }
  #plot(g1,layout=layout.reingold.tilford(g1, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

  ########2. Add leaf node into the nearest nonleaf node
  ## Keep the degree of leaf node group to its parent is one
  len = length(leaves)
  for(i in 1:len){
    lv = leaves[[i]]
    if(length(lv)>1){
      typeC = lv[length(lv)]
      sv = as.integer(lv[1:(length(lv)-1)])
    }else{
      sv = as.integer(lv)
    }
    ## 1. Find nearest nonleaf node of the leaf group and add the edge
    ls = sv
    ld = sv
    for(j in 1:length(sv)){
      ds = neigDist[sv[j], nonleaf]
      ls[j] = nonleaf[which.min(ds)]
      ld[j] = min(ds)
    }
    pid = which.min(ld)
    k = length(E(g1)$weight) + 1
    g1 = g1 + edge(sv[pid], ls[pid])
    E(g1)$weight[k] = cenDist[sv[pid], ls[pid]] # use cluster distance

    ##### add edge among leave using MST
    if(length(sv)>1){
      ## 2. Add nonleaf node to its children
      if(typeC=="parallel"){
        k = length(E(g1)$weight) + 1
        for(z in 1:length(sv)){
          g1 = g1 + edge(sv[z], ls[pid])
          E(g1)$weight[k] = cenDist[sv[z], ls[pid]] # use cluster distance
          k = k + 1
        }
      }else if(typeC=="linear"){
        B <- cenDist[sv,root]
        od = order(B)
        k = length(E(g1)$weight) + 1
        for(z in 1:(length(od)-1)){
          ks = sv[od[z]]
          ke = sv[od[z+1]]
          g1 = g1 + edge(ks, ke)
          E(g1)$weight[k] = cenDist[ks, ke] # use cluster distance
          k = k + 1
        }
      }else if(typeC=="mst"){
        B <- as.matrix(neigDist[sv,sv] )
        gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
        dp_mst = mst(gp)
        se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)

        k = length(E(g1)$weight) + 1
        for(z in 1:nrow(se)){
          g1 = g1 + edge(sv[se[z,1]] ,sv[se[z,2]])
          E(g1)$weight[k] = cenDist[sv[se[z,1]], sv[se[z,2]]] # use cluster distance
          #E(g1)$weight[k] = se[z,3]
          k = k + 1
        }
      }else{
        stop("Wrong leaf types specified")
        return()
      }

    }
  }
  g1 = simplify(g1)
  edList= cbind(get.edgelist(g1), E(g1)$weight)
  print(nrow(edList))
  # plot(g1,layout=layout.reingold.tilford(g1, circular=F))
  # pieTreeRoot(mst=g1, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)

  ####3.If the graph is not fully connected, connect the nonleaf nodes to its parent
  for(i in 1:length(nonleafroot)){
    nb = as.integer(neighbors(g1, nonleafroot[i]))
    tf = nb %in% allleaf
    #if the nonleaf node is connected with only leaf node, select its nearest node
    if(length(which(tf==TRUE)) == length(tf)){
      nk =  nonleafroot[i]
      localDS = cenDist[nk, ]
      localOrd = order(localDS)
      ef = !(localOrd %in% c(nk,allleaf))
      leftID = localOrd[ef]
      k = length(E(g1)$weight) + 1
      g1 = g1 + edge(nk,leftID[1])
      E(g1)$weight[k] = cenDist[nk,leftID[1]]
    }
  }

  ###If the component is larger than 1
  ng = components(g1, mode =  "strong" )
  if(ng$no>1){
    cn = ng$membership
    for(i in 2:ng$no){
      l1 = which(ng$membership==i)
      l1 = l1[-which(l1 %in% allleaf)]
      for(j in 1:(i-1)){
        l2 = which(ng$membership==j)
        l2 = l2[-which(l2 %in% allleaf)]
        p1 = l1
        w1 = l1
        if(length(l1) > 0 & length(l2) > 0){
          for(z in 1:length(l1)){
            km = which.min(cenDist[l1[z], l2])
            p1[z] = l2[km]
            w1[z] = min(cenDist[l1[z], l2])
          }
          se = l1[which.min(w1)]
          et = p1[which.min(w1)]
          k = length(E(g1)$weight) + 1
          g1 = g1 + edge(se, et)
          E(g1)$weight[k] = cenDist[se,et]
        }
      }
    }
  }
  # plot(g1,layout=layout.reingold.tilford(g1, circular=F))

  ###### 4. Find a path from root to leaf
  g4 = igraph::make_empty_graph(M, directed = FALSE)
  paths = list() #get.all.shortest.paths(g1,from = leaves, to=root, mode = "all")
  c=1
  for(i in 1:length(leaves)){
    lv = leaves[[i]]
    if(length(lv) > 1) {
      typeC = lv[length(lv)]
      sv = as.integer(lv[1:(length(lv)-1)])
    }else{
      sv = as.integer(lv)
    }
    ##4.1 Get the path from each leaf to root
    ds = distances(g1, v=sv, to= root)
    id = which.min(ds[,1])
    p1 = get.shortest.paths(g1, from=root, to=sv[id])
    if(length(p1$vpath)==0){
      message(Sys.time(), "Error: Please set larger L2 !!!")
      return()
    }
    ev = as.vector(p1$vpath[[1]])

    if(length(ev)>1){
      ## 4.2 Add the edge into the graph from root to leaf
      h = length(ev)-1
      for(j in 1:h){
        ei <- get.edge.ids(g4, c(ev[j],ev[j+1]))
        if(ei==0){
          g4 = g4 + edge(ev[j],ev[j+1])
          E(g4)$weight[c] = distances(g1, v=ev[j], to= ev[j+1])[1,1]
          c = c + 1
        }
      }

      if(length(sv) > 1){#4.3 Add edge among the leaf
        if(typeC=="parallel"){# add edge between parent and leaf
          pid = ev[length(ev)-1]
          for(z in 1:length(sv)){
            if(sv[z] != pid){
              g4 = g4 + edge(pid, sv[z] )
              E(g4)$weight[c] = cenDist[pid,sv[z] ]
              c = c + 1
            }
          }
        }else if(typeC=="linear"){
          B <- cenDist[sv,root]
          od = order(B)
          for(z in 1:(length(od)-1)){
            ks = sv[od[z]]
            ke = sv[od[z+1]]
            g4 = g4 + edge(ks, ke)
            E(g4)$weight[c] = cenDist[ks, ke] # use cluster distance
            c = c + 1
          }
        }else if(typeC=="mst"){
          B <- as.matrix(neigDist[sv,sv] )
          gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
          dp_mst = mst(gp)
          se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
          for(z in 1:nrow(se)){
            g4 = g4 + edge(sv[se[z,1]], sv[se[z,2]])
            E(g4)$weight[c] = cenDist[sv[se[z,1]], sv[se[z,2]]]  #se[z,3]
            c = c + 1
          }
        }else{
          stop("Wrong leaf types specified")
          return()
        }
      }
    }
  }
  g4 = simplify(g4)
  # plot(g4,layout=layout.reingold.tilford(g4, circular=F))
  # pieTreeRoot(mst=g4, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  ###### 4. Merge isolated nonleaf node into the nearest edges
  nonleafroot = nonleaf#[-which(nonleaf==root)]
  dg = degree(g4, v = nonleafroot)
  id = which(dg==0)

  g5 = g4
  while (length(id)>0) {
    isonode = nonleafroot[id]
    if(length(isonode)==1){

    }
    curnode = nonleaf[ !(nonleaf %in% isonode) ]
    #rds = cenDist[root,curnode] #distances(g0, v=root, to=curnode)
    #ids = cenDist[root,isonode] #distances(g0, v=root, to=isonode)
    ids =  cenDist[root,isonode]
    k = which.min(ids)

    v1 = isonode[k]
    ds = neigDist[v1,curnode] #distances(g0, v=v1, to=curnode)
    od = order(ds)
    curnode = curnode[od]
    j=0
    if(j==0){
      nb = curnode[1]
      ns = neighbors(g5,v=nb) # get its neighbors
      ns = as.vector(ns)
      ns = ns[!(ns %in% c(isonode,v1))]
      rs = neigDist[v1, ns]  #distances(g1, v=v1, to= ns )
      od = order(rs)
      ns = ns[od]
      j = 1
    }

    if(length(ns)>0){
      ## add edge between v1 and ns[j]
      g5 = g5 + edge(v1, nb)
      # E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
      E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nb]
      g5 = g5 + edge(v1, ns[j])
      E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,ns[j]] # distances(g0,v=v1,to= ns[j])[1,1]
      #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,ns[j]]
      edges = paste(nb,"|",ns[j],sep = "")
      g5 = delete_edges(g5, edges)
      es = distances(g5,v=nb,to=ns[j])
    }else{
      print("Error")
      if(nb %in% allleaf){## add nonleaf node into nearest leaf
        #zs = distances(g5, nb, nonleaf)
        zs = cenDist[nb,nonleaf]

        nn = nonleaf[which.min(zs[1,])]
        pth =  get.shortest.paths(g5, from=nn, to=nb)
        pth = as.vector(pth$vpath[[1]])
        nb = pth[2]

        ## add edge between v1 and ns[j]
        g5 = g5 + edge(v1, nb)
        #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
        E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nb]
        g5 = g5 + edge(v1, nn)
        #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nn] # distances(g0,v=v1,to= ns[j])[1,1]
        E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nn]

        edges = paste(nb,"|",nn,sep = "")
        g5 = delete_edges(g5, edges)

      }else{
        print("Error")
        g5 = g5 + edge(v1, nb)
        #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
        E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nb]
      }
    }

    dg = degree(g5, v = nonleafroot)
    id = which(dg==0)

  }
  g5 = simplify(g5)
  # plot(g5,layout=layout.reingold.tilford(g5, circular=F), vertex.size=2,  vertex.label.cex=2)
  # pieTreeRoot(mst=g5, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  list(g5=g5,neigDist=neigDist,w=sum(E(g5)$weight))
  #gr=list(g5=g5,rootDist=rootDist)
}


#' getMSTLeafRootISOMulroot: create spanning tree with specified root and leaves
#' @param neigDist
#'     neighbor distance matrix of clusters
#' @param cenDist
#'     distance matrix of cluster centers
#' @param cluid
#'     vector, the cluster id
#' @param L2
#'    integer, default is 2, knn size
#' @param rootL
#'    integer vector, default is 1, the tree will connect the root cluster in linear mode
#'    For example, rootL = c(1,2,3), the root edges will be 1->2->3
#' @param leaves
#'    list, specify the leave groups, for each leave group, make it as a vector,
#'    For example,a=c(9,12,13,"linear"), b=c(2,5, "mst"), c=c(7,8,9,"parallel")
#' @return list
#'    list(g5=g5,neigDist=neigDist,w=sum(E(g5)$weight))
#' @export
#'
#' @examples
#'
getMSTLeafRootISOMulroot <- function(neigDist, cenDist, cluid, L2=2, rootL=1,
                                leaves = list(a=c(9,12,13,"linear"), b=17,c=c(2,6),d=3,e=c(4,5),f=14))
{
  allleaf = NULL
  for(i in 1:length(leaves)){
    cl = leaves[[i]]
    if(length(leaves[[i]]) > 1) {
      allleaf = c(allleaf, cl[1:(length(cl)-1)])
    }else{
      allleaf = c(allleaf, cl)
    }
  }
  root = rootL[1]# select the

  allleaf = as.integer(allleaf)
  nonleaf = (1:max(cluid))[-allleaf]
  nonleafroot = nonleaf[-which(nonleaf %in% rootL)]

  ###1.Build kNN graph using the distance matrix of all nodes
  M = length(unique(cluid))
  if(length(nonleaf)==1){
    ## convert into all nodes graph
    g1 <- igraph::make_empty_graph(nrow(cenDist), directed = FALSE)
  }else{
    ######## 1.1 Build kNN graph using the distance matrix of nonleaf nodes
    g0 = kNNDist(cenDist, ks=L2)
    dg = degree(g0, v=1:nrow(cenDist))
    dg
    n = components(g0, mode =  "strong" )
    mark=TRUE
    while(n$no>1){
      print("L2 is too small to build a connected graph!")
      L2 = L2 + 1
      g0 <- kNNDist(cenDist, ks=L2)
      n = components(g0, mode =  "strong" )
      print(L2)
    }
    #
    edList= cbind(get.edgelist(g0), E(g0)$weight)
    #plot(g0,layout=layout.reingold.tilford(g0, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

    for(i in 1:nrow(edList)){###1.2 Remove connection between leaf and nonleaf nodes
      if(edList[i,1] %in% allleaf  | edList[i,2] %in% allleaf){
        g0 = delete_edges(g0, edges=paste0(edList[i,1], "|", edList[i,2]))
      }
    }
    #plot(g0,layout=layout.reingold.tilford(g0, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

    ### Remove all edges among the root and keep only the linear edge
    if(length(rootL) > 1){
      for(i in 1:(length(rootL)-1)){
        nb = as.integer(neighbors(g0, v = rootL[i]))
        if(i > 1){
          id = which(nb %in% rootL[i-1])
          nb = nb[-id]
        }

        for(j in 1:length(nb)){
          g0 = delete_edges(g0, edges=paste0(nb[j], "|", rootL[i]))
        }
        ## add edge to connect root and its next root
        kn = length(E(g0)$weight) + 1
        g0 = g0 + edge(rootL[i], rootL[i+1])
        E(g0)$weight[kn] = cenDist[rootL[i], rootL[i+1]]
      }
    }
    g1 <- g0
  }
  #plot(g1,layout=layout.reingold.tilford(g1, circular=F),vertex.label.cex = 1.5,vertex.size=0.1)

  if(length(rootL) >1){
    nonleafRaw  =  nonleaf
    kn = length(rootL) -1
    sid = which(nonleaf %in% rootL[1:kn])
    nonleaf = nonleaf[-sid]
  }

  ########2. Add leaf node into the nearest nonleaf node
  ## Keep the degree of leaf node group to its parent is one
  len = length(leaves)
  for(i in 1:len){
    lv = leaves[[i]]
    if(length(lv)>1){
      typeC = lv[length(lv)]
      sv = as.integer(lv[1:(length(lv)-1)])
    }else{
      sv = as.integer(lv)
    }
    ## 1. Find nearest nonleaf node of the leaf group and add the edge
    ls = sv
    ld = sv
    for(j in 1:length(sv)){
      ds = neigDist[sv[j], nonleaf]
      ls[j] = nonleaf[which.min(ds)]
      ld[j] = min(ds)
    }
    pid = which.min(ld)
    k = length(E(g1)$weight) + 1
    g1 = g1 + edge(sv[pid], ls[pid])
    E(g1)$weight[k] = cenDist[sv[pid], ls[pid]] # use cluster distance

    ##### add edge among leave using MST
    if(length(sv)>1){
      ## 2. Add nonleaf node to its children
      if(typeC=="parallel"){
        k = length(E(g1)$weight) + 1
        for(z in 1:length(sv)){
          g1 = g1 + edge(sv[z], ls[pid])
          E(g1)$weight[k] = cenDist[sv[z], ls[pid]] # use cluster distance
          k = k + 1
        }
      }else if(typeC=="linear"){
        B <- cenDist[sv,root]
        od = order(B)
        k = length(E(g1)$weight) + 1
        for(z in 1:(length(od)-1)){
          ks = sv[od[z]]
          ke = sv[od[z+1]]
          g1 = g1 + edge(ks, ke)
          E(g1)$weight[k] = cenDist[ks, ke] # use cluster distance
          k = k + 1
        }
      }else if(typeC=="mst"){
        B <- as.matrix(neigDist[sv,sv] )
        gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
        dp_mst = mst(gp)
        se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)

        k = length(E(g1)$weight) + 1
        for(z in 1:nrow(se)){
          g1 = g1 + edge(sv[se[z,1]] ,sv[se[z,2]])
          E(g1)$weight[k] = cenDist[sv[se[z,1]], sv[se[z,2]]] # use cluster distance
          #E(g1)$weight[k] = se[z,3]
          k = k + 1
        }
      }else{
        stop("Wrong leaf types specified")
        return()
      }

    }
  }
  g1 = simplify(g1)
  edList= cbind(get.edgelist(g1), E(g1)$weight)
  print(nrow(edList))
  # plot(g1,layout=layout.reingold.tilford(g1, circular=F))
  # pieTreeRoot(mst=g1, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)

  ####3.If the graph is not fully connected, connect the nonleaf nodes to its parent
  for(i in 1:length(nonleafroot)){
    nb = as.integer(neighbors(g1, nonleafroot[i]))
    tf = nb %in% allleaf
    #if the nonleaf node is connected with only leaf node, select its nearest node
    if(length(which(tf==TRUE)) == length(tf)){
      nk =  nonleafroot[i]
      localDS = cenDist[nk, ]
      localOrd = order(localDS)
      ef = !(localOrd %in% c(nk,allleaf))
      leftID = localOrd[ef]
      k = length(E(g1)$weight) + 1
      g1 = g1 + edge(nk,leftID[1])
      E(g1)$weight[k] = cenDist[nk,leftID[1]]
    }
  }
  # plot(g1,layout=layout.reingold.tilford(g1, circular=F))

  ###### 4. Find a path from root to leaf

  g4 = igraph::make_empty_graph(M, directed = FALSE)
  paths = list() #get.all.shortest.paths(g1,from = leaves, to=root, mode = "all")
  c=1
  for(i in 1:length(leaves)){
    lv = leaves[[i]]
    if(length(lv) > 1) {
      typeC = lv[length(lv)]
      sv = as.integer(lv[1:(length(lv)-1)])
    }else{
      sv = as.integer(lv)
    }
    ##4.1 Get the path from each leaf to root
    ds = distances(g1, v=sv, to= root)
    id = which.min(ds[,1])
    p1 = get.shortest.paths(g1, from=root, to=sv[id])
    if(length(p1$vpath)==0){
      message(Sys.time(), "Error: Please set larger L2 !!!")
      return()
    }
    ev = as.vector(p1$vpath[[1]])

    if(length(ev)>1){
      ## 4.2 Add the edge into the graph from root to leaf
      h = length(ev)-1
      for(j in 1:h){
        ei <- get.edge.ids(g4, c(ev[j],ev[j+1]))
        if(ei==0){
          g4 = g4 + edge(ev[j],ev[j+1])
          E(g4)$weight[c] = distances(g1, v=ev[j], to= ev[j+1])[1,1]
          c = c + 1
        }
      }

      if(length(sv) > 1){#4.3 Add edge among the leaf
        if(typeC=="parallel"){# add edge between parent and leaf
          pid = ev[length(ev)-1]
          for(z in 1:length(sv)){
            if(sv[z] != pid){
              g4 = g4 + edge(pid, sv[z] )
              E(g4)$weight[c] = cenDist[pid,sv[z] ]
              c = c + 1
            }
          }
        }else if(typeC=="linear"){
          B <- cenDist[sv,root]
          od = order(B)
          for(z in 1:(length(od)-1)){
            ks = sv[od[z]]
            ke = sv[od[z+1]]
            g4 = g4 + edge(ks, ke)
            E(g4)$weight[c] = cenDist[ks, ke] # use cluster distance
            c = c + 1
          }
        }else if(typeC=="mst"){
          B <- as.matrix(neigDist[sv,sv] )
          gp <- graph.adjacency(B, mode = "undirected", weighted = TRUE)
          dp_mst = mst(gp)
          se = cbind(get.edgelist(dp_mst), E(dp_mst)$weight)
          for(z in 1:nrow(se)){
            g4 = g4 + edge(sv[se[z,1]], sv[se[z,2]])
            E(g4)$weight[c] = cenDist[sv[se[z,1]], sv[se[z,2]]]  #se[z,3]
            c = c + 1
          }
        }else{
          stop("Wrong leaf types specified")
          return()
        }
      }
    }
  }
  g4 = simplify(g4)
  # plot(g4,layout=layout.reingold.tilford(g4, circular=F))
  # pieTreeRoot(mst=g4, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  ###### 4. Merge isolated nonleaf node into the nearest edges
  nonleafroot = nonleaf#[-which(nonleaf==root)]
  dg = degree(g4, v = nonleafroot)
  id = which(dg==0)

  g5 = g4
  while (length(id)>0) {
    isonode = nonleafroot[id]
    if(length(isonode)==1){

    }
    curnode = nonleaf[ !(nonleaf %in% isonode) ]
    #rds = cenDist[root,curnode] #distances(g0, v=root, to=curnode)
    #ids = cenDist[root,isonode] #distances(g0, v=root, to=isonode)
    ids =  cenDist[root,isonode]
    k = which.min(ids)

    v1 = isonode[k]
    ds = neigDist[v1,curnode] #distances(g0, v=v1, to=curnode)
    od = order(ds)
    curnode = curnode[od]
    j=0
    if(j==0){
      nb = curnode[1]
      ns = neighbors(g5,v=nb) # get its neighbors
      ns = as.vector(ns)
      ns = ns[!(ns %in% c(isonode,v1))]
      rs = neigDist[v1, ns]  #distances(g1, v=v1, to= ns )
      od = order(rs)
      ns = ns[od]
      j = 1
    }

    if(length(ns)>0){
      ## add edge between v1 and ns[j]
      g5 = g5 + edge(v1, nb)
      # E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
      E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nb]
      g5 = g5 + edge(v1, ns[j])
      E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,ns[j]] # distances(g0,v=v1,to= ns[j])[1,1]
      #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,ns[j]]
      edges = paste(nb,"|",ns[j],sep = "")
      g5 = delete_edges(g5, edges)
      es = distances(g5,v=nb,to=ns[j])
    }else{
      print("Error")
      if(nb %in% allleaf){## add nonleaf node into nearest leaf
        #zs = distances(g5, nb, nonleaf)
        zs = cenDist[nb,nonleaf]

        nn = nonleaf[which.min(zs[1,])]
        pth =  get.shortest.paths(g5, from=nn, to=nb)
        pth = as.vector(pth$vpath[[1]])
        nb = pth[2]

        ## add edge between v1 and ns[j]
        g5 = g5 + edge(v1, nb)
        #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
        E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nb]
        g5 = g5 + edge(v1, nn)
        #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nn] # distances(g0,v=v1,to= ns[j])[1,1]
        E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nn]

        edges = paste(nb,"|",nn,sep = "")
        g5 = delete_edges(g5, edges)

      }else{
        print("Error")
        g5 = g5 + edge(v1, nb)
        #E(g5)$weight[length(E(g5)$weight)] = neigDist[v1,nb] # distances(g0,v=v1,to= nb)[1,1]
        E(g5)$weight[length(E(g5)$weight)] = cenDist[v1,nb]
      }
    }

    dg = degree(g5, v = nonleafroot)
    id = which(dg==0)

  }
  g5 = simplify(g5)
  # plot(g5,layout=layout.reingold.tilford(g5, circular=F), vertex.size=2,  vertex.label.cex=2)
  # pieTreeRoot(mst=g5, root=rt,ctyM, cellN, vs=0.5, ew=1, cex=1.3, seed=1)
  list(g5=g5,neigDist=neigDist,w=sum(E(g5)$weight))
  #gr=list(g5=g5,rootDist=rootDist)
}



#' pieTreeRoot2: Visaulize the trajectory in pie tree
#'
#' @param mst
#'     igraph, the trajectory
#' @param root
#'     integer,  the root cluster id
#' @param ctyM
#'     cell type matrix
#' @param cellN
#'     vector, the names of clusters
#' @param logSize
#'     Default is FALSE, to do log on the cell number to reduce the cirdius of pie
#' @param vs
#'     Default 1
#' @param ew
#'     Default 1 the label size
#' @param seed
#'     1
#' @param cex
#'     2
#' @param lbd
#'     1
#' @param lgs
#'     2
#' @param title
#'     string
#' @param colorname
#'     colorname: rownames(brewer.pal.info): "BrBG"     "PiYG"     "PRGn"     "PuOr"     "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"   "Spectral" "Accent"   "Dark2"
#'     "Paired"   "Pastel1"  "Pastel2"  "Set1"     "Set2"     "Set3"     "Blues"    "BuGn"     "BuPu"     "GnBu"     "Greens"
#'     "Greys"    "Oranges"  "OrRd"     "PuBu"     "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"     "YlGn"     "YlGnBu"
#'    "YlOrBr"   "YlOrRd"
#' @param showLab
#'     FALSE, show labels
#' @return
#'
#' @export
#'
#' @examples
#'
pieTreeRoot2 <- function(mst, root=1,ctyM, cellN, logSize=FALSE, vs=1, ew=1, seed=1, cex=2, lbd=1,
                         lgs=2,title="Cells",  colorname="Spectral", showLab=FALSE) {

  dmst = as.directed(mst, mode = "mutual" )
  ##delete edges direct to root
  nb = neighbors(dmst, v=root, mode =  "all")
  len = length(nb)
  nn = as.integer(as.character(nb))
  nn = unique(nn)
  dnb = bfs(dmst, root=V(dmst)[root], neimode="out")
  od = as.vector(dnb$order)
  for(i in 2:length(od)){
    ch = od[i]
    for(j in 1:(i-1)){
      fa = od[j]
      nb = as.vector(neighbors(dmst, v=V(dmst)[ch], mode =  "all"))
      if(fa %in% nb){
        if(is.null(V(mst)[[ch]]$name)){
          dmst = delete_edges(dmst, edges= paste( V(mst)[[ch]],"|", V(mst)[[fa]],sep="") )
        }else{
          dmst = delete_edges(dmst, edges= paste( V(mst)[[ch]]$name,"|", V(mst)[[fa]]$name,sep="") )
        }
      }
    }
  }
  # plot(dmst,layout=layout.reingold.tilford(dmst, circular=F))
  graph2 = dmst
  graph2$name =NULL
  graph2$size =NULL
  graph2$shape=NULL
  graph2$pie=NULL
  graph2$pie.color=NULL
  graph2$label=NULL
  graph2$weight=NULL
  graph2$color = NULL
  graph2$width = NULL

  rows = dim(ctyM)[1]
  cols = dim(ctyM)[2]

  labs = rep(" ", rows)
  vsize = rep(0,rows)
  s = NULL
  for(i in 1:rows){
    m=0
    for(j in 1:cols){
      s = paste(s,ctyM[i,j],sep = "\n")
      m = m + ctyM[i,j]
    }
    labs[i] =s
    vsize[i] = m
  }

  values =  lapply(1:rows, function(x) ctyM[x, ])
  if(logSize){
    vsize <- sqrt(vsize)*vs
  }else{
    vsize <- sqrt(vsize)*vs #
  }

  coid = which(rownames(brewer.pal.info)==colorname )
  ac = colorRampPalette(brewer.pal(brewer.pal.info[coid,1], colorname))(ncol(ctyM))
  lbs = list()
  for(i in 1:nrow(ctyM)){
    id = which(ctyM[i,]>10)
    prop = ctyM[i, id] / sum(ctyM[i,])
    str = paste(cellN[id],round(prop, 2), sep = '')
    od = order(prop, decreasing = TRUE)
    lbs[[i]] = paste(i, str[od], collapse = '\n')
  }
  color1 = list(ac)

  ########clear plot
  #dev.off(dev.list()["RStudioGD"])

  too_many_pops <-  1
  pops_correction <-  lgs
  yintersperse <- 0.5

  wg = ew
  set.seed(seed)
  par(mar=c(0, 0, 0, 14.1), bg=NA)
  if(showLab){
    plot(graph2, vertex.shape="pie",
         vertex.pie=values,
         vertex.pie.color=color1,
         vertex.size=vsize,vertex.label.cex = cex,vertex.label.dist= lbd,
         edge.width= wg, vertex.label= lbs,
         edge.arrow.width = 0.5, edge.arrow.size=0.5,
         layout=layout.reingold.tilford(graph2, circular=F)) #, vertex.label=NA
  }else{
    plot(graph2, vertex.shape="pie",
         vertex.pie=values,
         vertex.pie.color=color1,
         vertex.size=vsize,vertex.label.cex = cex,vertex.label.dist= lbd,
         edge.arrow.width = 0.5, edge.arrow.size=0.5,
         edge.width= wg, layout=layout.reingold.tilford(graph2, circular=F)) #, vertex.label=NA
  }
  legend("right", bty = "n", cex = 1.2^pops_correction, inset=c(-1.3,0),
         legend = cellN, fill = ac, bg="transparent", xpd=TRUE,
         border = NULL, ncol = too_many_pops, x.intersp = 0.2,
         y.intersp = yintersperse)
}


#' runShinyVis: run the shiny app
#' @export
runShinyVis <- function() {
  appDir <- system.file("shinyApp", "visSCgraph", package = "LISA2")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing 'LISA2'.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

