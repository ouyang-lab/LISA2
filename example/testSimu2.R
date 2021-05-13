library(FFT)
library(igraph)
load("data/simulation_res.RData")

###### New one
pcd = pcaL1(data=daC, pvalue=0.05, pcaSel = 1,  sel=2, center=TRUE, scale=FALSE, ith=1e-2)
nc = getPCAN(pcd$pca_out,pvalue=0.05, sel=2, ith=1e-2)
Z = pcd$pca_out$x[,1:nc]
dim(Z)
#Z = pcd$pca_out$x[,1:20]

library(umap)
set.seed(1)
emb2 <-  umap::umap(Z)
plotAlls(X=emb2$layout,draw=1, t=NULL, cluid= branch,
         type=2,size =6,opacity=1, titlename=paste0("UMAP: PC number = ", ncol(Z)))
uY = emb2$layout


tclu3 = proc.time()
ks = seq(4,20, by=2)
cluLeiden = list()
setLeignPath()
for(i in 1:length(ks)){
  cluLeiden[[i]] = cmcluster(Z, ksize = ks[i], cmmethod="leiden")
}
tclu4 = proc.time() - tclu3


library(fpc)
t1 = proc.time()
dsM = distances(cluLeiden[[k]]$knng$g)
proc.time() - t1
dim(dsM)
etc = cluster.stats(dsM, cluLeiden[[k]]$cluid)
msM = etc$ave.between.matrix
msM = etc$separation.matrix


cm = cmcluster(Z, ksize = 11, cmmethod="leiden")
cluid = cm$cluid
plotAlls(X=uY,draw=1, t=NULL, cluid=as.character(cluid), type=2,size =6,opacity=1,
         titlename=paste("UMAP-leiden,k=", 11, sep=""))

k=5
cluid = cluLeiden[[k]]$cluid
table(cluid)
plotAlls(X=uY,draw=1, t=NULL, cluid=as.character(cluid), type=2,size =6,opacity=1,
         titlename=paste("UMAP-leiden,k=", ks[k], sep=""))

#### 3. Do L-ISOMAP
#centid = getCluCenterID(Z, cluid)
peakid = searchPeaks(gg=cluLeiden[[k]]$knng$g, cluid)
X  <- autoLISOMAP(Z, peaks=peakid, ndim = 3, L1=NULL, L2=NULL, N=1)
PZ = rbind(uY, uY[peakid,])
pc = c(cluid, rep("peaks", length(peakid)))
plotAlls(X=PZ,draw=1, t=NULL, cluid=as.character(pc), type=2,size =6,opacity=1,
         titlename=paste("Peaks in UMAP,k=", ks[k], sep=""))

#### 4. Visualization
plotAlls(X=X,draw=3, t=NULL, cluid=as.character(branch), type=2,size =6,opacity=1, titlename="ISOMAP")
plotAlls(X=X,draw=3, t=NULL, cluid=as.character(cluid), type=2,size =5,opacity=1, titlename="ISOMAP")

#### 5. Build trajectory

#### 5.1 User specified root and leaves to build spanning tree

##1. Use neighbor distance of clusters in graph built by PCA, but not good
msM = getCluNeighbors(Z, k=ks[k], cluid=cluid)
##2.
dsM = distGraphClu(Z, ksize=ks[k], cluid=cluid)


rt = 3 # clu
g1 = getMSTLeafRootISO_1(neigDist = msM, cenDist=dsM, cluid=cluid, L2=3, root=rt,
                       leaves=list(a= 9 , b=10, c=5, d= 8, e=1, f=6,g= 11 ))


rt = 2 # clu
g1 = getMSTLeafRootISO_3(neigDist = msM, cenDist=dsM, cluid=cluid, L2=3, root=rt,
                         leaves=list(a=c(5,7,"parallel"), b=c(8,9,"parallel"),
                                     c=c(16,12,"linear"), d=10, f=4))

clucells = getCluCellTh(cluid=cluid, branch, th=0.2, subN=10)
g2 <- set.vertex.attribute(g1$g5, "name", value=clucells)
ctyM = getCellSM(cluid, branch)
cellN =  as.character(names(table(branch)))

g3 <- set.vertex.attribute(g1$g5, "name", value= paste("Y_",1:max(cluid),sep=""))
pieTreeRoot(mst=g3, root=rt, ctyM, cellN, vs=1, ew=1, cex=2, seed=1,colorname="Paired") # Plot the cluster tree

#### Plot 3D trajectory
plot3DTreeCells(X[,1:3], peaks =NULL, cluid= as.character(cluid), celltype = as.character(branch),
                colorname="Paired", mst=g1$g5, pointsize = 6,labelsize=15, lineW=10)
plot3DTreeCells(X[,1:3], peaks =NULL, cluid= as.character(cluid),colorname="Paired",
                celltype = as.character(cluid),mst=g1$g5, pointsize = 6,labelsize=15, lineW=10)

####   Plot the tree
Y <- getCenters(X, cluid)
cluidS = as.character(cluid)
anno =  cbind(cluidS, as.character(branch))
colnames(anno) = c("cluster", "celltype")
anno = as.data.frame(anno)

g3 <- set.vertex.attribute(g1$g5, "name", value= paste("Y_",1:max(cluid),sep=""))
#g3 <- set.vertex.attribute(g1, "name", value= paste("Y_",1:max(cluid),sep=""))
plotCellTree(X=X, uY, cluid=cluid, mst=g3, root=rt, colorby="celltype", anno=anno, height=0.2,
             width=0.3, pointSize=2, edgeW=0.5, textSize=6, legendTextS=20, legendColS=10, hv=0.5,
             vv=-1.2, colorName = "Paired")

save(a=pcd, b=Z, c=uY, d=cluLeiden, e=X, f=g1, g=branch, h=daC, i=cluid,  file = "data/simulation_res.RData")


##################################################################################
### PTA analysis
######################################################################################
load("data/simulation_res.RData")

library(FFT)
library(PTA)
library(doParallel)

lpath1 = c(2,15,11,9)
lpath2 = c(2,15,11,8)
lpath3 = c(2,14,13,6,4)
lpath4 = c(2,14,13,6,3,7)
lpath5 = c(2,14,13,6,3,5)
lpath6 = c(2,14,13,1,16,12)
lpath7 = c(2,14,13,1,10)

cluid = cluLeiden[[5]]$cluid
fpath = "example/simu_PTA_lineage_1.RData"
usePTA(sd=daC, cluid=cluid, path=lpath1, pt = res$pt, fname=fpath,
       it = 50, intv = 400, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)


cluid = cluLeiden[[5]]$cluid
fpath = "data/simu_PTA_lineage_1.RData"
usePTA(sd=daC, cluid=cluid, path=lpath1, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)

fpath = "data/simu_PTA_lineage_2.RData"
usePTA(sd=daC, cluid=cluid, path=lpath2, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)

fpath = "data/simu_PTA_lineage_3.RData"
usePTA(sd=daC, cluid=cluid, path=lpath3, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=4)

fpath = "data/simu_PTA_lineage_4.RData"
usePTA(sd=daC, cluid=cluid, path=lpath4, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)


fpath = "data/simu_PTA_lineage_5.RData"
usePTA(sd=daC, cluid=cluid, path=lpath5, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)


fpath = "data/simu_PTA_lineage_6.RData"
usePTA(sd=daC, cluid=cluid, path=lpath6, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)


fpath = "data/simu_PTA_lineage_7.RData"
usePTA(sd=daC, cluid=cluid, path=lpath7, pt = res$pt, fname=fpath,
       it = 25, intv = 200, feature.flag = TRUE,
       err_flag = 1, scale=FALSE, nCores=1)
