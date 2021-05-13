
#' plotDens: plot 3D density plot of tZ
#'
#' @param tZ
#'     a matrix of tSNE
#' @param pdf
#'     An S4 object of kepdf-class with slots, see details in kepdf {pdfCluster}
#' @param phi
#'    phi, to adjust the degree of the up and down view of the figure in [0,360] range.
#' @param theta
#'    theta, to adjust the degree of the left and right view of the figure in [0,360] range.
#' @return title
#'    The name of figure
#' @export
#'
#' @examples
#'
plotDens <- function(tZ, pdf,phi=30, theta=100, title="data"){
  xm = max(tZ[,1]) - min(tZ[,1])
  ym =  max(tZ[,2]) - min(tZ[,2])
  xg = ceiling(200*xm/(xm+ym))
  yg = ceiling(200*ym/(xm+ym))
  plot(pdf, main= paste(title, ": kenerl density estimate of data"),  ticktype="simple", col="green",
       n.grid=c(xg,yg), method="perspective", phi=phi, theta=theta, cex.sub = 1.3, font.sub = 2,
       axes=FALSE,ltheta = 150,
       xlim=c(min(tZ[,1])-ym*0.01, max(tZ[,1])+ym*0.01),
       ylim=c(min(tZ[,2])-ym*0.01, max(tZ[,2])+ym*0.01))
}


#' plotISOMST: plot 2D cell trajectory on ISOMAP or PCA space
#'
#' @param X
#'        a matrix of ISOMAP or PCA
#' @param type
#'      "ISOMAP" or "PCA". It just changes the axis labels.
#' @param color_by "State", "Time", "Marker"
#' @param memb
#'     The clusteirng id
#' @param discrete
#'     To show whether the value of memb is discreate or continuous.
#' @param markername
#'     This the marker name when setting colory_by="Marker"
#' @param x
#'     The dimension of ISOMAP/PCA to visualize, default 1
#' @param y
#'     The dimension of ISOMAP/PCA to visualize, default 2
#' @param cenTree
#'    The minimum spanning tree
#' @param show_tree
#'    logical, (default TRUE), whether to draw the tree
#' @param show_cell_names
#'    logical, (default FALSE),Show cell name for each cell
#' @param cell_name_size
#'    Adjust the cell name size
#' @param pointsize
#'    Adjust the points size
#' @param lengendsize
#'    Adjust the lengend size
#' @param labelsize
#'    Adjust the label size for the verticle in the tree
#' @return
#'
#' @export
#'
#' @examples
#'
plotISOMST <- function (X, type= c("ISOMAP", "PCA"), color_by = c("State", "Time", "Marker"),
                        memb, discrete=TRUE, markername ="",
                        x = 1, y = 2, Y,cenTree, show_tree = TRUE,
                        show_cell_names = FALSE, cell_name_size = 3,
                        pointsize=2, lengendsize = 20, labelsize = 7 )
{
  names(memb) <- rownames(X)
  if(color_by == "State"){
    lib_info_with_pseudo <- data.frame(State = memb,
                                       sample_name = rownames(X))
    lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
    markerexpr <- memb
  }else if(color_by == "Time"){
    lib_info_with_pseudo <- data.frame(Time = memb,
                                       sample_name = rownames(X))
    if(discrete){
      lib_info_with_pseudo$Time <- factor(lib_info_with_pseudo$Time)
    }else{
      lib_info_with_pseudo$Time <-  lib_info_with_pseudo$Time
    }
    markername <- "Time"
    markerexpr <- memb
  }else if(color_by == "Marker"){
    lib_info_with_pseudo <- data.frame(Marker = memb,
                                       sample_name = rownames(X))
    lib_info_with_pseudo$Marker <-  lib_info_with_pseudo$Marker
    markerexpr <- memb
  }


  S_matrix <- X
  pca_space_df <- data.frame(S_matrix[, c(x, y)])
  colnames(pca_space_df) <- colnames(S_matrix)[c(x,y)]
  dimnames <- colnames(S_matrix)
  pca_space_df$sample_name <- row.names(pca_space_df)

  edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name",
                   by.y = "sample_name")
  edge_df$markerexpr <- markerexpr[edge_df$sample_name]
  if(type=="ISOMAP"){
    g <- ggplot(data = edge_df, aes(x = iso1, y = iso2))
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE,
                        size = pointsize)
  }else if(type=="PCA"){
    g <- ggplot(data = edge_df, aes(x = PC1, y = PC2))
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE,
                        size = pointsize)
  }

  if (show_cell_names & discrete) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }

  if (show_tree) {
    clucenter <- Y[, c(x, y)]
    clulines <- NULL

    edMat = get.edgelist(cenTree)
    for (i in 1:nrow(edMat)) {
      st = as.integer( strsplit(edMat[i,1], split = "_" )[[1]][2]  )
      ed = as.integer( strsplit(edMat[i,2], split = "_" )[[1]][2] )
      clulines <- rbind(clulines, c(clucenter[st,], clucenter[ed, ]))
    }
    clulines <- data.frame(x = clulines[, 1], xend = clulines[, 3], y = clulines[, 2], yend = clulines[, 4])
    g <- g + geom_segment(aes_string(x = "x", xend = "xend",
                                     y = "y", yend = "yend", size = NULL), data = clulines,
                          size = 1)
    clucenter <- data.frame(x = clucenter[, 1], y = clucenter[, 2], id = 1:nrow(clucenter))
    g <- g + geom_text(aes_string(label = "id", x = "x",
                                  y = "y", size = NULL), data = clucenter, size = labelsize)
  }
  if(discrete){
    g <- g + guides(colour = guide_legend(override.aes = list(size = 8))) +
      theme(panel.border = element_blank(), axis.line = element_line()) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
      theme(legend.position = "top", legend.key.size = unit(0.3, "in"),
            legend.text = element_text(size = lengendsize), legend.title = element_text(size = lengendsize)) +
      theme(legend.key = element_blank()) + theme(panel.background = element_rect(fill = "white")) +
      theme(axis.text.x = element_text(size = 17, color = "black"),
            axis.text.y = element_text(size = 17, color = "black"),
            axis.title.x = element_text(size = 20, vjust = -1),
            axis.title.y = element_text(size = 20, vjust = 1),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
  }else{
    g <- g + scale_colour_continuous(name  = markername)
    theme(panel.border = element_blank(), axis.line = element_line()) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
      theme(legend.position = "top", legend.key.size = unit(0.3, "in"),
            legend.text = element_text(size = lengendsize), legend.title = element_text(size = lengendsize)) +
      theme(legend.key = element_blank()) + theme(panel.background = element_rect(fill = "white")) +
      theme(axis.text.x = element_text(size = 17, color = "black"),
            axis.text.y = element_text(size = 17, color = "black"),
            axis.title.x = element_text(size = 20, vjust = -1),
            axis.title.y = element_text(size = 20, vjust = 1),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
  }

  g
}


#' drawMarkerValue: plot marker genes expression among all clusters
#'
#' @param sd
#'        a matrix of gene expression: row: marker genes, column: cells
#' @param memb
#'     The clusteirng id
#' @param min.exp
#'     Default 0.05. Minimum proportion of expressing cells (0-1) to be shown on the plot
#' @param mean.expressing.only
#'     (Logical) Should mean expression value exclude cells with no expression
#' @param title
#'     The title of the plot
#' @return
#'
#' @export
#'
#' @examples
#'
drawMarkerValue <- function(sd, memb,  min.exp = 0.05,
                            mean.expressing.only = F, title = "") {
  genes = rownames(sd)
  cluster.ids <- memb
  data.use <- as.data.frame(as.matrix(t(sd)))

  if (mean.expressing.only)
    mean.fun <- "mean.of.logs.pos"
  else mean.fun <- "mean.of.logs"
  mean.expression <- aggregate(data.use, by = list(cluster.ids),
                               FUN = mean.fun)
  percent.expressed <- aggregate(data.use, by = list(cluster.ids),
                                 FUN = "prop.exp")
  mean.melt <- reshape2::melt(data = mean.expression, id.vars = c("Group.1"))
  names(mean.melt) <- c("Cluster", "Gene", "Mean")
  expr.melt <- reshape2::melt(data = percent.expressed, id.vars = c("Group.1"))
  names(expr.melt) <- c("Cluster", "Gene", "Prop.Exp")
  gg.data <- merge(mean.melt, expr.melt)
  gg.data <- gg.data[which(gg.data$Prop.Exp >= min.exp), ]
  the.plot <- ggplot(data = gg.data, aes(x = Cluster, y = Gene,
                                         color = Mean, size = Prop.Exp)) + geom_point() +
    scale_color_gradientn(colours = defaultURDContinuousColors()) +
    ggtitle(title) + scale_y_discrete(limits = rev(genes))

  the.plot <- the.plot + scale_x_discrete(limits = cluster.ids)
  return(the.plot)
}


#' plot3DTree: Visaulize trees on the 3D ISOMAP
#'
#' @param X
#'        3D ISOMAP:row samples, col:3D
#' @param cluid
#'     The clusteirng id
#' @param pointsize
#'     Default 2 for point size
#' @param labelsize
#'     Default 20. For enlarge the label size in the legend
#' @param lineW
#'     Default 5. line width in the MST
#' @param colorname
#'     colorname: rownames(brewer.pal.info): "BrBG"     "PiYG"     "PRGn"     "PuOr"     "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"   "Spectral" "Accent"   "Dark2"
#'     "Paired"   "Pastel1"  "Pastel2"  "Set1"     "Set2"     "Set3"     "Blues"    "BuGn"     "BuPu"     "GnBu"     "Greens"
#'     "Greys"    "Oranges"  "OrRd"     "PuBu"     "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"     "YlGn"     "YlGnBu"
#'    "YlOrBr"   "YlOrRd"
#' @param titles
#'     The title of the plot
#' @return
#'
#' @export
#'
#' @examples
#'
plot3DTree <- function(X, cluid,
                       pointsize=2, labelsize=20, lineW=5, colorname="Spectral", titles="3D ISOMAP"){

  library(RColorBrewer)
  cluidL <- c(cluid, rep("cluster center", length(unique(cluid))))

  Y <- getCenters(X, cluid) #   cluidL
  mst <- getCenMST(Y)
  rX <- rbind(X, Y)

  if(dim(X)[2]==2){
    da = data.frame(x = rX[,1], y = rX[,2], z = rep(0, nrow(rX)), clu=cluidL)
  }else if(dim(X)[2]==3){
    da = data.frame(x = rX[,1], y = rX[,2], z = rX[,3], clu=cluidL)
  }

  coid = which(rownames(brewer.pal.info)==colorname)
  colors = colorRampPalette(brewer.pal(brewer.pal.info[coid,1], colorname))(length(unique(cluidL)))

  p <- plot_ly() %>%
    add_trace(data = da, x= ~x, y= ~y, z= ~z,
              color = ~as.factor(cluidL),
              type="scatter3d", colors = colors,
              marker = list(opacity=1, size=2),showlegend = TRUE)

  p <- plotly_build(p)
  markersize <- pointsize
  markerlegendsize <- labelsize
  for(i in seq(1, length(sort(unique(cluidL) )))) {
    length.group <- nrow(da[which(cluidL == sort(unique(cluidL))[i]), ])
    p$x$data[[i]]$marker$size <- c(rep(markersize,length.group), rep(c(-markersize+2*markerlegendsize), length.group))
  }

  SY <- Y[1:2,]
  edMat = get.edgelist(mst)
  for (i in 1:nrow(edMat)) {
    st = as.integer(strsplit(edMat[i,1],'_')[[1]][2])
    ed = as.integer(strsplit(edMat[i,2],'_')[[1]][2])

    SY[1,] = Y[st,]
    SY[2,] = Y[ed,]
    p <- p %>%
      add_trace(x =SY[,1], y= SY[,2], z=SY[,3], mode = "lines",  line = list(color = "black", lineW = 5),
                name = paste("e",st, ed, sep = "_"), inherit=FALSE )
  }
  p %>%
    layout(title = paste( titles, sep=""), legend =list(font = list(size=labelsize)))

}


#' plotAlls: Visaulize 3D ISOMAP with different annotation
#'
#' @param X
#'     matrix
#' @param draw
#'     integer, 1:draw 2D in 3D space; 2: draw 2D 3: draw 3D
#' @param t
#'     Default for pseudotime annotation
#' @param cluid
#'     cluster id or other  annotation
#' @param type
#'     integer, 1: using t annotation 2: use cluid annotation 3. use cluid annotation and mcs colors
#' @param size
#'     Default 4. the cell size
#' @param opacity
#'     Default 1
#' @param labelsize
#'     Default 20. the label size
#' @param titlename
#'     string, the title
#' @param colorname
#'     colorname: rownames(brewer.pal.info): "BrBG"     "PiYG"     "PRGn"     "PuOr"     "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"   "Spectral" "Accent"   "Dark2"
#'     "Paired"   "Pastel1"  "Pastel2"  "Set1"     "Set2"     "Set3"     "Blues"    "BuGn"     "BuPu"     "GnBu"     "Greens"
#'     "Greys"    "Oranges"  "OrRd"     "PuBu"     "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"     "YlGn"     "YlGnBu"
#'    "YlOrBr"   "YlOrRd"
#' @param mcs
#'     The color vector specified by users
#' @return
#'
#' @export
#'
#' @examples
#'
plotAlls <- function(X=NULL,draw=2, t=NULL, cluid=NULL, type=1,size =4,opacity=1,labelsize=20,
                     titlename="ISOMAP", colname= rownames( brewer.pal.info)[12], mcs=NULL){

  if(draw==1){
    dc = data.frame(x = X[,1], y = X[,2], z = rep(0, nrow(X)))
  }else if(draw==3){
    dc = data.frame(x = X[,1], y = X[,2], z = X[,3])
  }else if(draw==2){
    dc = data.frame(x = X[,1], y = X[,2])
  }

  if(!is.null(cluid)){
    ud = unique(cluid)
    cluid = factor(cluid, levels = ud[order(ud)])
    scluid = factor(cluid, levels = as.character(ud[order(ud)]))
  }

  if(ncol(dc)==2){
    if(is.null(t) & is.null(cluid)){
      p <- plot_ly(dc, x=~x, y=~y,
                   marker = list(size = size,opacity=opacity))
    }else{
      if(type==1){
        p <- plot_ly(dc, x=~x, y=~y,
                     color =  t ,
                     marker = list(size = size,opacity=opacity))
      }else if(type==2){

        N= length(unique(cluid))
        if(N > brewer.pal.info[colname,1]){
          cs = colorRampPalette(brewer.pal(brewer.pal.info[colname,1],colname))(N)
        }else{
          cs = brewer.pal(N,colname)
        }
        cu = unique(cluid)
        cu = cu[order(cu)]
        css = sapply(cluid, function(x){ cs[which(cu==x)] })
        rid = order(as.character(cu))
        nod = order(rid)
        #cs = cs[nod]

        p <-  plot_ly(dc, x=~x, y=~y,   mode="markers",
                      color = scluid, colors=cs,
                      marker = list(size = size,opacity=opacity ))
      }else if(type==3){
        p <-  plot_ly(dc, x=~x, y=~y,   mode="markers",
                      color = scluid, colors=mcs,
                      marker = list(size = size,opacity=opacity ))
      }
    }
  }else{
    if(is.null(t) & is.null(cluid)){
      plot_ly(dc, x=~x, y=~y,z=~z,  mode="markers", type="scatter3d",
              marker = list(size = size,opacity=opacity))
    }else{
      if(type==1){
        p <- plot_ly(dc, x=~x, y=~y,  z=~z,   mode="markers",
                     color =  t , type="scatter3d",
                     marker = list(size = size,opacity=opacity))
      }else if(type==2){
        N= length(unique(cluid))
        if(N > brewer.pal.info[colname,1]){
          cs = colorRampPalette(brewer.pal(brewer.pal.info[colname,1],colname))(N)
        }else{
          cs = brewer.pal(N,colname)
        }
        cu = unique(cluid)
        cu = cu[order(cu)]
        css = sapply(cluid, function(x){ cs[which(cu==x)] })
        rid = order(as.character(cu))
        nod = order(rid)
        #cs = cs[nod]

        p <-  plot_ly(dc, x=~x, y=~y, z=~z,  mode="markers",
                      color = scluid,
                      colors= cs, type="scatter3d",
                      marker = list(size = size,opacity=opacity))
      }else if(type==3){
        p <-  plot_ly(dc, x=~x, y=~y,z=~z, mode="markers", type="scatter3d",
                      color = scluid, colors=mcs,
                      marker = list(size = size,opacity=opacity ))
      }
    }
  }

  if(type>1){
    p <- plotly_build(p)
    markersize <- size
    markerlegendsize <- labelsize
    kid = sort(unique(cluid) )
    len = length(p$x$data)
    for(i in 1:len) {
      length.group <- nrow(dc[which(cluid == kid[i]), ])
      p$x$data[[i]]$marker$size <- c(rep(markersize,length.group), rep(c(-markersize+2*markerlegendsize), length.group))
    }
  }

  p %>%
    layout(title = paste( titlename, sep=""), legend =list(font = list(size=labelsize)))
}



#' pieTreeRoot: Visaulize the trajectory in pie tree
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
#' @param legendFigWidth
#'     Default is c(2,4), 2 is the width of legend, 4 is the width of figure.
#'     User can change them to adjust the layout of legend and figure
#' @param vs
#'     Default 100
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
pieTreeRoot <- function(mst, root=1,ctyM, cellN, logSize=FALSE, legendFigWidth = c(2,4),
                        vs=100, ew=1, seed=1, cex=2, lbd=1,
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

  #coid = which(rownames(brewer.pal.info)==colorname )
  #ac = colorRampPalette(brewer.pal(brewer.pal.info[coid,1], colorname))(ncol(ctyM))
  N = ncol(ctyM)
  if(N > brewer.pal.info[colorname,1]){
    ac = colorRampPalette(brewer.pal(brewer.pal.info[colorname,1],colorname))(N)
  }else{
    ac = brewer.pal(N,colorname)
  }

  lbs = list()
  for(i in 1:nrow(ctyM)){
    id = which(ctyM[i,]>10)
    prop = ctyM[i, id] / sum(ctyM[i,])
    str = paste(cellN[id],round(prop, 2), sep = '')
    od = order(prop, decreasing = TRUE)
    lbs[[i]] = paste(i, str[od], collapse = '\n')
  }
  color1 = list(ac)

  ########
  # Plotting parameters.
  def.par <- par(no.readonly = TRUE)
  # Return top level plot to defaults.
  on.exit({
    graphics::layout(matrix(1, ncol = 1, byrow = TRUE))
    graphics::par(def.par)
  })

  layout(matrix(c(1, 2, 1, 3), ncol = 2, byrow = TRUE),
         widths = legendFigWidth, heights = c(4.5, 0.5))

  par(mar = c(0, 0, 1, 0) + 0.5, bg=NA,xpd=TRUE)
  too_many_pops <-  1
  pops_correction <-  lgs
  yintersperse <- 0.62

  graphics::plot(c(0, 2), c(0, 1), type = "n", axes = F, xlab = "", ylab = "", main = title, bg="transparent")
  graphics::legend("topright", bty = "n", cex = 1.2^pops_correction,
                   legend = cellN, fill = ac, bg="transparent",
                   border = NULL, ncol = too_many_pops, x.intersp = 0.45,
                   y.intersp = yintersperse)

  wg = ew
  set.seed(seed)
  par(bg=NA)
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

}



#' plot3DTreeCells: Visaulize the trajectory in 3D ISOMAP
#'
#' @param X
#'     matrix, the 3d ISOMAP
#' @param peaks
#'     Default is NULL, use the cluster centers as the vertexes.
#' @param cluid
#'     vector, the cluster id
#' @param celltype
#'     vector, the cell type annotation
#' @param mst
#'     igraph, the trajectory
#' @param opcacity
#'     Default 0.7
#' @param pointsize
#'     Default 2
#' @param labelsize
#'     Default 20
#' @param lineW
#'     5
#' @param colorname
#'     colorname: rownames(brewer.pal.info): "BrBG"     "PiYG"     "PRGn"     "PuOr"     "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"   "Spectral" "Accent"   "Dark2"
#'     "Paired"   "Pastel1"  "Pastel2"  "Set1"     "Set2"     "Set3"     "Blues"    "BuGn"     "BuPu"     "GnBu"     "Greens"
#'     "Greys"    "Oranges"  "OrRd"     "PuBu"     "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"     "YlGn"     "YlGnBu"
#'    "YlOrBr"   "YlOrRd"
#' @param titles
#'     FALSE, show labels
#' @return
#'
#' @export
#'
#' @examples
#'
plot3DTreeCells <- function(X,  peaks=NULL, cluid, celltype, mst=NULL, opcacity=0.7,
                            pointsize=2, labelsize=20, lineW=5, colorname="Spectral", titles="ISOMAP"){

  cluidA <- c(celltype, rep("cluster center", length(unique(cluid))))
  if(is.null(peaks)){
    Y <- getCenters(X, cluid) #   cluidA
  }else{
    Y <- X[peaks,]
    rownames(Y) <- paste("Y_", 1:nrow(Y), sep = "")
  }
  if(is.null(mst)){
    mst <- getCenMST(Y)
  }
  rX <- rbind(X, Y)

  if(dim(X)[2]==2){
    da = data.frame(x = rX[,1], y = rX[,2], z = rep(0, nrow(rX)), clu=cluidA)
  }else if(dim(X)[2]>=3){
    da = data.frame(x = rX[,1], y = rX[,2], z = rX[,3], clu=cluidA)
  }

  coid = which(rownames(brewer.pal.info)==colorname)
  colors = colorRampPalette(brewer.pal(brewer.pal.info[coid,1], colorname))(length(unique(cluidA)))

  p <- plot_ly() %>%
    add_trace(data = da, x= ~x, y= ~y, z= ~z,
              color = ~as.factor(cluidA),
              type="scatter3d", colors = colors,
              marker = list(opacity=opcacity, size=2),showlegend = TRUE)

  p <- plotly_build(p)
  markersize <- pointsize
  markerlegendsize <- labelsize
  for(i in seq(1, length(sort(unique(cluidA) )))) {
    length.group <- nrow(da[which(cluidA == sort(unique(cluidA))[i]), ])
    p$x$data[[i]]$marker$size <- c(rep(markersize,length.group), rep(c(-markersize+2*markerlegendsize), length.group))
  }

  SY <- Y[1:2,]
  edMat = get.edgelist(mst)
  ecs = colorRampPalette(brewer.pal(brewer.pal.info[coid,1], colorname))(nrow(edMat))
  for (i in 1:nrow(edMat)) {
    if(substr(edMat[i,1],1,1)=="Y"){
      st = as.integer(strsplit(edMat[i,1],'_')[[1]][2])
      ed = as.integer(strsplit(edMat[i,2],'_')[[1]][2])
    }else{
      st = as.integer(edMat[i,1])
      ed = as.integer(edMat[i,2])
    }

    SY = Y[1:2,]
    SY[1,] = Y[st,]
    SY[2,] = Y[ed,]
    if(ncol(SY)==2){
      SY = cbind(SY, rep(0, nrow(SY)))
    }
    p <- p %>%
      add_trace(x =SY[,1], y= SY[,2], z=SY[,3], mode = "lines",
                line = list(color = ecs[i], width=lineW),
                name = paste("e",st, ed, sep = "_"), inherit=FALSE )
  }
  p %>%
    layout(title = paste( titles, sep=""), legend =list(font = list(size=labelsize)))

}


#' plotCellTree: Visaulize the trajectory in 2D plot
#'
#' @param X
#'     matrix, the 3d ISOMAP
#' @param Y
#'     Default is the cluster centers as the vertexes.
#' @param cluid
#'     vector, the cluster id
#' @param mst
#'     igraph, the trajectory
#' @param root
#'     Default 1
#' @param colorby
#'     Default is "cluster"
#' @param anno
#'     other annotations
#' @param colorname
#'     colorname: rownames(brewer.pal.info): "BrBG"     "PiYG"     "PRGn"     "PuOr"     "RdBu"     "RdGy"     "RdYlBu"   "RdYlGn"   "Spectral" "Accent"   "Dark2"
#'     "Paired"   "Pastel1"  "Pastel2"  "Set1"     "Set2"     "Set3"     "Blues"    "BuGn"     "BuPu"     "GnBu"     "Greens"
#'     "Greys"    "Oranges"  "OrRd"     "PuBu"     "PuBuGn"   "PuRd"     "Purples"  "RdPu"     "Reds"     "YlGn"     "YlGnBu"
#'    "YlOrBr"   "YlOrRd"
#' @param height
#'     Default 0.2
#' @param width
#'     0.2
#' @param pointSize
#'     Default 0.2
#' @param edgeW
#'     0.5
#' @param textSize
#'     Default 0.8
#' @param legendTextS
#'     20
#' @param legendColS
#'     10
#' @param hv
#'     0.5
#' @param vv
#'     -1
#' @return
#'
#' @export
#'
#' @examples
#'
plotCellTree <- function(X, Y, cluid, mst=NULL, root=1, colorby="cluster", anno=NULL, colorName=rownames(brewer.pal.info)[22],
                         height=0.2,width=0.2, pointSize=0.2, edgeW=0.5, textSize=8, legendTextS=20, legendColS=10,
                         hv=0.5,vv=-1){
  #: display.brewer.all()
  if(is.null(mst)){
    mst <- getCenMST(Y)
  }

  leaves= which(degree(mst, v = V(mst), mode = "all")==1, useNames = T)

  layD = getLayout(mst, root=root)
  minx = 0.9*(max(layD[,1]) - min(layD[,1]))/length(leaves)
  edMat = get.edgelist(mst)
  vMat = matrix(0, nrow = nrow(edMat), ncol = 8)
  #1.get the middle points
  for(i in 1:nrow(edMat)){

    if(distances(mst, V(mst)[root], edMat[i,1]) <
       distances(mst, V(mst)[root], edMat[i,2]) ){
      O = edMat[i,1]
      B = edMat[i,2]
    }else{
      O = edMat[i,2]
      B = edMat[i,1]
    }

    vO = c(layD[B,1], layD[O,2])

    co = layD[O,1:2]
    cb = layD[B,1:2]
    vMat[i,1] = O
    vMat[i,2] = B
    vMat[i,3:4] = as.numeric(layD[O,1:2])
    vMat[i,5:6] = as.numeric(layD[B,1:2])
    vMat[i,7:8] = as.numeric(vO)
  }
  #2.get point coordination
  newD = X[,1:2]
  cc = unique(cluid)
  for(i in 1:length(cc)){
    cid = which(cluid == i)
    v = as.numeric(layD[V(mst)[[i]]$name,1:2])
    vd =  rep(v, length(cid))
    vd = t(matrix(vd, nrow = 2, ncol = length(cid)))
    newD[cid,] =vd
  }
  #plot(newD)

  #3. create edge_df
  edge_df = matrix("", nrow(vMat)*2, 6)
  for(i in 1:nrow(vMat)){
    j = (i-1)*2 + 1
    edge_df[j,1] = vMat[i,1]
    edge_df[j,2] = paste("v",vMat[i,1],sep="")
    edge_df[j,3:4] = vMat[i, 3:4]
    edge_df[j,5:6] = vMat[i, 7:8]

    j = j + 1
    edge_df[j,1] = paste("v",vMat[i,1],sep="")
    edge_df[j,2] = vMat[i,2]
    edge_df[j,3:4] = vMat[i, 7:8]
    edge_df[j,5:6] = vMat[i, 5:6]
  }
  colnames(edge_df) = c("target", "sample_name", "source_prin_graph_dim_1",
                        "source_prin_graph_dim_2", "target_prin_graph_dim_1",
                        "target_prin_graph_dim_2")

  edge_df = as.data.frame(edge_df, stringsAsFactors=FALSE)
  for(i in 3:6){
    edge_df[,i] = as.numeric(edge_df[,i])
  }

  #4.branch point
  #branchNode= names(which(degree(mst, v = V(mst), mode = "all")>=3, useNames = T))
  #branch_point_df = subset(edge_df, target %in% branchNode)[, c("sample_name", "source_prin_graph_dim_1",
  #                                                              "source_prin_graph_dim_2")]
  #branch_point_df$branch_point_idx = 1

  #5.cluster name
  cluNode = cbind(row.names(layD), layD[,1:2])
  colnames(cluNode) = c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")
  clu_df = as.data.frame(cluNode, stringsAsFactors=FALSE)

  #6. annotation
  lib_info_with_pseudo = anno
  rownames(lib_info_with_pseudo) = rownames(X)

  #7.data
  data_df <- data.frame(newD)
  row.names(data_df) <-   rownames(X)
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name",
                   by.y = "row.names")


  #8.draw
  g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1",
                                   y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
                                   yend = "target_prin_graph_dim_2"), size = edgeW,
                        linetype = "solid", na.rm = TRUE, data = edge_df)

  if(  class(data_df[, colorby])=="numeric"){
    g <- g + geom_jitter(aes_string(color = colorby), width=width,
                         size = I(pointSize), na.rm = TRUE, height = height) +
      scale_colour_distiller(palette=colorName,guide = "colourbar",
                             aesthetics = "colour", direction = 1)


  }else{
    #g <- g + geom_jitter(aes_string(color = colorby), width=width,
    #                     size = I(pointSize), na.rm = TRUE, height = height)  +
    #  scale_colour_brewer(palette = colorName)

    cN = length(unique(anno[[colorby]]))
    #pt = brewer.pal.info[colorName,1]
    #cols = colorRampPalette(brewer.pal(pt, colorName))(cN)

    if(cN > brewer.pal.info[colorName,1]){
      cols = colorRampPalette(brewer.pal(brewer.pal.info[colorName,1],colorName))(cN)
    }else{
      cols = brewer.pal(cN, colorName)
    }

    g <- g + geom_jitter(aes_string(color = colorby), width=width,
                         size = I(pointSize), na.rm = TRUE, height = height) +
      scale_color_manual(values = cols)
  }

  #g

  g <- g + geom_text(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2",
                                label = "sample_name"), hjust = hv, vjust = vv,
                     size = textSize, color = "black", na.rm = TRUE,
                     data = clu_df)

  g <- g + theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + xlab("") + ylab("") +
    theme(legend.text=element_text(size=legendTextS), legend.title=element_text(size=legendTextS))+
    theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(line = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank())
  if(class(data_df[, colorby])=="numeric"){

  }else{
    g = g +  guides(colour = guide_legend(override.aes = list(size=legendColS)))
  }
  g
}



#' getLayout: get layout for the trajectory in 2D
#'
#' @param mst
#'     igraph, the trajectory
#' @param root
#'     Default 1
#' @return
#'     igraph.layout
#' @export
#'
#' @examples
#'
getLayout <- function(mst, root=1){
  #plot(mst, layout=layout_as_tree(mst, root=root))
  ####1. convert undirect to direct based root
  dmst = as.directed(mst, mode = "mutual" )

  ##delete edges direct to root
  nb = neighbors(dmst, v=root, mode =  "all")
  len = length(nb)
  nn = rep("",len)
  for(i in 1:len){
    nn[i] = nb[[i]]$name
  }
  nn = unique(nn)
  dnb = bfs(dmst, root=root, neimode="out")
  od = as.vector(dnb$order)
  for(i in 2:length(od)){
    ch = od[i]
    for(j in 1:(i-1)){
      fa = od[j]
      if( V(dmst)[fa] %in% neighbors(dmst, ch)){
        dmst = delete_edges(dmst, edges= paste(V(dmst)[ch]$name, "|", V(dmst)[fa]$name, sep="") )
      }
    }
  }
  #plot(dmst)
  #plot(dmst, layout=layout_as_tree(dmst, root=root), vertex.size=1,vertex.label.cex = 1, edge.arrow.width = 0.5, edge.arrow.size=0.5,)
  ########2.get layout

  igraph.layout <- create_layout(dmst, layout = "dendrogram")
  rownames(igraph.layout) =igraph.layout$name
  igraph.layout
  # plot(dmst, layout= as.matrix(igraph.layout[,1:2]) )
}

#g <- make_ring(10) + make_full_graph(5)
#coords <- layout_(g, as_star())
#plot(g, layout = coords)
