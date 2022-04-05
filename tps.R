#' tps
#'
#'
#' shade.scale,
#' reference grids, optional grid plotting,
#' different shading methods (contours, arrows)
#' convex hulls...
#' groups of landmarks
#' ellipses, individual landmark variation
#' smooth (smooth colors under grid)
#'


tps <- function(target.lm,
                reference.lm,
                A=NULL,
                shade.scale=FALSE,
                scale.lm = 1, # what to scale the whole landmark configuration by
                mag = 1,
                add = FALSE,
                at = c(0,0), # where to add the TPS grids, if anywhere
                snap = TRUE,
                #snap.buffer = c(0,0,0,0) , # how many rows or columns to add onto the snap as a buffer... left, right, down, up (xmin, xmax, ymin, ymax)
                shade = FALSE,
                shade.trans = 0.7,
                palette = NULL,
                n.grid.col = 24,
                n.grid.row = NULL,
                grid.aes = list(col = "black", lwd = 1, trans = 1),
                links = NULL,
                link.aes = list(col = "black", lwd = 3, trans = 1),
                plot.lm = TRUE,
                lm.aes = list(col = "black", cex = 1, trans = 1),
                plot.ref.links = FALSE,
                ref.link.aes = list(col = "gray", lwd = 3, trans = 0.85),
                plot.ref.lm = FALSE,
                ref.lm.aes = list(col = "gray", cex = 1, trans = 0.5),
                plot.arrows = FALSE,
                arrow.aes = list(col = "black", lwd = 1, trans = 1, code = 2, length = 0.1),
             #   smooth = FALSE,
                buffer.x = 0.3,
                buffer.y = 0.3){

require(pracma) ; require(colorBlindness)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


# from 'Morphometrics with R' by Claude
# computes deformation of grid
  tps2d<-function(M, matr, matt){p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
  P<-matrix(NA, p, p)
  for (i in 1:p)
     {for (j in 1:p){
         r2<-sum((matr[i,]-matr[j,])^2)
         P[i,j]<- r2*log(r2)}}
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
  m2<-rbind(matt, matrix(0, 3, 2))
  coefx<-solve(L)%*%m2[,1]
  coefy<-solve(L)%*%m2[,2]
  fx<-function(matr, M, coef)
     {Xn<-numeric(q)
      for (i in 1:q)
           {Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
           Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
      Xn}
  matg<-matrix(NA, q, 2)
  matg[,1]<-fx(matr, M, coefx)
  matg[,2]<-fx(matr, M, coefy)
  matg}

  # my own internal function to make an array of the cells
  # cells -> array ==== celery  :-)
celery <- function(X){
  for(i in 1:n.row){
  rows <- array(dim=c(4,2,n.col))
  for(j in 1:n.col){
  rows[1,,j] <- X[1+(j-1)+((i-1)*n),]
  rows[2,,j] <- X[1+(j-1)+n+((i-1)*n),]
  rows[3,,j] <- X[1+j+n+((i-1)*n),]
  rows[4,,j] <- X[1+j+((i-1)*n),]
  }
if(i==1){cell.array<-rows}else{
cell.array <- bindArr(cell.array,rows,along=3)}
  }
  return(cell.array)
}

# adds transparency to a color code
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}






  reference.lm <- reference.lm*scale.lm

    # rotating the target landmark configuration onto the reference
  target.lm <- rotonto(x = reference.lm, y = target.lm, scale = TRUE)$yrot

  # incorporating the magnitude of the shape deformation
  # by modifying the 'target.lm' object
  if(mag != 1){
    target.lm <- reference.lm+((target.lm-reference.lm)*mag)
  }


  #
  if(!add){at <- c(0,0)}
  n <- n.grid.col
  s <- 1/1#zoom
  matr <- translate.lm(reference.lm, at)
  matt <- translate.lm(target.lm, at)
  M2 <- matt
  xm<-min(matt[,1])*s # x min
  ym<-min(matt[,2])*s # x max
  xM<-max(matt[,1])*s # y min
  yM<-max(matt[,2])*s # y max
  rX<-xM-xm; rY<-yM-ym # ranges


  pad.x <- buffer.x*rX
  pad.y <- buffer.y*rY

  a<-seq(xm-pad.x, xM+pad.x, length=n) # x values of grid points

  # if the number of rows isn't defined, this calculates it
  # if it is, this assigns n.grid.row to object 'm'
  if(is.null(n.grid.row)){m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))}else{(m <- n.grid.row)}


  b<-seq(ym-pad.y, yM+pad.y,length=m) # y values of grid points
  M<-as.matrix(expand.grid(a,b))

  ngrid<-tps2d(M,matr,matt)


  n.col <- (n-1)
  n.row <- length(b)-1
  n.cell <- n.col*n.row


# warped array of cells
  # storing the grid-points that define each cell in an array
cell.array <- celery(X=ngrid)

# reference cells
cell.array.ref <- celery(X=M)

# calculating the individual area of a cell in the reference (unwarped) configuration
# for use later, but it feels right to compute it here
area.reference <- abs(polyarea(cell.array.ref[,1,1],cell.array.ref[,2,1]))




# if 'snap' is true, this section trims the grid to only fit cells with landmarks inside
  # probably a much faster way of doing this
if(snap){
xmr <- min(matr[,1])
ymr <- min(matr[,2])
xMr <- max(matr[,1])
yMr <- max(matr[,2])

grid.id <- expand.grid(1:n.col,1:n.row) # making row and columns id's of each grid cell

xm.idr <- which(matr[,1] == xmr)
ym.idr <- which(matr[,2] == ymr)
xM.idr <- which(matr[,1] == xMr)
yM.idr <- which(matr[,2] == yMr)

dists.xmr <- c()
dists.ymr <- c()
dists.xMr<- c()
dists.yMr <- c()
for(i in 1:n.cell){
  dists.xmr[[i]] <- euc.dist(matr[xm.idr,],centroid1(cell.array.ref[3:4,,i]))
  dists.ymr[[i]] <- euc.dist(matr[ym.idr,],centroid1(cell.array.ref[2:3,,i]))
  dists.xMr[[i]] <- euc.dist(matr[xM.idr,],centroid1(cell.array.ref[1:2,,i]))
  dists.yMr[[i]] <- euc.dist(matr[yM.idr,],centroid1(cell.array.ref[c(1,4),,i]))
  }
which.xmr <- which.min(unlist(dists.xmr))
which.ymr <- which.min(unlist(dists.ymr))
which.xMr <- which.min(unlist(dists.xMr))
which.yMr <- which.min(unlist(dists.yMr))

col.min.idr <- grid.id[which.xmr,1]#-snap.buffer[1]
col.max.idr <- grid.id[which.xMr,1]#+snap.buffer[2]
row.min.idr <- grid.id[which.ymr,2]#-snap.buffer[3]
row.max.idr <- grid.id[which.yMr,2]#+snap.buffer[4]

min.cell <- cell.array.ref[1,,which(grid.id[,1]==col.min.idr & grid.id[,2]==row.min.idr)]
max.cell <- cell.array.ref[3,,which(grid.id[,1]==col.max.idr & grid.id[,2]==row.max.idr)]

grid.to.drop <- which(M[,1] < min.cell[1] | M[,1] > max.cell[1] | M[,2] < min.cell[2] | M[,2] > max.cell[2])
cells.to.drop <- which(grid.id[,1] < col.min.idr | grid.id[,1] > col.max.idr | grid.id[,2] < row.min.idr | grid.id[,2] > row.max.idr)

cell.array <- cell.array[,,-cells.to.drop]
n.cell <- n.cell-length(cells.to.drop)

ngrid <- ngrid[-grid.to.drop,]
n <- length(col.min.idr:col.max.idr)+1
m <- length(row.min.idr:row.max.idr)+1
}







# if there's an array:
  # need to consolidate this into a couple of clean functions.
if(!is.null(A)){
  if(is.matrix(A)){arrayspecs(A=A,p=dim(target.lm)[1],k=dim(target.lm)[2])} # if a matrix is provided, this switches it to an array

# landmark variation
target.array <- A
all.diff <- c()

for(i in 1:dim(target.array)[3])
  {

  target.lm1 <- rotonto(x = reference.lm, y = target.array[,,i], scale = TRUE)$yrot

  # incorporating the magnitude of the shape deformation
  # by modifying the 'target.lm' object
  if(mag != 1){
    target.lm1 <- reference.lm+((target.lm1-reference.lm)*mag)
  }

  #
  if(!add){at <- c(0,0)}
  matt1 <- target.lm1*scale.lm + at

  ngrid1<-tps2d(M,matr,matt1)


  # warped array of cells
  # storing the gridpoints that define each cell in an array
cell.array1 <- celery(X=ngrid1)

cell.array1 <- cell.array1[,,-cells.to.drop]

diff1 <- c()
for(i in 1:n.cell){
  area.target1 <- abs(polyarea(cell.array1[,1,i],cell.array1[,2,i]))
  diff1[[i]] <- area.target1-area.reference
  }
diff1 <- as.vector(unlist(diff1))


if(i == 1){all.diff <- diff1}
if(i > 1){all.diff <- c(all.diff, diff1)}
}

cs <- range(unlist(all.diff)) # the color scale, taking into account all the variation in the landmark data provided

}
if(is.null(A)){}




# setting the color palette according to area of the cells
if(!is.null(palette) | shade){

  # setting the pallette if none is given:
  if(is.null(palette)){shade.cols <- colorRampPalette(Blue2DarkRed18Steps)(n.cell)}else{shade.cols <- palette(n.cell)}
  idcolor <- addTrans(shade.cols, 255*shade.trans)

diff <- numeric(n.cell)
for(i in 1:n.cell){
  area.target <- abs(pracma::polyarea(cell.array[,1,i],cell.array[,2,i]))
    diff[i] <- (area.target-area.reference)
  }
#if(ABS){diff <- abs(rescale.numeric(diff,to=c(-1,1)))} # if you want small cells to be the same color as larger cells


if(shade.scale){nm <- max(abs(cs))*1.01 ; brk <- seq(-nm,nm,length.out = n.cell)}else{
  nm <- max(abs(range(diff)))*1.01 ; brk <- seq(-nm,nm,length.out = n.cell)
}


dit2<-cut(diff,breaks = length(brk))
}


# plotting:

# if you aren't adding to an existing plot, this just plots a lil canvas
if(!add){
plot(NA, asp=1,axes=F,xlab="",ylab="",xlim = range(ngrid[,1]), ylim = range(ngrid[,2]))
  }

# move this to be below everything else in the plotting section !
if(!is.null(palette) | shade){

# plotting the shade
for(i in 1:n.cell){
  if(shade.scale){j <- i+1} else{j <- i}
  polygon(x = cell.array[,,i], col = idcolor[as.numeric(dit2[i])], border = NA)
}

}


# plotting the grids
for (i in 1:m){lines(ngrid[(1:n)+(i-1)*n,],lwd=grid.aes$lwd, col=addTrans(grid.aes$col, grid.aes$trans*255))}
for (i in 1:n){lines(ngrid[(1:m)*n-i+1,],lwd=grid.aes$lwd, col=addTrans(grid.aes$col, grid.aes$trans*255))}


# plotting the landmarks. reference landmarks under the target landmarks
if(plot.ref.lm){points(matr,pch = 19,cex=ref.lm.aes$cex, col = addTrans(ref.lm.aes$col, ref.lm.aes$trans*255))}
if(plot.lm){points(M2,pch = 19,cex=lm.aes$cex, col = addTrans(lm.aes$col, lm.aes$trans*255))}


if(plot.ref.links){
for (i in 1:nrow(links)) {
segments(matr[links[i, 1], 1], matr[links[i, 1], 2], matr[links[i, 2], 1], matr[links[i, 2], 2],
         lwd = ref.link.aes$lwd, col = addTrans(ref.link.aes$col, ref.link.aes$trans*255))}
  }
if(!is.null(links)){
for (i in 1:nrow(links)) {
segments(M2[links[i, 1], 1], M2[links[i, 1], 2], M2[links[i, 2], 1], M2[links[i, 2], 2],
         lwd = link.aes$lwd, col = addTrans(link.aes$col, link.aes$trans*255))}
  }

if(plot.arrows){
  for(i in 1:nrow(M2)){
    arrows(x0 = matr[i,1], y0 = matr[i,2], x1 = M2[i,1], y1 = M2[i,2], code = arrow.aes$code,
           col = addTrans(arrow.aes$col, arrow.aes$trans*255), lwd = arrow.aes$lwd, length = arrow.aes$length)}
  }

}



# examples
library(geomorph)
data("plethodon")
LM <- plethodon$land
proc <- gpagen(LM)
A <- proc$coords

# basic
tps(A[,,1], mshape(A), links = plethodon$links)

# with shade
tps(A[,,1], mshape(A), links = plethodon$links, shade = TRUE)

# prettier
tps(A[,,1], mshape(A), links = plethodon$links, plot.ref.links = TRUE, plot.ref.lm = TRUE,
    shade = T, shade.trans = 0.33, grid.aes = list(col = 'darkgray', trans = 0.8),
    ref.lm.aes = list(col = "darkgray", trans = 1), ref.link.aes = list(col = "darkgray", trans = 1, lwd = 3))

# with arrows
tps(A[,,1], mshape(A), mag = 1.5, links = plethodon$links, plot.ref.links = TRUE, plot.ref.lm = TRUE,
    ref.lm.aes = list(col = "darkgray", trans = 1), ref.link.aes = list(col = "darkgray", trans = 1, lwd = 3),
    plot.arrows = TRUE, arrow.aes = list(code = 2, col = "red", length = 0.1, lwd = 2))


# on a PCA plot of the salamanders
pca <- gm.prcomp(A = A)

g1 <- which(pca$x[,1] > 0 & pca$x[,2] < 0.02)
g2 <- which(pca$x[,1] < 0 & pca$x[,2] < 0)
g3 <- which(pca$x[,2] > 0.02)
gs <- c(g1,g2,g3)
col.gp <- numeric(dim(A)[3])
col.gp[g1] <- "dodgerblue" ; col.gp[g2] <- "sienna2" ; col.gp[g3] <- "gray"
G.x <- rbind(centroid1(pca$x[g1,1:2]),centroid1(pca$x[g2,1:2]),centroid1(pca$x[g3,1:2]))
ats <- rbind(c(0,-0.03), c(-0.07, -0.01), c(0.03,0.06))
ramps <- c(colorRampPalette(c("dodgerblue",'white',"dodgerblue")),
           colorRampPalette(c("sienna2",'white',"sienna2")),
           colorRampPalette(c("black",'white',"black")))

plot(pca$x[,1:2], pch = 19, col = col.gp, cex = 1.33, lwd = 2, asp = T, xlab = "PC1", ylab = "PC2") ; grid()
points(G.x, pch = 22, bg = c('dodgerblue','sienna2','gray'), cex = 3, lwd = 2)
#tps(mshape(A),mshape(A),at =c(0,0), add=T, scale.lm = 0.05, links=plethodon$links)
for(i in 1:NROW(G.x)){
  M <- mshape(A[,,gs[[i]]])
  tps(target.lm = M,
      reference.lm = mshape(A),
      A = A, shade.scale = TRUE,
      palette = ramps[[i]],
      scale.lm = 0.067, mag = 1.2,
      shade = T, shade.trans = 0.25,
      add = T, at = ats[i,], plot.lm = FALSE,
      links = plethodon$links,
      plot.ref.links = T,
      grid.aes = list(col = 'gray', trans = 0.5))
}




# grid of theoretical shapes in a moprhospace
nx <- 5 ; ny <- 5
xs <- seq(range(pca$x[,1])[1], range(pca$x[,1])[2], length.out = nx)
ys <- seq(range(pca$x[,2])[1], range(pca$x[,2])[2], length.out = ny)
scores <- expand.grid(xs,ys) ; #plot(scores)
shapes <- matrix(NA, ncol = dim(A)[1]*dim(A)[2], nrow = nrow(scores))
PC <- pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(proc$consensus))+pca$rotation[,1]*scores[i,1]+pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 2, p = 12)


# \/ time consuming:

# plot(NA, xlim = range(pca$x[,1])*1.2, ylim=range(pca$x[,2])*1.2, asp = T,
#      xlab = "PC1", ylab = "PC2") ; grid()
# for(i in 1:NROW(scores)){
#   tps(target.lm = shapes1[,,i],
#       reference.lm = mshape(A),
#       n.grid.col = 18,
#       A = A, shade.scale = TRUE,
#       scale.lm = 0.033, mag = 1,
#       shade = T, shade.trans = 0.25,
#       add = T, at = as.matrix(scores)[i,],
#       links = plethodon$links,
#       link.aes = list(lwd = 1, col = "black"),
#       grid.aes = list(col = 'gray', trans = 0.5),
#       plot.lm = FALSE)
# }

