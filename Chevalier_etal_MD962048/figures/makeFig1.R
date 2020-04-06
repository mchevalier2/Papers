## Figure 1: The CREST method
##
## Loading necessary data
load(url('https://github.com/mchevalier2/Papers/raw/master/Chevalier_Limpopo/figures/Fig1.RData'))
##
## Parameterizing
NREP=500
BIN_WIDTH=2
OUTPUT_FOLDER=getwd()
##
## Fornmtting data
spclim=c(as.data.frame(clim)[sp1==1,1],as.data.frame(clim)[sp2==1,1],as.data.frame(clim)[sp3==1,1],as.data.frame(clim)[sp4==1,1])
spname=c(rep('sp1', sum(sp1)),rep('sp2', sum(sp2)),rep('sp3', sum(sp3)),rep('sp4', sum(sp4)))
spclim_unique=as.data.frame(clim)[(sp1+sp2+sp3+sp4)>=1,1]
##
xx=seq(-8, 10, length.out=500)
ccs=calib_clim_space(clim_space=values(raster::raster(clim)), bin_width=BIN_WIDTH)
##
weights=sqrt(c(sum(sp1),sum(sp2),sum(sp3),sum(sp4)))
pdfsp=as.matrix(do.call(cbind, tapply(spclim, spname, function(x) return(pdf_sp(x, bin_width=BIN_WIDTH, ccs=ccs, xx=xx, shape='normal')))))
pdfpol=((pdfsp) %*% weights)/sum(weights)
##
## Fitting randomized pdfs
pdfrand=array(0, dim=c(NREP, length(xx)))
for(i in 1:NREP){
  d=as.data.frame(rclim[i])
  sprclim=c(d[sp1==1,1],d[sp2==1,1],d[sp3==1,1],d[sp4==1,1])
  rpdfsp=as.matrix(do.call(cbind, tapply(sprclim, spname, function(x) return(pdf_sp(x, bin_width=BIN_WIDTH, ccs=ccs, xx=xx, shape='normal')))))
  pdfrand[i,]=(rpdfsp %*% weights)/sum(weights)
}
##
## Calculating indices
output=matrix(0,ncol=5,nrow=NREP)
colnames(output)=c("optimum", "tolerance", "skewness", "value", "pvalue")
for(i in 1:NREP){
  output[i,1:3]=calculate_indices(pdfrand[i,], xx)
}
test=calculate_indices(pdfpol, xx)
##
fhat1=ks::kde(output[,1:3], gridsize=c(100,100,100))
##
h=hist(sortedclim,plot=FALSE, breaks=seq(-6, 9,1))
h2=hist(spclim_unique,plot=FALSE, breaks=h$breaks)
##
##
pdf(paste0(OUTPUT_FOLDER, '/Chevalier_etal_MD962048_Fig1.pdf", width=7.54, height=7.54/8*2.46, useDingbats=FALSE)  ;  {
  par(mar=c(0,0,0,0))
  par(oma=c(0,0,0,0))
  par(ps=7*1.5)
  layout(matrix(c(12,12,13,13,14,14,15,15,
                  1,1,2,3,6,7,8,8,
                  1,1,4,5,6,6,8,8,
                  9,9,9,9,10,10,11,11), byrow=TRUE, ncol=8),
      height=c(0.08, 0.5,0.5,0.15),
      width=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5))
  ##
  ##-- PANE A
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(9.5,35.5), ylim=c(-35.5,-9.5), xaxs="i",yaxs="i")
  image(clim,add=TRUE, col=BlaYel, zlim=c(-6,9))
  rect(10,-10, 35,-35, lwd=0.8)
  ##
  XY=coordinates(clim)
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(9,35.5), ylim=c(-35.5,-9), xaxs="i",yaxs="i")
  image(clim,add=TRUE, col=makeTransparent(BlaYel,alpha=0.5), zlim=c(-6,9))
  for(w in which(sp1==1)) rect(XY[w,1]-0.125, XY[w,2]-0.125, XY[w,1]+0.125, XY[w,2]+0.125, col=SPCOL[1], border=NA)
  rect(10,-10, 35,-35, lwd=0.8)
  rect(22.5,-34.5,34.5,-31.5, col=makeTransparent('white', alpha=0.8), lwd=0.4)
  text(28.5,-33, 'Species 1', cex=6/7, col=SPCOL[1])
  ##
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(9.5,36), ylim=c(-35.5,-9), xaxs="i",yaxs="i")
  image(clim,add=TRUE, col=makeTransparent(BlaYel,alpha=0.5), zlim=c(-6,9))
  for(w in which(sp2==1)) rect(XY[w,1]-0.125, XY[w,2]-0.125, XY[w,1]+0.125, XY[w,2]+0.125, col=SPCOL[2], border=NA)
  rect(10,-10, 35,-35, lwd=0.8)
  rect(22.5,-34.5,34.5,-31.5, col=makeTransparent('white', alpha=0.8), lwd=0.4)
  text(28.5,-33, 'Species 2', cex=6/7, col=SPCOL[2])
  ##
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(9,35.5), ylim=c(-36,-9.5), xaxs="i",yaxs="i")
  image(clim,add=TRUE, col=makeTransparent(BlaYel,alpha=0.5), zlim=c(-6,9))
  for(w in which(sp3==1)) rect(XY[w,1]-0.125, XY[w,2]-0.125, XY[w,1]+0.125, XY[w,2]+0.125, col=SPCOL[3], border=NA)
  rect(10,-10, 35,-35, lwd=0.8)
  rect(22.5,-34.5,34.5,-31.5, col=makeTransparent('white', alpha=0.8), lwd=0.4)
  text(28.5,-33, 'Species 3', cex=6/7, col=SPCOL[3])
  ##
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(9.5,36), ylim=c(-36,-9.5), xaxs="i",yaxs="i")
  image(clim,add=TRUE, col=makeTransparent(BlaYel,alpha=0.5), zlim=c(-6,9))
  for(w in which(sp4==1)) rect(XY[w,1]-0.125, XY[w,2]-0.125, XY[w,1]+0.125, XY[w,2]+0.125, col=SPCOL[4], border=NA)
  rect(10,-10, 35,-35, lwd=0.8)
  rect(22.5,-34.5,34.5,-31.5, col=makeTransparent('white', alpha=0.8), lwd=0.4)
  text(28.5,-33, 'Species 4', cex=6/7, col=SPCOL[4])
  ##
  XXrange=range(xx)  ;  XXrange[1]=XXrange[1]-0.15*diff(XXrange)  ;  XXrange[2]=XXrange[2]+0.02*diff(XXrange)
  YYrange=c(0,max(pdfsp)*1.02)  ;  YYrange[1]=YYrange[1]-0.02*diff(YYrange)  ;  YYrange[2]=YYrange[2]+0.02*diff(YYrange)
  plot(0,0, type='n', axes=FALSE,frame=FALSE, xlim=XXrange, ylim=YYrange, xaxs="i",yaxs="i")
  for(sp in c(4,2,3,1)){
    polygon(xx,pdfsp[,sp], col=makeTransparent(SPCOL[sp],alpha=1), border=NA)
  }
  for(sp in c(4,2,3,1)){
    polygon(xx,pdfsp[,sp], col=NA, border='black', lwd=0.8)
  }
  polygon(xx,pdfpol,col=makeTransparent('black', alpha=0.8), lwd=1, border='white')
  segments(xx[1],0,rev(xx)[1],0, lwd=1)
  polygon(c(-8.5, -8.8, -9.3, -8.5, -7.7, -8.2, -8.5),c(0.03, .39, 0.39, 0.47,0.39, 0.39,0.03), col='black')
  text(-9.8, 0.23, 'Density of Probability', adj=c(0.5, 0.5), srt=90)
  ##
  XXrange=range(h$breaks)  ;  XXrange[1]=XXrange[1]-0.1*diff(XXrange)  ;  XXrange[2]=XXrange[2]+0.02*diff(XXrange)
  YYrange=c(0,max(h$counts))  ;  YYrange[1]=YYrange[1]-0.15*diff(YYrange)  ;  YYrange[2]=YYrange[2]+0.02*diff(YYrange)
  plot(0,0, type='n', axes=FALSE,frame=FALSE, xlim=XXrange, ylim=YYrange)
  for(i in 1:length(h$density)) rect(h$breaks[i],0,h$breaks[i+1], h$counts[i], lwd=0.8 )
  polygon(rep(h2$breaks,each=2),c(0,rep(h2$counts, each=2),0), col='black' )
  polygon(c(-5,6,6,8,6,6,-5), c(-80,-55,-20,-80,-140,-105,-80), col='black', border=NA)
  text(1.5, -210, "Environmental variable", cex=6/7, adj=c(0.5,0.5))
  polygon(c(-6, -6.3, -6.7, -6, -5.3, -5.7, -6), c(0.03, .39, 0.39, 0.47,0.39, 0.39,0.03)*3500, col='black', border=NA)
  text(-7.5, 800, 'Abundance', cex=6/7, adj=c(0.5, 0.5), srt=90)
  ##
  XXrange=range(xx[which(pdfpol>0.0001)])  ;  XXrange[1]=XXrange[1]-0.15*diff(XXrange)  ;  XXrange[2]=XXrange[2]+0.02*diff(XXrange)
  YYrange=c(0,max(pdfpol)*1.05)  ;  YYrange[1]=YYrange[1]-0.02*diff(YYrange)  ;  YYrange[2]=YYrange[2]+0.02*diff(YYrange)
  plot(0,0, type='n', axes=FALSE,frame=FALSE, xlim=XXrange, ylim=YYrange, xaxs="i",yaxs="i")
  polygon(xx,pdfpol,col=makeTransparent('black', alpha=0.9), lwd=1)
  shape::Arrows(test[1]-10*test[3], max(pdfpol)*1.11, test[1], max(pdfpol)*1.11, code=3, arr.type='triangle', arr.adj=1, arr.length=0.1)
  segments(test[1],0,test[1], max(pdfpol), lwd=1, col='white')
  segments(test[1]-3*test[3],0,test[1]-3*test[3], pdfpol[which(xx>test[1]-3*test[3])[1]-1], lwd=1, col='white')
  shape::Arrows(xx[rev(which(pdfpol>=max(pdfpol)/4))[1]], max(pdfpol)/4, xx[which(pdfpol>=max(pdfpol)/4)[1]], max(pdfpol)/4, code=3, arr.type='triangle', arr.adj=1, arr.length=0.1, col="white")
  text(test[1]-5*test[3], max(pdfpol)*1.16, 'Skewness', cex=1, adj=c(0.5,0))
  text(test[1]+2*test[3], max(pdfpol)*2/3, 'Optimum', cex=1, adj=c(0.5,0), srt=90, col="white")
  text(test[1]-5*test[3], pdfpol[which(xx>test[1]-3*test[3])[1]-1]/2, 'Mean', cex=1, adj=c(0.5,0.9), srt=90, col="white")
  rect(xx[rev(which(pdfpol>=max(pdfpol)/4))[1]], max(pdfpol)*0.22, xx[which(pdfpol>=max(pdfpol)/4)[1]], max(pdfpol)*0.16, col=makeTransparent('black', alpha=0.8))
  text(mean(c(xx[rev(which(pdfpol>=max(pdfpol)/4))[1]], xx[which(pdfpol>=max(pdfpol)/4)[1]])), max(pdfpol)*0.22, 'Uncertainty', cex=1, adj=c(0.5,1), col="white")
  polygon(c(-8.2, -8.4, -8.7, -8.2, -7.7, -8, -8.2),c(0.03, .39, 0.39, 0.47,0.39, 0.39,0.03)*0.572, col='black')
  text(-9, 0.13, 'Density of Probability', adj=c(0.5, 0.5), srt=90)
  ##
  par(mar=c(0,0,0,0))
  plot(0,0, type='n', axes=FALSE,frame=FALSE, xlim=c(-0.5,1.5), ylim=c(0.5,1.4), xaxs="i",yaxs="i")
  for(i in 1:20){
    rect((i-1)/20, 0.7, i/20,1, col=BlaYel[i], border=NA)
  }
  rect(0,0.7,1,1, col=NA, lwd=0.3)
  text(0.5, 1.1, 'Environmental gradient', adj=c(0.5, 0), cex=1)
  text(0, 1.1, '-', adj=c(0.5, 0), cex=8/7, font=2)
  text(1, 1.1, '+', adj=c(0.5, 0), cex=8/7, font=2)
  ##
  plot(0,0, type='n', axes=FALSE,frame=FALSE, xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  polygon(c(0.15, 0.80,0.80,0.95,0.80,0.80,0.15), c(0.8, 0.87,0.98,0.8,0.62, 0.73, 0.8), col='black')
  text(0.45, 0.5, "Environmental variable", cex=1, adj=c(0.5,0.5))
  plot(0,0, type='n', axes=FALSE,frame=FALSE, xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  polygon(c(0.15, 0.80,0.80,0.95,0.80,0.80,0.15), c(0.8, 0.87,0.98,0.8,0.62, 0.73, 0.8), col='black')
  text(0.45, 0.5, "Environmental variable", cex=1, adj=c(0.5,0.5))
  ##
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  text(0.5, 0.45, 'A. Environmental gradient', cex=1, font=2, adj=c(0.5, 0.5))
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  text(0.5, 0.45, 'B. Species distributions', cex=1, font=2, adj=c(0.5, 0.5))
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  text(0.5, 0.45, "C. Species' and pollen taxon's pdfs", cex=1, font=2, adj=c(0.5, 0.5))
  plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  text(0.5, 0.45, "D. Posterior distribution of climate", cex=1, font=2, adj=c(0.5, 0.5))
}  ;  dev.off()
##
##
##
##
##
