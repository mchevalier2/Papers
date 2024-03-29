""## Figure 2: Study area and calibration dataset
##
## Loading necessary data
load(url('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/figures/Fig2.RData'))

OUTPUT_FOLDER=getwd()
s <- readline(prompt=paste0("Where should the figure be saved?\nDefault is current workin directory (",OUTPUT_FOLDER,"): "))
if(s != '') OUTPUT_FOLDER <- s

pkg2install=c()
if (! ("sp" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'sp')
if (! ("raster" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'raster')


makePlot <- TRUE

if (length(pkg2install) > 0){
    s=''
    while (! s %in% c('y', 'yes', 'Y', 'YES', 'n', 'no', 'N', 'NO')){
        s <- readline(prompt=paste0("The following are required: ", paste(pkg2install, collapse=', '),". Do you want to install them? [yes/no] "))
    }
    if(s %in% c('y', 'yes', 'Y', 'YES')){
        install.packages(pkg2install)
    }else{
        print("Script aborded.")
        makePlot <- FALSE
    }
}

if (makePlot) {
    library(raster)

    EXT <- c(10,42,-40,5)
    EXT_CATCH <- c(24.575, 35.5, -28, -20)

    makeTransparent <- function(..., alpha=0.5) {
        if(alpha>1) alpha=1
        if(alpha<0) alpha=0
        alpha = floor(255*alpha)
        newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
        .makeTransparent = function(col, alpha) {
          rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
        }
        newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
        return(newColor)
    }

    ## Colours
    RdYlBu=rev(colorRampPalette(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))(12))
    heatcol4=c(rgb(254,204,92,maxColorValue=255),
               rgb(253,141,60,maxColorValue=255),
               rgb(240,59,32,maxColorValue=255),
               rgb(189,0,38,maxColorValue=255))
    DIV.10col=c(rgb(255,255,179,maxColorValue=255),
                rgb(141,211,199,maxColorValue=255),
                rgb(190,186,218,maxColorValue=255),
                rgb(251,128,114,maxColorValue=255),
                rgb(128,177,211,maxColorValue=255),
                rgb(253,180,98,maxColorValue=255),
                rgb(179,222,105,maxColorValue=255),
                rgb(252,205,229,maxColorValue=255),
                rgb(188,128,189,maxColorValue=255)
              )[c(3,8,8,8,8,8,7,7,5,4,4,6,1,3,9)]

    ## Plot
    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_Fig2.pdf"), width=7.54, height=7.54/3*(EXT[4]-EXT[3])/(EXT[2]-EXT[1])+1.9)  ;  {
        layout(matrix(c(1,3,5,2,4,6), ncol=3, byrow=TRUE), width=1, height=c(7.54/3*(EXT[4]-EXT[3])/(EXT[2]-EXT[1]), 1.9))
        par(mar=c(0,0,0,0))

        ##-- PANE A
        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT[1:2]+c(-1.2,-1.2), ylim=EXT[3:4])
        for(i in seq(-35,5,5)) {
            if(i%%10 == 0) segments(EXT[1], i, EXT[2], i, lwd=0.3, col='grey50')
            text(EXT[1]-0.3, i, i, adj=c(1,0.4), cex=0.8)
        }
        for(i in seq(10,40,5)) {
            text(i, EXT[4]+0.3, i, adj=c(0,0), cex=0.8)
        }
        image(BIO1, col=RdYlBu, zlim=c(8,32), breaks=seq(8,32,2), add=TRUE, legend=FALSE, interpolate=FALSE)
        plot(LAKES, add=TRUE, col="grey50", border=NA)
        plot(M1, lwd=0.5, add=TRUE)
        plot(WATERSHED, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='black')
        plot(WATERSHED2, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='black')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=1.5)
        text(33, -27.5, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        points(34+26/60, -11-18/60,pch=23, col='deeppink4', bg='deeppink', cex=2, lwd=1.2)
        text(35.5, -11-18/60, 'Lake\nMalawi', cex=1, adj=c(0, 0.6), font=2)
        rect(EXT[1], EXT[3]-20, EXT[2], EXT[4], col=NA,border='black', lwd=0.3)
        for(i in seq(2,11,1)){  rect(EXT[1]+6+(i-2)*2, -37.5, EXT[1]+6+(i-2)*2+2, -38.5, col=RdYlBu[i], border=NA)  }
        rect(EXT[1]+6, -37.5, EXT[1]+26, -38.5, col=NA, border='black', lwd=0.3)
        for(i in seq(2,12,1)){
            segments(EXT[1]+6+(i-2)*2, -37.5, EXT[1]+6+(i-2)*2, -38.5, lwd=0.4)
            text(EXT[1]+6+(i-2)*2, -38.9, i*2+6, adj=c(0.5, 1), cex=0.8)
        }
        polygon(c(EXT[1]+4, EXT[1]+5.5, EXT[1]+5.5), c(-38, -37.5, -38.5), lwd=0.4, col=RdYlBu[1])
        polygon(c(EXT[2]-4, EXT[2]-5.5, EXT[2]-5.5), c(-38, -37.5, -38.5), lwd=0.4, col=RdYlBu[12])
        text(EXT[1]+(EXT[2]-EXT[1])/2, -36.5, 'Mean Annual Temperature (°C)', adj=c(0.5, 0.5), cex=0.9)
        rect(EXT[1]+0.8, EXT[4]-0.8, EXT[1]+3.5, EXT[4]-3.8, border="black", lwd=0.4, col="white")
        text(EXT[1]+0.8+(3.5-0.8)/2, EXT[4]-0.8+(-3.8+0.8)/2, "A", cex=2, adj=c(0.5,0.5), font=2)

        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT_CATCH[1:2]+c(-0.41,-0.41), ylim=EXT_CATCH[3:4])
        for(i in seq(-27,-21,2)) {
            segments(EXT_CATCH[1], i, EXT_CATCH[2], i, lwd=0.3, col='grey50')
            text(EXT_CATCH[1]-0.1, i, i, adj=c(1,0.4), cex=0.8)
        }
        for(i in seq(25,35,2)) {
            text(i, EXT_CATCH[3]-0.1, i, adj=c(0,1), cex=0.8)
        }
        image(mask(BIO1, WATERSHED), col=RdYlBu, zlim=c(8,32), breaks=seq(8,32,2), add=TRUE, legend=FALSE, interpolate=FALSE)
        image(mask(BIO1, WATERSHED2), col=RdYlBu, zlim=c(8,32), breaks=seq(8,32,2), add=TRUE, legend=FALSE, interpolate=FALSE)
        plot(WATERSHED, add=TRUE, col=NA, border='black')
        plot(WATERSHED2, add=TRUE, col=NA, border='black')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=3, lwd=1.5)
        text(32.6, -26.9, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        rect(EXT_CATCH[1], EXT_CATCH[3], EXT_CATCH[2], EXT_CATCH[4]+20, col=NA,border='black', lwd=0.3)

        ##-- PANE B
        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT[1:2], ylim=EXT[3:4])
        for(i in seq(-35,5,5)) {
            if(i%%10 == 0) segments(EXT[1], i, EXT[2], i, lwd=0.3, col='grey50')
        }
        for(i in seq(10,40,5)) {
            text(i, EXT[4]+0.3, i, adj=c(0,0), cex=0.8)
        }
        plot(VEG,add=TRUE, lwd=0.1, col=DIV.10col[VEG$BIOME], border=NA)
        plot(LAKES, add=TRUE, col="grey50", border=NA)
        plot(M1, lwd=0.5, add=TRUE)
        plot(WATERSHED, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='black')
        plot(WATERSHED2, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='black')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=1.5)
        text(33, -27.5, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        points(34+26/60, -11-18/60,pch=23, col='deeppink4', bg='deeppink', cex=2, lwd=1.2)
        text(35.5, -11-18/60, 'Lake\nMalawi', cex=1, adj=c(0, 0.6), font=2)
        rect(EXT[1], EXT[3]-20, EXT[2], EXT[4], col=NA,border='black', lwd=0.3)
        LABELS=c('1', '2', '7', '9', '10', '12', '13')
        for(i in seq(1,7,1)){
            rect(EXT[1]+5.5+(i-1)*3, -37.5, EXT[1]+5.5+(i-1)*3+2.5, -38.5, col=DIV.10col[as.numeric(LABELS[i])], border='black', lwd=0.4)
            text(EXT[1]+5.5+(i-1)*3+1.25, -38.9, LABELS[i], adj=c(0.5,1), cex=0.8)
        }
        text(EXT[1]+(EXT[2]-EXT[1])/2, -36.5, 'Terrestrial biomes (cf. caption)', adj=c(0.5, 0.5), cex=0.9)
        rect(EXT[1]+0.8, EXT[4]-0.8, EXT[1]+3.5, EXT[4]-3.8, border="black", lwd=0.4, col="white")
        text(EXT[1]+0.8+(3.5-0.8)/2, EXT[4]-0.8+(-3.8+0.8)/2, "B", cex=2, adj=c(0.5,0.5), font=2)

        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT_CATCH[1:2], ylim=EXT_CATCH[3:4])
        for(i in seq(-27,-21,2)) {
            segments(EXT_CATCH[1], i, EXT_CATCH[2], i, lwd=0.3, col='grey50')
        }
        for(i in seq(25,35,2)) {
            text(i, EXT_CATCH[3]-0.1, i, adj=c(0,1), cex=0.8)
        }
        plot(raster::intersect(VEG, WATERSHED),add=TRUE, lwd=0.1, col=DIV.10col[raster::intersect(VEG, WATERSHED)$BIOME], border=NA)
        plot(raster::intersect(VEG, WATERSHED2[1,]),add=TRUE, lwd=0.1, col=DIV.10col[raster::intersect(VEG, WATERSHED2[1,])$BIOME], border=NA)
        plot(crop(VEG, extent(WATERSHED2[2,])),add=TRUE, lwd=0.1, col=DIV.10col[crop(VEG, extent(WATERSHED2[2,]))$BIOME], border=NA)
        plot(raster::intersect(VEG, WATERSHED2[3,]),add=TRUE, lwd=0.1, col=DIV.10col[raster::intersect(VEG, WATERSHED2[3,])$BIOME], border=NA)
        plot(WATERSHED, add=TRUE, col=NA, border='black')
        plot(WATERSHED2, add=TRUE, col=NA, border='black')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=3, lwd=1.5)
        text(32.6, -26.9, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        rect(EXT_CATCH[1], EXT_CATCH[3], EXT_CATCH[2], EXT_CATCH[4]+20, col=NA,border='black', lwd=0.3)

        ##-- PANE C
        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT[1:2]+c(1.2,1.2), ylim=EXT[3:4])
        for(i in seq(-35,5,5)) {
            if(i%%10 == 0) segments(EXT[1], i, EXT[2], i, lwd=0.3, col='grey50')
            text(EXT[2]+0.3, i, i, adj=c(0,0.4), cex=0.8)
        }
        for(i in seq(10,40,5)) {
            text(i, EXT[4]+0.3, i, adj=c(0,0), cex=0.8)
        }
        plot(CONTINENT, lwd=0.5, add=TRUE, col='black')
        image(GBIF, col=heatcol4, zlim=c(1,4), add=TRUE, legend=FALSE, interpolate=FALSE)
        plot(LAKES, add=TRUE, col="grey50", border=NA)
        plot(M1, lwd=0.5, add=TRUE)
        plot(WATERSHED, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.7), border='black')
        plot(WATERSHED2, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='black')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=1.5)
        text(33, -27.5, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        points(34+26/60, -11-18/60,pch=23, col='deeppink4', bg='deeppink', cex=2, lwd=1.2)
        text(35.5, -11-18/60, 'Lake\nMalawi', cex=1, adj=c(0, 0.6), font=2, col='white')
        rect(EXT[1], EXT[3]-20, EXT[2], EXT[4], col=NA,border='black', lwd=0.3)
        ##
        LABELS=c('1-9', '10-99', '100-999', '>1000')
        for(i in seq(1,4)){
            rect(EXT[1]+7.5+(i-1)*4.5, -37.5, EXT[1]+7.5+(i-1)*4.5+3.5, -38.5, col=heatcol4[i], border='black', lwd=0.4)
            text(EXT[1]+7.5+(i-1)*4.5+1.75, -38.9, LABELS[i], adj=c(0.5,1), cex=0.8)
        }
        text(EXT[1]+(EXT[2]-EXT[1])/2, -36.5, 'Number of plant occurence per grid cell', adj=c(0.5, 0.5), cex=0.9)
        rect(EXT[1]+0.8, EXT[4]-0.8, EXT[1]+3.5, EXT[4]-3.8, border="black", lwd=0.4, col="white")
        text(EXT[1]+0.8+(3.5-0.8)/2, EXT[4]-0.8+(-3.8+0.8)/2, "C", cex=2, adj=c(0.5,0.5), font=2)

        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT_CATCH[1:2]+c(0.41,0.41), ylim=EXT_CATCH[3:4])
        for(i in seq(-27,-21,2)) {
            segments(EXT_CATCH[1], i, EXT_CATCH[2], i, lwd=0.3, col='grey50')
            text(EXT_CATCH[2]+0.1, i, i, adj=c(0,0.4), cex=0.8)
        }
        for(i in seq(25,35,2)) {
            text(i, EXT_CATCH[3]-0.1, i, adj=c(0,1), cex=0.8)
        }
        plot(WATERSHED, add=TRUE, col='black', border='black')
        image(mask(GBIF, WATERSHED), col=heatcol4, zlim=c(1,4), add=TRUE, legend=FALSE, interpolate=FALSE)
        image(mask(GBIF, WATERSHED2), col=heatcol4, zlim=c(1,4), add=TRUE, legend=FALSE, interpolate=FALSE)
        plot(WATERSHED, add=TRUE, col=NA, border='black')
        plot(WATERSHED2, add=TRUE, col=NA, border='black')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=3, lwd=1.5)
        text(32.6, -26.9, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        rect(EXT_CATCH[1], EXT_CATCH[3], EXT_CATCH[2], EXT_CATCH[4]+20, col=NA,border='black', lwd=0.3)

    }  ;  dev.off()
}

#-;
