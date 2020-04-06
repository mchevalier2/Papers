## Figure 2: Study area and calibration dataset
##
## Loading necessary data
load(url('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/figures/Fig2.RData'))

OUTPUT_FOLDER=getwd()
s <- readline(prompt=paste0("Where should the figure be saved?\nDefault is current workin directory (",OUTPUT_FOLDER,"): "))
if(s != '') OUTPUT_FOLDER <- s

pkg2install=c()
if (! ("sp" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'sp')
if (! ("raster" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'raster')
if (! ("shape" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'shape')


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
    EXT <- c(10,42,-40,5)

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_Fig2.pdf"), width=7.54, height=7.54/3*(EXT[4]-EXT[3])/(EXT[2]-EXT[1]))  ;  {
        par(mar=c(0,0,0,0), mfrow=c(1,3))

        ##-- PANE A
        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT[1:2], ylim=EXT[3:4])
        image(BIO1, col=RdYlBu, zlim=c(8,32), breaks=seq(8,32,2), add=TRUE, legend=FALSE, interpolate=FALSE)
        plot(LAKES, add=TRUE, col="grey50", border=NA)
        plot(M1, lwd=0.5, add=TRUE)
        plot(WATERSHED, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='aquamarine4')
        plot(WATERSHED2, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='aquamarine4')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=1.5)
        text(33, -27.5, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        rect(EXT[1], EXT[3], EXT[2], EXT[4], col=NA,border='black', lwd=0.3)
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

        ##-- PANE B
        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT[1:2], ylim=EXT[3:4])
        plot(VEG,add=TRUE, lwd=0.1, col=DIV.10col[VEG$BIOME], border=NA)
        plot(LAKES, add=TRUE, col="grey50", border=NA)
        plot(M1, lwd=0.5, add=TRUE)
        plot(WATERSHED, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='aquamarine4')
        plot(WATERSHED2, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='aquamarine4')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=1.5)
        text(33, -27.5, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        rect(EXT[1], EXT[3], EXT[2], EXT[4], col=NA,border='black', lwd=0.3)

        LABELS=c('1', '2', '7', '9', '10', '12', '13')
        for(i in seq(1,7,1)){
            rect(EXT[1]+5.5+(i-1)*3, -37.5, EXT[1]+5.5+(i-1)*3+2.5, -38.5, col=DIV.10col[as.numeric(LABELS[i])], border='black', lwd=0.4)
            text(EXT[1]+5.5+(i-1)*3+1.25, -38.9, LABELS[i], adj=c(0.5,1), cex=0.8)
        }
        text(EXT[1]+(EXT[2]-EXT[1])/2, -36.5, 'Terrestrial biomes (cf. caption)', adj=c(0.5, 0.5), cex=0.9)
        rect(EXT[1]+0.8, EXT[4]-0.8, EXT[1]+3.5, EXT[4]-3.8, border="black", lwd=0.4, col="white")
        text(EXT[1]+0.8+(3.5-0.8)/2, EXT[4]-0.8+(-3.8+0.8)/2, "B", cex=2, adj=c(0.5,0.5), font=2)

        ##-- PANE C
        plot(0,0, type='n', axes=FALSE,frame=FALSE, asp=1, xlim=EXT[1:2], ylim=EXT[3:4])
        plot(CONTINENT, lwd=0.5, add=TRUE, col='black')
        plot(GBIF, col=heatcol4, zlim=c(1,4), add=TRUE, legend=FALSE, interpolate=FALSE)
        plot(LAKES, add=TRUE, col="grey50", border=NA)
        plot(M1, lwd=0.5, add=TRUE)
        plot(WATERSHED, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.7), border='aquamarine4')
        plot(WATERSHED2, add=TRUE, col=makeTransparent('aquamarine1', alpha=0.5), border='aquamarine4')
        plot(LIMPOPO, add=TRUE, col='aquamarine4',lwd=((LIMPOPO$NAME=='Limpopo')+1))
        plot(LIMPOPO2, add=TRUE, col='aquamarine4',lwd=((LIMPOPO2$NAME=='Limpopo')+1))
        points(34.016700, -26.166700,pch=23, col='aquamarine4', bg='aquamarine1', cex=1.5)
        text(33, -27.5, 'MD96-2048', cex=1, adj=c(0, 1), font=2)
        rect(EXT[1], EXT[3], EXT[2], EXT[4], col=NA,border='black', lwd=0.3)
        ##
        LABELS=c('1-9', '10-99', '100-999', '>1000')
        for(i in seq(1,4)){
            rect(EXT[1]+7.5+(i-1)*4.5, -37.5, EXT[1]+7.5+(i-1)*4.5+3.5, -38.5, col=heatcol4[i], border='black', lwd=0.4)
            text(EXT[1]+7.5+(i-1)*4.5+1.75, -38.9, LABELS[i], adj=c(0.5,1), cex=0.8)
        }
        text(EXT[1]+(EXT[2]-EXT[1])/2, -36.5, 'Number of plant occurence per grid cell', adj=c(0.5, 0.5), cex=0.9)
        rect(EXT[1]+0.8, EXT[4]-0.8, EXT[1]+3.5, EXT[4]-3.8, border="black", lwd=0.4, col="white")
        text(EXT[1]+0.8+(3.5-0.8)/2, EXT[4]-0.8+(-3.8+0.8)/2, "C", cex=2, adj=c(0.5,0.5), font=2)
    }  ;  dev.off()
}

#-;
