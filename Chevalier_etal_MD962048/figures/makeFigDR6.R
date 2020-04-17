## Figure DR6: Morlet nalysis of pollen diversity
##

OUTPUT_FOLDER=getwd()
s <- readline(prompt=paste0("Where should the figure be saved?\nDefault is current workin directory (",OUTPUT_FOLDER,"): "))
if(s != '') OUTPUT_FOLDER <- s

pkg2install=c()
if (! ("rio" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'rio')
if (! ("dplR" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'dplR')
if (! ("plot3D" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'plot3D')


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

    ## Calculate the Gaussian density of Probability
    ## defined by xbar and sigma, at x
    gauss=function(x, xbar, sigma){return(1/sqrt(2*pi*sigma**2)*exp(-(x-xbar)**2/sigma**2))}

    ## Apply the Gaussian smoothing kernel on dat, using sigma
    ## as a kernel width. xout defines the output axis
    gausmooth=function(dat, xout, sigma, interp=TRUE){
        yout=rep(NA,length(xout))
        for(i in 1:length(xout)){
          if((xout[i] >= min(dat[,1]) & xout[i] <= max(dat[,1])) | interp){
            yout[i]=sum(dat[,2]*gauss(dat[,1], xout[i], sigma))/sum(gauss(dat[,1], xout[i], sigma))
          }
        }
        return(yout)
    }

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

    MAT=rio::import('https://github.com/mchevalier2/ClimateReconstructions/raw/master/MD96-2048_MAT_01.xlsx', which=2)[1:181,]
    POLLEN=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,-c(2,3)]
    POLLENSUM=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,c(1,2,3)]
    ECC=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=2)[1:800,c(1,2)]

    ## Species richness (S)
    S <- vegan::specnumber(POLLEN[,-1]) ## rowSums(BCI > 0) does the same... # Richness
    ## Pielou's evenness
    DMG=(S-1) / log(POLLENSUM[,2])


    DMG.interp=approx(MAT[,1],DMG, xout=seq(0,790,1))
    morlet=dplR::morlet(DMG.interp$y, DMG.interp$x, siglvl=0.99, p2=8.7, dj=0.1)
    morletP=log2(morlet$Power)[,ncol(morlet$Power):1]
    morletP[morletP < -4] = -4

    Signif <- t(matrix(morlet$Signif, dim(morlet$Power)[2], dim(morlet$Power)[1]))
    Signif <- morlet$Power/Signif

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_FigDR6.pdf"), width=7.54, height=7.54/2, useDingbats=FALSE)  ;  {
        par(ps=7,bg=makeTransparent("white",alpha=0),mar=rep(0,4),cex=1,cex.main=1)
        layout(matrix(1:2, ncol=2, byrow=TRUE), width=c(1.2, 0.8), height=1)

        COL='black'
        COL2='red'
        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(DMG)+diff(range(DMG))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], DMG, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(DMG)+diff(range(DMG))/2, "Margalef's Index", adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(DMG)-0.02*diff(range(DMG)),809,max(DMG)+0.02*diff(range(DMG)),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(DMG)-0.02*diff(range(DMG)),i,min(DMG)-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(DMG)), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(DMG)-0.04*diff(range(DMG)), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(DMG)-0.1*diff(range(DMG)), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
            }
            plot.window(xlim=c(-100,900),ylim=range(ECC[,2])+diff(range(ECC[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(ECC, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(0.005,0.05,0.005)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(0.005,0.05,0.01)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(ECC[,2])+diff(range(ECC[,2]))/2, 'Eccentricty', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        par(mar=c(3.5,2.2,3,.2))
        plot3D::image2D(z=morletP[,1:70],y=rev(morlet$period)[1:70], x=morlet$x, ylim=rev(range(rev(morlet$period)[1:70])), col = plot3D::jet.col(100), cex.axis=7/7, colkey=FALSE, resfac=2, tck=-.013, mgp=c(1.3, .3, 0), las=1, hadj=c(1,1), xlab='Age (calendar yr BP x1000)', ylab='Periods (in thousand of years)', cex.lab=8/7, contour=FALSE, log='y', lwd=1.5)
        contour(morlet$x, morlet$period, Signif, levels = 1, labels = morlet$siglvl, drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, add = TRUE, lwd = 1, col = "black")
        polygon(c(0,morlet$x, 792,0),c(max(morlet$Scale),2**log2(morlet$coi), max(morlet$period),max(morlet$period)), col=makeTransparent('white', alpha=0.6), lwd=0.2)
        plot3D::colkey(side=3, length=0.8, dist=-0.01, lwd=0.1, cex.axis=8/7, clim=range(morletP), col=plot3D::jet.col(100), clab='log2(power)', font.clab=1, line.clab=1.3, adj.clab=0.5, add=TRUE, tck=-0.4, mgp=c(3, .25, 0), lwd.tick=0.7)
    }  ; dev.off()
}






#-;
