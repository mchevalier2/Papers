## Figure 3: Presentation of the reconstruction
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

    XX.interp=0:800
    MAT=rio::import('https://github.com/mchevalier2/ClimateReconstructions/raw/master/MD96-2048_MAT_01.xlsx', which=2)[1:181,]
    pdf=rio::import('https://github.com/mchevalier2/ClimateReconstructions/raw/master/MD96-2048_MAT_01.xlsx', which=3)
    colnames(pdf) = pdf[1,]
    pdf = pdf[-1,]
    MAT.ysmooth=gausmooth(MAT[,c(1,2)], XX.interp, mean(diff(MAT[,1])))
    pdfter=pdf

    for(i in 2:ncol(pdf)){
        oo=rev(order(pdf[,i]))
        tmp1=pdf[oo,1]
        tmp2=pdf[oo,i]
        oo=order(tmp1)
        pdfter[,i]=cumsum(tmp2/sum(tmp2))[oo]
    }


    MAT.interp=approx(MAT[,1:2], xout=seq(0,790,1))
    morlet=dplR::morlet(MAT.interp$y, MAT.interp$x, siglvl=0.99, p2=8.7, dj=0.1)
    morletP=log2(morlet$Power)[,ncol(morlet$Power):1]
    morletP[morletP < -4] = -4

    Signif <- t(matrix(morlet$Signif, dim(morlet$Power)[2], dim(morlet$Power)[1]))
    Signif <- morlet$Power/Signif

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_Fig3.pdf"), width=7.54, height=7.54/2, useDingbats=FALSE)  ;  {
        par(mar=c(2.3,2.2,3,0.5))
        layout(matrix(1:2, ncol=2), width=c(4,2))
        plot3D::image2D(z=(1-as.matrix(t(pdfter[,-1]))),y=pdfter[,1], x=MAT[,1], xlim=c(0,790), ylim=c(15,24), zlim=c(0,1), col = plot3D::gg2.col(200)[1:100], cex.axis=6/7, colkey=FALSE, resfac=2, tck=-.013, mgp=c(1.3, .3, 0), las=1, hadj=c(1,1), xlab='Age (calendar yr BP x1000)', ylab='Mean Annual Temperature (°C)', cex.lab=6/7)
        segments(430,15,430,24, lty=2)
        text(425,23.8, 'MBT', cex=6/7, adj=c(1,0), srt=90)
        text(435,15.2, 'MBT', cex=6/7, adj=c(0,1), srt=90)
        points(MAT[,1:2], pch=18, col='white', cex=0.8)
        points(MAT[,1:2], col='white', cex=0.3, type='l')
        points(XX.interp, MAT.ysmooth, pch=15, col='black', cex=0.3, type='l')
        plot3D::colkey(side=3, length=0.8, dist=-0.01, lwd=0.1, cex.axis=6/7, clim=c(1,0), col=plot3D::gg2.col(200)[1:100], clab='A - Confidence level', font.clab=1, line.clab=1.3, adj.clab=0.5, add=TRUE, tck=-0.4, mgp=c(3, .25, 0), lwd.tick=0.7)

        par(mar=c(2.3,2.2,3,.2))
        plot3D::image2D(z=morletP[,1:65],y=rev(morlet$period)[1:65], x=morlet$x, ylim=rev(range(rev(morlet$period)[1:65])), col = plot3D::jet.col(100), cex.axis=6/7, colkey=FALSE, resfac=2, tck=-.013, mgp=c(1.3, .3, 0), las=1, hadj=c(1,1), xlab='Age (calendar yr BP x1000)', ylab='Periods (in thousand of years)', cex.lab=6/7, contour=FALSE, log='y', lwd=1.5)
        contour(morlet$x, morlet$period, Signif, levels = 1, labels = morlet$siglvl, drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, add = TRUE, lwd = 1, col = "black")
        polygon(c(0,morlet$x, 792,0),c(max(morlet$Scale),2**log2(morlet$coi), max(morlet$period),max(morlet$period)), col=makeTransparent('white', alpha=0.6), lwd=0.2)
        plot3D::colkey(side=3, length=0.8, dist=-0.01, lwd=0.1, cex.axis=6/7, clim=range(morletP), col=plot3D::jet.col(100), clab='B - log2(power)', font.clab=1, line.clab=1.3, adj.clab=0.5, add=TRUE, tck=-0.4, mgp=c(3, .25, 0), lwd.tick=0.7)
    }  ; dev.off()
}

#-;
