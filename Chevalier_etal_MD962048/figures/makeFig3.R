## Figure 3: Presentation of the reconstruction
##
## Loading necessary data
load(url('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/figures/Fig3.RData'))

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
    CREST_folder='/Users/mchevali1/Research/CREST/Sites/MD96-2048/MD96-2048_2019_18.06.19_at_16h20_1'

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

    XX.interp=1:800
    PANN=read.csv(paste0(CREST_folder, "/Reconstructions/bio1.csv"))[1:181,]
    PANN.ysmooth=gausmooth(PANN[,c(1,2)], XX.interp, mean(diff(PANN[,1])))
    pdf=read.csv(paste0(CREST_folder, "/Densities/Numerical_Values_PdfVar/bio1.csv"))[,1:182]
    pdfter=pdf

    for(i in 2:ncol(pdf)){
        oo=rev(order(pdf[,i]))
        tmp1=pdf[oo,1]
        tmp2=pdf[oo,i]
        oo=order(tmp1)
        pdfter[,i]=cumsum(tmp2/sum(tmp2))[oo]
    }

    MORLET=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/Data_790kyr.xlsx', which=7)
    MORLET.log2=MORLET[,1]
    MORLET.sig=unique(MORLET[,2])
    MORLET.x=1:396
    MORLET=log2(t(MORLET[,-c(1,2)])[,106:1])
    MORLET[MORLET < -4] = -4
    wave= dplR::morlet(gausmooth(PANN[,1:2], seq(1,792,2), mean(diff(PANN[,1]))), x1=1:396, siglvl=0.95)

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_Fig3.pdf"), width=7.54, height=7.54/2, useDingbats=FALSE)  ;  {
        par(mar=c(2.3,2.2,3,0.5))
        layout(matrix(1:2, ncol=2), width=c(4,2))
        plot3D::image2D(z=(1-as.matrix(t(pdfter[,-1]))),y=pdfter[,1], x=PANN[,1], xlim=c(0,790), ylim=c(15,24), zlim=c(0,1), col = plot3D::gg2.col(200)[1:100], cex.axis=6/7, colkey=FALSE, resfac=2, tck=-.013, mgp=c(1.3, .3, 0), las=1, hadj=c(1,1), xlab='Age (calendar yr BP x1000)', ylab='Mean Annual Temperature (°C)', cex.lab=6/7)
        segments(430,15,430,24, lty=2)
        text(425,23.8, 'MBT', cex=6/7, adj=c(1,0), srt=90)
        text(435,15.2, 'MBT', cex=6/7, adj=c(0,1), srt=90)
        points(PANN[,1:2], pch=18, col='white', cex=0.8)
        points(PANN[,1:2], col='white', cex=0.3, type='l')
        points(1:800, PANN.ysmooth, pch=15, col='black', cex=0.3, type='l')
        plot3D::colkey(side=3, length=0.8, dist=-0.01, lwd=0.1, cex.axis=6/7, clim=c(1,0), col=plot3D::gg2.col(200)[1:100], clab='A - Confidence level', font.clab=1, line.clab=1.3, adj.clab=0.5, add=TRUE, tck=-0.4, mgp=c(3, .25, 0), lwd.tick=0.7)

        par(mar=c(2.3,2.2,3,.2))
        plot3D::image2D(z=MORLET,y=2**(rev(MORLET.log2)), x=MORLET.x*2, ylim=rev(range(2**(rev(MORLET.log2)))), col = plot3D::jet.col(100), cex.axis=6/7, colkey=FALSE, resfac=2, tck=-.013, mgp=c(1.3, .3, 0), las=1, hadj=c(1,1), xlab='Age (calendar yr BP x1000)', ylab='Periods (in thousand of years)', cex.lab=6/7, contour=list(levels=log2(MORLET.sig), drawlabels=FALSE), log='y', lwd=1.5)
        polygon(c(0,wave$x*2, 792,0),c(396,2*2**log2(wave$coi), 396,396), col=makeTransparent('white', alpha=0.6), lwd=0.2)
        plot3D::colkey(side=3, length=0.8, dist=-0.01, lwd=0.1, cex.axis=6/7, clim=range(MORLET), col=plot3D::jet.col(100), clab='B - log2(power)', font.clab=1, line.clab=1.3, adj.clab=0.5, add=TRUE, tck=-0.4, mgp=c(3, .25, 0), lwd.tick=0.7)
    }  ; dev.off()



}


#-;