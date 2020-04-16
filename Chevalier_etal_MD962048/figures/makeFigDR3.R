## Figure 4: Comparison of the reconstruction with independent records
##

OUTPUT_FOLDER=getwd()
s <- readline(prompt=paste0("Where should the figure be saved?\nDefault is current workin directory (",OUTPUT_FOLDER,"): "))
if(s != '') OUTPUT_FOLDER <- s

pkg2install=c()
if (! ("rio" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'rio')

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

    cat(">>> Loading data.\n")
    MAT=rio::import('https://github.com/mchevalier2/ClimateReconstructions/raw/master/MD96-2048_MAT_01.xlsx', which=2)[1:181,]
    CO2=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=7)[62:1883,c(2,3)]
    SSTs=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=5)[1:306, c(23, 24)]
    LR04=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=10)[1:701,c(1,2)]
    DomeC=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=6)[1:5784,c(1,2)]
    MALAWI=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=8)[1:295,c(1,3)]
    LEAFWAX=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=9)[1:177,c(1,10)]
    LEAFWAX.detrended=cbind(LEAFWAX[,1], LEAFWAX[,2] - LEAFWAX[,1]*coef(lm(LEAFWAX[,2]~ LEAFWAX[,1]))[2] - coef(lm(LEAFWAX[,2]~ LEAFWAX[,1]))[1])

    COL='black'
    COL2=rgb(217,95,2,maxColorValue=255)

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_FigDR3.pdf"), width=7.48, height=7, useDingbats=FALSE)  ;  {
        par(ps=7,bg=makeTransparent("white",alpha=0),mar=rep(0,4),cex=1,cex.main=1)
        layout(matrix(1:6, ncol=3, byrow=TRUE), width=1, height=1)

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(10, max(MAT[,2]), 'A', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(SSTs[,2])+diff(range(SSTs[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(SSTs, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(-3,3,1)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(-3,3,2)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(SSTs[,2])+diff(range(SSTs[,2]))/2, 'Mozambique Channel SSTs PC1', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## LEAFWAX.detrended
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(10, max(MAT[,2]), 'B', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=rev(range(LEAFWAX.detrended[,2]))-diff(range(LEAFWAX.detrended[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(LEAFWAX.detrended, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(-0.04,0.05,0.01)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(-0.04,0.05,0.02)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(LEAFWAX.detrended[,2])+diff(range(LEAFWAX.detrended[,2]))/2, 'Ratio of long-chain n-alkanes C31/(C29+C31) [detrended]', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## MALAWI
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(10, max(MAT[,2]), 'C', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(MALAWI[,2])+diff(range(MALAWI[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MALAWI, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(17,27,1)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(17,27,2)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(MALAWI[,2])+diff(range(MALAWI[,2]))/2, 'Malawi lake surface temperature (°C)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## DomeC
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(10, max(MAT[,2]), 'D', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(DomeC[,2])+diff(range(DomeC[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(DomeC, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(-10,5,2.5)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(-10,5,5)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(DomeC[,2])+diff(range(DomeC[,2]))/2, 'Dome C, Antarctica Temperature Reconstruction (°C)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## CO2
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(10, max(MAT[,2]), 'E', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(CO2[,2])+diff(range(CO2[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(CO2, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(180,315,15)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(180,315,30)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(CO2[,2])+diff(range(CO2[,2]))/2, 'Dome C, Antarctica pCO2 (ppm)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## LR04
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent(COL, alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(10, max(MAT[,2]), 'F', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=rev(range(LR04[,2]))-diff(range(LR04[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(LR04, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(5,3.2,-0.25)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(5,3.2,-0.5)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(LR04[,2])+diff(range(LR04[,2]))/2, 'Global ice volume (LR04) d18Obenthic (permil VPDB)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }
    dev.off()  ;  }

}


#-;
