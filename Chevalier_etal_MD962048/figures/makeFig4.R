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
    INSOL=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=1)[1201:2001,c(3,7)]
    ECC=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=2)[1:800,c(1,2)]
    CO2=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=7)[62:1883,c(2,3)]
    SSTs=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=5)[1:306, c(23, 24)]
    LR04=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=10)[1:701,c(1,2)]
    POLLENSUM=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,c(1,2,3)]
    POLLENTAXA=apply(ifelse(rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,-c(1,2,3,219,220)]>0,1,0),1,sum)
    POLLEN=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,-c(2,3)]
    MARGALEF=(apply(POLLEN[,-1] > 0, 1, sum) - 1) / log(POLLENSUM[,2])


    XX.interp=MAT[,1]
    MAT.smooth=gausmooth(MAT[,c(1,2)], XX.interp, mean(diff(XX.interp)))
    CO2.smooth=gausmooth(CO2, XX.interp, mean(diff(XX.interp)))
    SSTs.smooth=gausmooth(SSTs, XX.interp, mean(diff(XX.interp)))
    LR04.smooth=gausmooth(LR04, XX.interp, mean(diff(XX.interp)))


    TIME=c()
    for(i in 1:180){
      TIME=c(TIME, (POLLENSUM[i,1]+POLLENSUM[i+1,1])/2)
    }
    TIME=c(0, TIME, 790)


    MIS.col=makeTransparent(c(rgb(205,91,65,maxColorValue=255),rgb(255,127,80,maxColorValue=255),rgb(122,197,205,maxColorValue=255),rgb(152,245,255,maxColorValue=255), 'white'), alpha=0.2)

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_Fig4.pdf"), width=3.54, height=9.45, useDingbats=FALSE)  ;  {
      par(ps=7,bg=makeTransparent("white",alpha=0),mar=rep(0,4),cex=1,cex.main=1)
      plot.new()

      plot.window(xlim=c(-150,950),ylim=c(-0.0125,1.0125),main='',ylab='',xlab='') ;  {
          ## Top part
          for(i in c(11.7, 29,57,71,130,191,243,300,337,374,424,478,533,565,621,676,712,761,790)){  segments(i,0,i,1.01, col="grey50", lwd=0.5, lty=2)  }
          rect(-9,0,11.7,1, col=MIS.col[5], border=NA)  ;  rect(-9,1.02,11.7,1, col="white", border='black', lwd=0.5)  ;  text(1.375,1.01,'1',cex=0.8, adj=c(0.5,0.5))
          rect(11.7,0,29,1, col=MIS.col[4], border=NA)  ;  rect(29,1.02,11.7,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(20.35,1.01,'2',cex=0.8, adj=c(0.5,0.5))
          rect(29,0,57,1, col=MIS.col[3], border=NA)  ;  rect(29,1.02,57,1, col=MIS.col[3], border='black', lwd=0.5)  ;  text(43,1.01,'3',cex=0.8, adj=c(0.5,0.5))
          rect(57,0,71,1, col=MIS.col[4], border=NA)  ;  rect(57,1.02,71,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(64,1.01,'4',cex=0.8, adj=c(0.5,0.5))
          rect(71,0,130,1, col=MIS.col[5], border=NA)  ;  rect(71,1.02,130,1, col="white", border='black', lwd=0.5)  ;  text(100.5,1.01,'5',cex=0.8, adj=c(0.5,0.5))
          rect(130,0,191,1, col=MIS.col[4], border=NA)  ;  rect(130,1.02,191,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(160.5,1.01,'6',cex=0.8, adj=c(0.5,0.5))
          rect(191,0,243,1, col=MIS.col[5], border=NA)  ;  rect(191,1.02,243,1, col="white", border='black', lwd=0.5)  ;  text(217,1.01,'7',cex=0.8, adj=c(0.5,0.5))
          rect(243,0,300,1, col=MIS.col[4], border=NA)  ;  rect(243,1.02,300,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(271.5,1.01,'8',cex=0.8, adj=c(0.5,0.5))
          rect(300,0,337,1, col=MIS.col[5], border=NA)  ;  rect(300,1.02,337,1, col="white", border='black', lwd=0.5)  ;  text(318.5,1.01,'9',cex=0.8, adj=c(0.5,0.5))
          rect(337,0,374,1, col=MIS.col[4], border=NA)  ;  rect(337,1.02,374,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(355.5,1.01,'10',cex=0.8, adj=c(0.5,0.5))
          rect(374,0,424,1, col=MIS.col[5], border=NA)  ;  rect(374,1.02,424,1, col="white", border='black', lwd=0.5)  ;  text(399,1.01,'11',cex=0.8, adj=c(0.5,0.5))
          rect(424,0,478,1, col=MIS.col[4], border=NA)  ;  rect(424,1.02,478,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(451,1.01,'12',cex=0.8, adj=c(0.5,0.5))
          rect(478,0,533,1, col=MIS.col[5], border=NA)  ;  rect(478,1.02,533,1, col="white", border='black', lwd=0.5)  ;  text(505.5,1.01,'13',cex=0.8, adj=c(0.5,0.5))
          rect(533,0,563,1, col=MIS.col[4], border=NA)  ;  rect(533,1.02,563,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(548,1.01,'14',cex=0.8, adj=c(0.5,0.5))
          rect(563,0,621,1, col=MIS.col[5], border=NA)  ;  rect(563,1.02,621,1, col="white", border='black', lwd=0.5)  ;  text(592,1.01,'15',cex=0.8, adj=c(0.5,0.5))
          rect(621,0,676,1, col=MIS.col[4], border=NA)  ;  rect(621,1.02,676,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(648.5,1.01,'16',cex=0.8, adj=c(0.5,0.5))
          rect(676,0,712,1, col=MIS.col[5], border=NA)  ;  rect(676,1.02,712,1, col="white", border='black', lwd=0.5)  ;  text(694,1.01,'17',cex=0.8, adj=c(0.5,0.5))
          rect(712,0,761,1, col=MIS.col[4], border=NA)  ;  rect(712,1.02,761,1, col=MIS.col[4], border='black', lwd=0.5)  ;  text(736.5,1.01,'18',cex=0.8, adj=c(0.5,0.5))
          rect(761,0,790,1, col=MIS.col[5], border=NA)  ;  rect(761,1.02,790,1, col="white", border='black', lwd=0.5)  ;  text(775.5,1.01,'19',cex=0.8, adj=c(0.5,0.5))
          rect(790,0,809,1, col=MIS.col[4], border=NA)  ;  rect(790,1.02,809,1, col=MIS.col[4], border='black', lwd=0.5)  #;  text(790,1.015,'20',cex=0.8, adj=c(0.5,0.5))

          ## Bottom part
          rect(-9,0,809,1, col=NA, border="black", lwd=0.5)
          text(400,1.035, "Marine Isotope Stages", cex=8/7, adj=c(0.5,0.5))
          for(i in seq(0,800,25)){  segments(i,0,i,-0.004-ifelse(i%%50 == 0, 0.002,0), lwd=0.5)  }
          for(i in seq(0,800,100)){  text(i,-0.01, i, cex=1, adj=c(0.5, 1))  }
          text(400, -0.03, 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
      }

      ## Insolation
      plot.window(xlim=c(-150,950),ylim=c(425,515)+90*c(-5.7,0)/0.7+90*c(-6.4*0.025, 6.4*0.025)/0.7,main='',ylab='',xlab='')  ;  {
          COL='black'
          points(INSOL, col='grey70', type='l', lwd=1.2)
          for(i in seq(430,510,10)) segments(-20,i,-9,i, lwd=0.5, col=COL)
          for(i in seq(430,510,20)) text(-25,i,i, adj=c(1,0.5), col=COL)
          text(-175, 465, '(a) DJF Insolation\nat 20°S (W/m2)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
      }

      ## Eccentricity
      plot.window(xlim=c(-150,950),ylim=c(-0.04,0.06)+0.1*c(-6,0)/0.7+0.1*c(-6.4*0.025, 6.4*0.025)/0.7,main='',ylab='',xlab='')  ;  {
          COL='black'
          points(ECC, col='grey20', type='l', lwd=1.2)
          for(i in seq(0.004, 0.04,0.008)) segments(820,i,809,i, lwd=0.5, col=COL)
          for(i in seq(0.004, 0.04,0.016)) text(825,i,i, adj=c(0,0.5), col=COL)
          text(975, 0.022, '(b) Eccentricity', adj=c(0.5,0), srt=90, col=COL, cex=8/7)
      }

      ## Pollen counts
      plot.window(xlim=c(-150,950),ylim=c(0,450)+450*c(-5,0.7)/0.7+450*c(-6.4*0.025, 6.4*0.025)/0.7,main='',ylab='',xlab='')  ;  {
          COL=makeTransparent('cornflowerblue', alpha=0.9)
          polygon(c(0,POLLENSUM[,1],POLLENSUM[181,1]), c(0,POLLENSUM[,3],0), col=COL, border=NA)
          COL=makeTransparent('forestgreen', alpha=0.9)
          polygon(c(0,POLLENSUM[,1],POLLENSUM[181,1]), c(0,POLLENSUM[,2],0), col=COL, border=NA)
          for(i in seq(0,400,50)) segments(820,i,809,i, lwd=0.5, col='black')
          for(i in seq(0,400,100)) text(825,i,i, adj=c(0,0.5), col='black')
          text(975, 225, '(c) Pollen counts\nTerrestrial & Aquatics', adj=c(0.5,0), srt=90, col='black', cex=8/7)
      }

      ## Pollen diversity
      plot.window(xlim=c(-150,950),ylim=c(-5,55)+60*c(-4,1.4)+60*c(-6.4*0.025, 6.4*0.025),main='',ylab='',xlab='')  ;  {
          COL='black'
          polygon(rep(TIME, each=2), c(0,rep(POLLENTAXA,each=2),0), col='forestgreen', border='forestgreen', lwd=1, lty=1, angle=-60, density=25)
          for(i in seq(0,50,5)) segments(-20,i,-9,i, lwd=0.5, col='black')
          for(i in seq(0,50,10)) text(-25,i,i, adj=c(1,0.5), col='black')
          text(-175, 25, '(d) Number of\npollen taxa per sample', adj=c(0.5,1), srt=90, col='black', cex=8/7)
          points(XX.interp, MARGALEF*5, col='grey20', type='l', lwd=1.2)
          for(i in seq(0, 9, 1.5)*5) segments(820,i,809,i, lwd=0.5, col=COL)
          for(i in seq(0, 9, 3)*5) text(825,i,i/5, adj=c(0,0.5), col=COL)
          text(975, 25, "(e) Margalef's diversity\nindex", adj=c(0.5,0), srt=90, col=COL, cex=8/7)
      }

      ## SST reconstruction
      plot.window(xlim=c(-150,950),ylim=c(-3.3,3.2)+6.5*c(-3,2.4)+6.5*c(-6.4*0.025, 6.4*0.025),main='',ylab='',xlab='')  ;  {
      COL='black'
          points(SSTs, col=makeTransparent('grey70', alpha=1), type='l', cex=0.3)
          points(XX.interp, SSTs.smooth, lwd=1.2, col='darkorchid3', type='l')
          for(i in seq(-3,3,1)) segments(-20,i,-9,i, lwd=0.5, col=COL)
          for(i in seq(-3,3,2)) text(-25,i,i, adj=c(1,0.5), col=COL)
          text(-175, -0.05, '(f) Mozambique Channel\nSSTs PC1', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
      }

      ## MAT reconstruction
      plot.window(xlim=c(-150,950),ylim=c(16.9,21)+4.1*c(-2,3.4)+4.1*c(-6.4*0.025, 6.4*0.025),main='',ylab='',xlab='')  ;  {
          COL='black'
          points(MAT, col=makeTransparent('grey70', alpha=1), type='l', cex=0.3)
          points(XX.interp, MAT.smooth, lwd=1.2, col='darkorchid3', type='l')
          for(i in seq(16.5,21.5,0.5)) segments(820,i,809,i, lwd=0.5, col=COL)
          for(i in seq(17,21.5,1)) text(825,i,i, adj=c(0,0.5), col=COL)
          text(975, 19, '(g) MD96-2048 Pollen-based\nMAT Reconstruction (°C)', adj=c(0.5,0), srt=90, col=COL, cex=8/7)
      }

      ## Come C CO2
      plot.window(xlim=c(-150, 950),ylim=c(175,300)+125*c(-1,4.4)+125*c(-6.4*0.025, 6.4*0.025),main='',ylab='',xlab='')  ;  {
          COL='black'
          points(CO2 , col=makeTransparent('grey70', alpha=1), type='l', cex=0.3)
          points(XX.interp, CO2.smooth, lwd=1.2, col='darkorchid3', type='l')
          for(i in seq(180,300,20)) segments(-20,i,-9,i, lwd=0.5, col=COL)
          for(i in seq(180,300,40)) text(-25,i,i, adj=c(1,0.5), col=COL)
          text(-175, 237.5, '(h) Dome C, Antarctica\npCO2 (ppm)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
      }

      ## LR04
      plot.window(xlim=c(-150,950),ylim=c(5,3.2)-1.8*c(0,5.4)-1.8*c(-6.4*0.025, 6.4*0.025),main='',ylab='',xlab='')  ;  {
      COL='black'
          points(LR04, col=makeTransparent('grey70', alpha=1), type='l', cex=0.3)
          points(XX.interp, LR04.smooth, lwd=1.2, col='darkorchid3', type='l')
          for(i in seq(5,3.2,-0.25)) segments(820,i,809,i, lwd=0.5, col=COL)
          for(i in seq(5,3.2,-0.5)) text(825,i,i, adj=c(0,0.5), col=COL)
          text(975, 4.1, '(i) Global ice volume (LR04)\nd18Obenthic (permil VPDB)', adj=c(0.5,0), srt=90, col=COL, cex=8/7)
      }
    dev.off()  ;  }

}


#-;
