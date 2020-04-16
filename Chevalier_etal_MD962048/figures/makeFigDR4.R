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

    cat(">>> Loading data.\n")
    MAT=rio::import('https://github.com/mchevalier2/ClimateReconstructions/raw/master/MD96-2048_MAT_01.xlsx', which=2)[1:181,]
    POLLEN=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,-c(2,3)]
    POLLENSUM=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=4)[1:181,c(1,2,3)]

    ## Diversity Shannon-Weiver
    H <- vegan::diversity(POLLEN[,-1], "shannon")
    ## Diversity Simpson
    D1 <- vegan::diversity(POLLEN[,-1], "simpson")
    ## Diversity Fisher's alpha
    alpha <- vegan::fisher.alpha(round(POLLEN[,-1]))
    ## Species richness (S)
    S <- vegan::specnumber(POLLEN[,-1]) ## rowSums(BCI > 0) does the same... # Richness
    ## Pielou's evenness
    J <- H/log(S)
    ## Margalef’s Index
    DMG=(S-1) / log(POLLENSUM[,2])


    COL='black'
    COL2=rgb(27,158,119,maxColorValue=255)

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_FigDR4.pdf"), width=7.48, height=6.4, useDingbats=FALSE)  ;  {
        par(ps=7,bg=makeTransparent("white",alpha=0),mar=rep(0,4),cex=1,cex.main=1)
        layout(matrix(1:6, ncol=2, byrow=TRUE), width=1, height=1)

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent('black', alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(20, max(MAT[,2]), 'A', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(H)+diff(range(H))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], H, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(0.75,3,0.25)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(0.75,3,0.5)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(H)+diff(range(H))/2, 'Shannon-Weiver Index (H)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent('black', alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(20, max(MAT[,2]), 'B', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(D1)+diff(range(D1))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], D1, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(-0.25,0.9,0.05)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(-0.3,0.9,0.1)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(D1)+diff(range(D1))/2, 'Simpson Index (D1)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent('black', alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(20, max(MAT[,2]), 'C', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(alpha)+diff(range(alpha))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], alpha, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(4,20,2)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(4,20,4)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(alpha)+diff(range(alpha))/2, "Fischer's alpha", adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent('black', alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(20, max(MAT[,2]), 'D', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(S)+diff(range(S))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], S, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(10,50,5)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(10,50,10)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(S)+diff(range(S))/2, 'Number of taxa (S)', adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent('black', alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(20, max(MAT[,2]), 'E', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(J)+diff(range(J))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], J, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(0.25,0.85,0.05)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(0.3,0.8,0.1)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(J)+diff(range(J))/2, "Pielou's evennes (J)", adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }

        plot.new()  ;  { ## SSTs
            plot.window(xlim=c(-100,900),ylim=range(MAT[,2])+diff(range(MAT[,2]))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT, col=makeTransparent('black', alpha=1), type='l', cex=0.3)
                for(i in seq(16.5,21.5,0.5)) segments(-20,i,-9,i, lwd=0.5, col=COL)
                for(i in seq(17,21.5,1)) text(-25,i,i, adj=c(1,0.5), col=COL)
                text(-115, min(MAT[,2])+diff(range(MAT[,2]))/2, 'MD96-2048 Pollen-based MAT Reconstruction (°C)', adj=c(0.5,1), srt=90, col=COL, cex=8/7)
                rect(-9,min(MAT[,2])-0.02*diff(range(MAT[,2])),809,max(MAT[,2])+0.02*diff(range(MAT[,2])),lwd=0.5)
                for(i in seq(0,800,25)){  segments(i,min(MAT[,2])-0.02*diff(range(MAT[,2])),i,min(MAT[,2])-ifelse(i%%50 == 0, 0.03,0.025)*diff(range(MAT[,2])), lwd=0.5)  }
                for(i in seq(0,800,100)){  text(i,min(MAT[,2])-0.04*diff(range(MAT[,2])), i, cex=1, adj=c(0.5, 1))  }
                text(400, min(MAT[,2])-0.1*diff(range(MAT[,2])), 'Age (calendar yr BP x1000)', adj=c(0.5,0.5), cex=8/7)
                text(20, max(MAT[,2]), 'F', cex=2.5, font=2, adj=c(0,1))
            }
            plot.window(xlim=c(-100,900),ylim=range(DMG)+diff(range(DMG))*c(-0.1,0.02),main='',ylab='',xlab='')  ;  {
                points(MAT[,1], DMG, col=makeTransparent(COL2, alpha=1), type='l', cex=0.3)
                for(i in seq(3,9,0.5)) segments(820,i,809,i, lwd=0.5, col=COL2)
                for(i in seq(3,9,1)) text(825,i,i, adj=c(0,0.5), col=COL2)
                text(920, min(DMG)+diff(range(DMG))/2, "Margelef's index (DMG)", adj=c(0.5,0), srt=90, col=COL2, cex=8/7)
            }
        }
    dev.off()  ;  }

}


#-;
