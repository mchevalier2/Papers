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

    get_indices_from_pdf <- function(clim, pdf){
        CS=cumsum(pdf)
        CI2.5 = clim[max(which(CS <= 0.025))]
        CI97.5 = clim[min(which(CS >= 0.975))]
        opt = clim[which.max(pdf)]
        return(c(opt, CI2.5, CI97.5))
    }

    MAT=rio::import('https://github.com/mchevalier2/ClimateReconstructions/blob/master/MD96-2048_MAT_01.xlsx?raw=true', which=2)[1:181,]
    pdfpol=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/MAT_pollen_responses.csv')
    pdfpol[,-1] = pdfpol[,-1] * (pdfpol[2,1] - pdfpol[1,1])
    IDX=matrix(0, ncol=ncol(pdfpol), nrow=3)
    for(i in 2:ncol(pdfpol)){
        IDX[,i] = get_indices_from_pdf(pdfpol[,1], pdfpol[,i])
    }
    TAXA = colnames(pdfpol)[order(IDX[1,])]
    IDX = IDX[, order(IDX[1,]) ]


    RdYlBu=rev(colorRampPalette(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))(ncol(pdfpol)))

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_FigDR2.pdf"), width=5.54, height=9.5, useDingbats=FALSE)  ;  {
        par(ps=7,bg=makeTransparent("white",alpha=0),mar=rep(0,4),cex=1,cex.main=1)
        plot(0,0, type='n', ylim=c(2,ncol(pdfpol)), xlim=c(-2,29.2), axes=FALSE, frame=FALSE)
        rect(IDX[1,2],1,IDX[1,ncol(IDX)],ncol(pdfpol)+1, border='grey50', lty=2, lwd=0.8, col='grey90')
        rect(min(MAT[,2]),1,max(MAT[,2]),ncol(pdfpol)+1, border='grey50', lty=2, lwd=0.8, col='grey60')
        for(i in 2:ncol(pdfpol)){
            rect(IDX[2,i],i+0.3,IDX[3,i], i-0.3,lwd=0.1, col=RdYlBu[i])
            segments(IDX[1,i], i+0.25, IDX[1,i], i-0.25, lwd=0.8, col=ifelse(IDX[1,i] >= min(MAT[,2]) & IDX[1,i] <= max(MAT[,2]), 'black', 'white'))
            text(IDX[2,i]-1,i, TAXA[i], cex=5/7, srt=180, adj=c(0,0.5))
        }
        segments(5,1,30,1)
        for(i in seq(5,30,2.5)) segments(i,1,i,0.4)
        for(i in seq(5,30,5)) text(i,0.2, i, adj=c(1,0.5), srt=90, cex=1)
        text(17.5, -2.5, 'Mean Annual Temperature (°C)', cex=1, srt=180, adj=c(0.5,0))

        segments(5,ncol(pdfpol)+1,30,ncol(pdfpol)+1)
        for(i in seq(5,30,2.5)) segments(i,ncol(pdfpol)+1,i,ncol(pdfpol)+1.6)
        for(i in seq(5,30,5)) text(i,ncol(pdfpol)+1.8, i, adj=c(0,0.5), srt=90, cex=1)
        text(17.5, ncol(pdfpol)+4.5, 'Mean Annual Temperature (°C)', cex=1, srt=180, adj=c(0.5,1))

        polygon(c(0,-1.5,-0.7, -0.7,-1.5,0,1.5,0.7,0.7,1.5,0)-1, c(3,10,10,157,157,166,157,157,10,10,3), lwd=0.8)
        text(0.5, 11, 'Cold indicator taxa', cex=1.2, adj=c(0,1), srt=90, font=2)
        text(0.5, 156, 'Warm indicator taxa', cex=1.2, adj=c(1,1), srt=90, font=2)
        text(-1, 2+165/2, 'Taxa indicator range', cex=1.2, adj=c(0.5,0.5), srt=90, font=2)
    dev.off()  ;  }

}


#-;
