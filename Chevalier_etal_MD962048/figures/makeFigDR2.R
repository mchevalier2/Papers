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


    #makeTransparent <- function(..., alpha=0.5) {
    #    if(alpha>1) alpha=1
    #    if(alpha<0) alpha=0
    #    alpha = floor(255*alpha)
    #    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    #    .makeTransparent = function(col, alpha) {
    #      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    #    }
    #    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    #    return(newColor)
    #}

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_FigDR2.pdf"), width=7.48, height=4, useDingbats=FALSE)  ;  {
        par(ps=7,bg=makeTransparent("white",alpha=0),mar=rep(0,4),cex=1,cex.main=1)

    dev.off()  ;  }

}


#-;
