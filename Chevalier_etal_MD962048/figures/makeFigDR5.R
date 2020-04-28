## Figure 5: Boxplots of the effect of glacials and interglacials on pollen diversity
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
    ## Rarefraction analysis
    Rarefraction <- vegan::rarefy(round(POLLEN[,-1]), min(apply(round(POLLEN[,-1]),1,sum)))
    ## Species richness (S)
    S <- vegan::specnumber(POLLEN[,-1]) ## rowSums(BCI > 0) does the same... # Richness
    ## Pielou's evenness
    J <- H/log(S)
    ## Margalef’s Index
    DMG=(S-1) / log(POLLENSUM[,2])

    GLACIAL=cbind(c(0, 11.7, 29, 57, 71, 130, 191, 243, 300, 337, 374, 424, 478, 533, 563, 621, 676, 712, 761), c(11.7, 29, 57, 71, 130, 191, 243, 300, 337, 374, 424, 478, 533, 563, 621, 676, 712, 761, 800), c(F, T, T, T, F, T, F, T, F, T, F ,T, F, T, F, T, F, T, F))
    GIG=rep('Glacials', 181)
    for(i in 1:181){
      for(j in 1:nrow(GLACIAL)){
        if(GLACIAL[j,3] == 0){
          if(MAT[i,1] <= GLACIAL[j,2] & MAT[i,1] >= GLACIAL[j,1]) GIG[i]='Interglacials'
        }
      }
    }

    COL1=rgb(27,158,119,maxColorValue=255)
    COL2=rgb(117,112,179,maxColorValue=255)

    pdf(paste0(OUTPUT_FOLDER, "/Chevalier_etal_MD962048_FigDR5.pdf"), width=7.48, height=3, useDingbats=FALSE)  ;  {
        par(ps=7,bg=makeTransparent("white",alpha=0),mar=c(2,2,1,0.3),cex=1,cex.main=1,mgp=c(3,0.5,0))
        layout(matrix(1:6, ncol=6, byrow=TRUE), width=1, height=1)

        boxplot(H~GIG, main='Shannon index', col=c(COL2,COL1), notch=TRUE)
        boxplot(D1~GIG, main='Simpson index', col=c(COL2,COL1), notch=TRUE)
        boxplot(Rarefraction~GIG, main="Rarefraction analysis", col=c(COL2,COL1), notch=TRUE)
        boxplot(S~GIG, main='Number of taxa (Richness)', col=c(COL2,COL1), notch=TRUE)
        boxplot(J~GIG, main="Pielou's evenness", col=c(COL2,COL1), notch=TRUE)
        boxplot(DMG~GIG, main="Margalef's index", col=c(COL2,COL1), notch=TRUE)
    dev.off()  ;  }

}


#-;
