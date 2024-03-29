## Figure 1: The CREST method
##
## Loading necessary data

pkg2install=c()
if (! ("vegan" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'vegan')


makeAnalysis <- TRUE

if (length(pkg2install) > 0){
    s=''
    while (! s %in% c('y', 'yes', 'Y', 'YES', 'n', 'no', 'N', 'NO')){
        s <- readline(prompt=paste0("The following are required: ", paste(pkg2install, collapse=', '),". Do you want to install them? [yes/no] "))
    }
    if(s %in% c('y', 'yes', 'Y', 'YES')){
        install.packages(pkg2install)
    }else{
        print("Script aborded.")
        makeAnalysis <- FALSE
    }
}

if (makeAnalysis) {
    ##

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
    load(url('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/figures/Fig2.RData'))

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
    DomeC=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=6)[1:5784,c(1,2)]
    MALAWI=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=8)[1:295,c(1,3)]
    LEAFWAX=rio::import('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/data/IndependentRecords.xlsx', which=9)[1:177,c(1,10)]
    LEAFWAX.detrended=cbind(LEAFWAX[,1], LEAFWAX[,2] - LEAFWAX[,1]*coef(lm(LEAFWAX[,2]~ LEAFWAX[,1]))[2] - coef(lm(LEAFWAX[,2]~ LEAFWAX[,1]))[1])

    SELECTEDTAXA=rio::import('https://github.com/mchevalier2/ClimateReconstructions/raw/master/MD96-2048_MAT_01.xlsx', which=6)
    SELECTEDTAXA=SELECTEDTAXA[SELECTEDTAXA[,2] == 'Yes', 1]

    XX.interp=1:800
    MAT.smooth=gausmooth(MAT[,c(1,2)], XX.interp, mean(diff(MAT[,1])))
    CO2.smooth=gausmooth(CO2, XX.interp, mean(diff(MAT[,1])))
    SSTs.smooth=gausmooth(SSTs, XX.interp, mean(diff(MAT[,1])))
    LR04.smooth=gausmooth(LR04, XX.interp, mean(diff(MAT[,1])))
    DomeC.smooth=gausmooth(DomeC, XX.interp, mean(diff(MAT[,1])))
    MALAWI.smooth=gausmooth(MALAWI, XX.interp, mean(diff(MAT[,1])))
    LEAFWAX.smooth=gausmooth(LEAFWAX.detrended, XX.interp, mean(diff(MAT[,1])))

    cat('\n\n>>> Modern temperatures across the catchment of MD96-2048\n')
    clim=c(
        raster::values(raster::intersect(BIO1, WATERSHED)),
        raster::values(raster::intersect(BIO1, WATERSHED2[1,])),
        raster::values(raster::intersect(BIO1, WATERSHED2[2,])),
        raster::values(raster::intersect(BIO1, WATERSHED2[3,]))
    )
    cat(paste0('Mean temperature value = ', round(mean(clim, na.rm=TRUE),1), '°C\n'))
    cat(paste0('Standard deviation = ', round(sd(clim, na.rm=TRUE),1), '°C\n'))


    cat('\n\n>>> Pollen diversity and climate reconstructions\n')
    cat('\nNumber of pollen taxa used for each sample:\n')
    print(summary(apply(POLLEN[, colnames(POLLEN) %in% SELECTEDTAXA], 1, function(x) return(sum(ifelse(x>0,1,0))))))
    cat('\nPercentage of terrestrial pollen taxa used for each sample:\n')
    POLLEN_PERC = POLLEN[,-1] * 100 / apply(POLLEN[,-1], 1, sum)
    print(summary(apply(POLLEN_PERC[, colnames(POLLEN_PERC) %in% SELECTEDTAXA], 1, sum)))



    cat('\n\n>>> Correlation analysis (last 342 kyrs / last 790 kyrs)\n')
    cat(paste0('cor(MAT, CO2) = ', round(cor(MAT.smooth[1:342], CO2.smooth[1:342]),3), ' / ', round(cor(MAT.smooth, CO2.smooth),3)), '\n')
    cat(paste0('cor(MAT, SSTs) = ', round(cor(MAT.smooth[1:342], SSTs.smooth[1:342]),3), ' / ', round(cor(SSTs.smooth, CO2.smooth),3)), '\n')
    cat(paste0('cor(MAT, LR04) = ', round(cor(MAT.smooth[1:342], LR04.smooth[1:342]),3), ' / ', round(cor(MAT.smooth, LR04.smooth),3)), '\n')
    cat(paste0('cor(MAT, DomeC) = ', round(cor(MAT.smooth[1:342], DomeC.smooth[1:342]),3), ' / ', round(cor(MAT.smooth, DomeC.smooth),3)), '\n')
    cat(paste0('cor(MAT, MALAWI) = ', round(cor(MAT.smooth[1:342], MALAWI.smooth[1:342]),3), ' / ', round(cor(MAT.smooth, MALAWI.smooth),3)), '\n')
    cat(paste0('cor(MAT, LEAFWAX) = ', round(cor(MAT.smooth[1:342], LEAFWAX.smooth[1:342]),3), ' / ', round(cor(MAT.smooth, LEAFWAX.smooth),3)), '\n')

    invisible(readline(prompt="\nPress [enter] to plot the cross-correlation between the temperature indicators"))
    pairs(cbind(MAT.smooth, CO2.smooth, SSTs.smooth, LR04.smooth, DomeC.smooth, MALAWI.smooth, LEAFWAX.smooth), pch="+", col="blue")



    cat('\n\n>>> Glacial / Interglacial characteristics\n')
    GLACIAL=cbind(c(0, 11.7, 29, 57, 71, 130, 191, 243, 300, 337, 374, 424, 478, 533, 563, 621, 676, 712, 761), c(11.7, 29, 57, 71, 130, 191, 243, 300, 337, 374, 424, 478, 533, 563, 621, 676, 712, 761, 800), c(F, T, T, T, F, T, F, T, F, T, F ,T, F, T, F, T, F, T, F))

    ig_tmp=g_tmp=c()
    for(i in 1:nrow(GLACIAL)){
      if(GLACIAL[i,3] == 0){
        ig_tmp=c(ig_tmp, max(MAT[which(MAT[,1] >= GLACIAL[i,1] & MAT[,1] < GLACIAL[i,2]),2]))
      }else{
        g_tmp=c(g_tmp, min(MAT[which(MAT[,1] >= GLACIAL[i,1] & MAT[,1] < GLACIAL[i,2]),2]))
      }
      if(GLACIAL[i,1]>342) break()
    }
    cat(paste0('Mean interglacial temperature (last 342 kyrs) = ', round(mean(ig_tmp),3)), '\n')
    cat(paste0('Mean glacial temperature (last 342 kyrs) = ', round(mean(g_tmp),3)), '\n')

    ig_tmp=g_tmp=c()
    for(i in 1:nrow(GLACIAL)){
      if(GLACIAL[i,3] == 0){
        ig_tmp=c(ig_tmp, max(MAT[which(MAT[,1] >= GLACIAL[i,1] & MAT[,1] < GLACIAL[i,2]),2]))
      }else{
        g_tmp=c(g_tmp, min(MAT[which(MAT[,1] >= GLACIAL[i,1] & MAT[,1] < GLACIAL[i,2]),2]))
      }
      if(GLACIAL[i,1]>800) break()
    }
    cat(paste0('Mean interglacial temperature (last 790 kyrs) = ', round(mean(ig_tmp),3)), '\n')
    cat(paste0('Mean glacial temperature (last 790 kyrs) = ', round(mean(g_tmp),3)), '\n')



    cat('\n\n>>> Effect of temperature on pollen diversity\n')
    GIG=rep('Glacial', 181)
    for(i in 1:181){
      for(j in 1:nrow(GLACIAL)){
        if(GLACIAL[j,3] == 0){
          if(MAT[i,1] <= GLACIAL[j,2] & MAT[i,1] >= GLACIAL[j,1]) GIG[i]='Interglacial'
        }
      }
    }

  ## Diversity Shannon-Weiver
    H <- vegan::diversity(POLLEN[,-1], "shannon")
    ## Diversity Simpson
    D1 <- vegan::diversity(POLLEN[,-1], "simpson")
    ## Diversity inverse Simpson
    D2 <- vegan::diversity(POLLEN[,-1], "inv")
    ## Rarefraction analysis
    Rarefraction <- vegan::rarefy(round(POLLEN[,-1]), min(apply(round(POLLEN[,-1]),1,sum)))
    ## Species richness (S)
    S <- vegan::specnumber(POLLEN[,-1]) ## rowSums(BCI > 0) does the same... # Richness
    ## Pielou's evenness
    J <- H/log(S)
    ## Margalef’s Index
    DMG=(S-1) / log(POLLENSUM[,2])
    MAT=MAT[,2]

    cat(paste0('cor(MAT, Richness) = ', round(cor(MAT, S), 3)), '\n')
    cat(paste0('cor(MAT, Pielou s Evenness) = ', round(cor(MAT, J), 3)), '\n')
    cat(paste0('cor(MAT, Rarefraction) = ', round(cor(MAT, Rarefraction), 3)), '\n')
    cat(paste0('cor(MAT, Diversity-Shannon) = ', round(cor(MAT, H), 3)), '\n')
    cat(paste0('cor(MAT, Diversity-Simpson) = ', round(cor(MAT, D1), 3)), '\n')
    cat(paste0('cor(MAT, Diversity-Simpson inverse) = ', round(cor(MAT, D2), 3)), '\n')
    cat(paste0('cor(MAT, Margalef s index) = ', round(cor(MAT, DMG), 3)), '\n')

    invisible(readline(prompt="\nPress [enter] to plot the cross-correlation between the diversity indices"))
    pairs(cbind(MAT, H, D1, D2, Rarefraction, S, J, DMG), pch="+", col="blue")

    invisible(readline(prompt="\nPress [enter] for the next plot"))
    par(mfrow=c(1,7), mar=c(5,2,1,0))
    boxplot(S~GIG, main='Richness')
    boxplot(J~GIG, main='Evenness')
    boxplot(Rarefraction~GIG, main='Rarefraction')
    boxplot(H~GIG, main='Diversity-Shannon')
    boxplot(D1~GIG, main='Diversity-Simpson')
    boxplot(D2~GIG, main='Diversity-Simpson inverse')
    boxplot(DMG~GIG, main='Margalef s index')




}



#-;
