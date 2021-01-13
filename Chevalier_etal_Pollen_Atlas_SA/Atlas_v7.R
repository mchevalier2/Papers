library(maptools)
library(raster)
library(plotrix)
library(ncdf4)
library(gdata)
library(rgdal)
library(png)
library(gridExtra)


WORKING_FOLDER = getwd()




# ------------------------------------------------------------------------------
setwd(WORKING_FOLDER)
source('Atlas_functions.R')
source('Atlas_load_data.R')


YMAX=list()  ;  for(v in variables){  YMAX[[v]]=max(table(VARIABLES[[v]][,3]%/%CLASS_WIDTH[[v]]))  }
TEXT_SIZE=1.5
PAGE_NB=20


pdf('Atlas.pdf',width=8.27,height=11.69,useDingbats=FALSE)
  {
      par(mar=c(0,0,0,0),ps=7,cex=1.5)
      plot(0,type='n',xlim=c(0,1),ylim=c(0,1),frame=FALSE,axes=FALSE,xlab="",ylab="")
      text(0.5,0.77,"ATLAS OF",adj=c(0.5,0.5),font=2, cex=2)
      text(0.5,0.7,"SOUTHERN AFRICAN POLLEN TAXA",adj=c(0.5,0.5),font=2, cex=2)
      text(0.5,0.63,"Distributions and Climatic Affinities",adj=c(0.5,0.5),font=1, cex=2)
      text(0.5,0.58,"Manuel Chevalier",adj=c(0.5,0.5),font=1, cex=1.2)
      text(0.5,0.55,"Brian M. Chase",adj=c(0.5,0.5),font=1, cex=1.2)
      text(0.5,0.52,"Lynne J. Quick",adj=c(0.5,0.5),font=1, cex=1.2)
      text(0.5,0.490,"Louis Scott",adj=c(0.5,0.5),font=1, cex=1.2)

      rect(0.1,0.25,0.9,0.17, col='grey85')
      text(0.5,0.21,'To cite this work: Chevalier, M., Chase, B.M., Quick, L.J., Scott, L., 2021. An atlas of\nsouthern African pollen types and their climatic affinities, Palaeocology of Africa,\nXX, XXX-XXX. doi: https://www.doi.org/XXXXXXXXXXXX', adj=c(0.5,0.5), font=1, cex=1)
      #rect(0.02,0.02,0.98,0.98)

      # Empty page;back of the cover
      plot(0,type='n',xlim=c(0,1),ylim=c(0,1),frame=FALSE,axes=FALSE,xlab="",ylab="")


      { # Some Descriptive text
          layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          #PAGE_NB = addPageNumber(PAGE_NB)
          plot(0,type='n',xlim=c(0,1),ylim=c(0,1),frame=FALSE,axes=FALSE,xlab="",ylab="")
          plot(0,type='n',xlim=c(0,1),ylim=c(0,1),frame=FALSE,axes=FALSE,xlab="",ylab="")
          text(0.5,0.7,"Insert paper here.",adj=c(0.5,0.5),cex=2)
      }


      { # Large climate maps
          layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          PAGE_NB = addPageNumber(PAGE_NB)
          plot(0,type='n',xlim=c(0,1),ylim=c(0,1),frame=FALSE,axes=FALSE,xlab="",ylab="")

          layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          PAGE_NB = addPageNumber(PAGE_NB)
          plot(0,type='n',xlim=c(0,18),ylim=c(0,26.7),frame=FALSE,axes=FALSE,xlab="",ylab="")
          text(18/2,19.5,"CLIMATE MAPS",adj=c(0.5,0.5),font=2, cex=4)
          txt=paste0( 'Spatial distribution of the five climate variables selected for this atlas\n\n',
                      '           ->  Temperature of the Wettest Quarter (TWetQ; [\u00B0C])\n\n',
                      '           ->  Minimum Temperature of the Coldest Month (TColdM; [\u00B0C])\n\n',
                      '           ->  Precipitation of the Warmest Quarter (PWarmQ; [mm])\n\n',
                      '           ->  Precipitation of the Coldest Quarter (PColdQ; [mm])\n\n',
                      '           ->  Aridity Index (Aridity; [no units])\n\n',
                      '\n\nThe map of the resulting biome completes the panel of maps.'
                    )
          text(4,12.5,txt,adj=c(0,0.5),font=1, cex=2)
          rect(0,0,18,26.7,cex=2)
          rect(0.2,0.2,17.8,26.5,cex=2)

          for(v in variables) {
              layout(matrix(c( c(4,4,4,4,4),
                               c(4,2,2,2,4),
                               c(4,3,3,3,4),
                               c(1+3*(PAGE_NB%%2),4,4,4,4-3*(PAGE_NB%%2))),
                            ncol=5, byrow=TRUE), height=c(1.5, 15, 11.7, 1.5), width=c(1.5,3.5,11,3.5,1.5))
              PAGE_NB = addPageNumber(PAGE_NB)
              plotmap(v, rel=2)

              segments(extent(M1)@xmin-0.25, extent(M1)@ymin - 0.5, extent(M1)@xmin-0.25, extent(M1)@ymax)
              for(i in seq(-35,-15,5)) text(extent(M1)@xmin-0.40, i, paste0(abs(i),'S'), adj=c(1,0.4), cex=1.5)
              segments(extent(M1)@xmin-0.25, extent(M1)@ymax, extent(M1)@xmax + 0.5, extent(M1)@ymax)
              for(i in seq(15,40,5)) text(i, extent(M1)@ymax+0.2, paste0(abs(i),'E'), adj=c(0.5,0), cex=1.5)


              makeHistograms(v, rel=2, title=FALSE)
          }

          v='biomes'
          layout(matrix(c( c(4,4,4,4,4),
                           c(4,2,2,2,4),
                           c(4,3,3,3,4),
                           c(1+3*(PAGE_NB%%2),4,4,4,4-3*(PAGE_NB%%2))),
                        ncol=5, byrow=TRUE), height=c(1.5, 15, 11.7, 1.5), width=c(1.5,3.5,11,3.5,1.5))
          PAGE_NB = addPageNumber(PAGE_NB)


          par(mar=rep(0,4))
          plot(0,0,type="n",xlim=extent(M1)[1:2]+c(-2 + 1.7, -0.7),ylim=extent(M1)[3:4]+c(0,1.5),asp=1,frame=FALSE,axes=FALSE,xlab="",ylab="")
          image(RASTERS[[v]], col=unlist(COL.BIOMES), add=TRUE, interpolate=FALSE)
          plot(M1,add=TRUE, lwd=0.5)
          text(mean(extent(M1)[1:2]),extent(M1)[4]+1.8,'Southern African Biomes',adj=c(0.5,1),font=2,cex=2*8/5)

          segments(extent(M1)@xmin-0.25, extent(M1)@ymin - 0.5, extent(M1)@xmin-0.25, extent(M1)@ymax)
          for(i in seq(-35,-15,5)) text(extent(M1)@xmin-0.40, i, paste0(abs(i),'S'), adj=c(1,0.4), cex=1.5)
          segments(extent(M1)@xmin-0.25, extent(M1)@ymax, extent(M1)@xmax + 0.5, extent(M1)@ymax)
          for(i in seq(15,40,5)) text(i, extent(M1)@ymax+0.2, paste0(abs(i),'E'), adj=c(0.5,0), cex=1.5)

          val=table(VARIABLES[[v]][VARIABLES[[v]][,3] != 'Lakes' ,3])
          plot(0,0, type='n', xlim=c(-3200,4000), ylim=c(0.0,9.5), axes=FALSE, frame=FALSE)
          segments(0,0.5,0,8.5, lwd=0.5)
          for(b in 1:8) {
              text(-100, b, paste(strwrap(names(val)[b],width=40), collapse='\n'), adj=c(1,0.5), cex=1.5*8/5)
              rect(0, b-0.4, val[b], b+0.4, border='black', col=COL.BIOMES[[names(val)[b]]], lwd=0.5)
          }
          segments(0, 8.5, 3804, 8.5, lwd=0.5)
          for(i in seq(0, 3500, 750)) {
              text(i, 8.8, i, adj=c(0.5,0), cex=8/5)
              segments(i,8.5,i,8.65, lwd=0.5)
          }
          text(3804/2,9.5,'Number of quarter-degree grid cells', adj=c(0.5,0), cex=1.8*8/5)
      }



      { # Introduction to the atlas, Explanation of the elements
          PAGE_NB = addEmptyPage(PAGE_NB)

          layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          PAGE_NB = addPageNumber(PAGE_NB)
          plot(0,type='n',xlim=c(0,18),ylim=c(0,26.7),frame=FALSE,axes=FALSE,xlab="",ylab="")
          text(18/2,19.5,"DISTRIBUTIONS AND CLIMATE AFFINITIES",adj=c(0.5,0.5),font=2, cex=4)
          rect(0,0,18,26.7,cex=2)
          rect(0.2,0.2,17.8,26.5,cex=2)

          layout(matrix(c( c(3,3,3),
                           c(3,2,3),
                           c(1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2))),
                         ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          PAGE_NB = addPageNumber(PAGE_NB)
          png = readPNG('page_example_1.png')
          plot(0,0,type='n', xlim=c(0,19), ylim=c(0,27.7), axes=FALSE, frame=FALSE)
          addImg(png, x = 19-14/2, y = 27.7/2, width = 15)
          rect(5,4, 19, 23.8)

          arrows(6.5,25,6.5,23.4, code=2, length=0.1)
          text(6.5, 25.3,paste(strwrap('Name of the pollen taxon and number of composing species', width=27),collapse='\n'), adj=c(0.5,0), cex=2)

          arrows(4,21,5.4,21, code=2, length=0.1)
          text(3.8, 21,paste(strwrap('Number of species with different distribution sizes', width=27),collapse='\n'), adj=c(1,0.5), cex=2)

          arrows(4,19,5.4,19, code=2, length=0.1)
          text(3.8, 19,paste(strwrap('Statistics of the climate response (Optimum, Mean, Standard Deviation, left/right Skewness (-/+), Range).', width=27), collapse='\n'), adj=c(1,0.5), cex=2)

          arrows(4,16,5.4,16, code=2, length=0.1)
          text(3.8, 16,paste(strwrap('Distribution of the pollen taxon over the five climate variables.', width=27),collapse='\n'), adj=c(1,0.5), cex=2)

          arrows(4,11,5.4,11, code=2, length=0.1)
          text(3.8, 11,paste(strwrap('Occupation of the modern climate space (coloured bars) bythe pollen taxon (black bars).', width=27),collapse='\n'), adj=c(1,0.5), cex=2)

          arrows(4,6,5.4,6, code=2, length=0.1)
          text(3.8, 6,paste(strwrap('Modelled response of the plant species (grey) and the pollen (black). The optimum is indicated by a dashed line.', width=27),collapse='\n'), adj=c(1,0.5), cex=2)


          layout(matrix(c( c(3,3,3),
                           c(3,2,3),
                           c(1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2))),
                         ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          PAGE_NB = addPageNumber(PAGE_NB)
          png = readPNG('page_example_2.png')
          plot(0,0,type='n', xlim=c(0,19), ylim=c(0,27.7), axes=FALSE, frame=FALSE)
          addImg(png, x = 14/2, y = 27.7/2, width = 15)
          rect(0,4, 14, 23.8)

          arrows(7,25,7,23.4, code=2, length=0.1)
          segments(1,23.4,13,23.4)
          text(7, 25.3,paste(strwrap('Proportion of presence records in each biome. One grid cell can count multiple times if multiple species have been recorded.', width=70),collapse='\n'), adj=c(0.5,0), cex=2)

          arrows(15,21,9.6,21, code=2, length=0.1)
          text(15.2, 21,paste(strwrap('Distribution of the pollen taxon colour-coded by biome.', width=27),collapse='\n'), adj=c(0,0.5), cex=2)

          arrows(15,18,13.6,18, code=2, length=0.1)
          text(15.2, 18,paste(strwrap('Photographes of the pollen grains.', width=27),collapse='\n'), adj=c(0,0.5), cex=2)

          arrows(15,9.7,13.6,9.7, code=2, length=0.1)
          text(15.2, 9.7,paste(strwrap('Scatterplot of the entire climate space (black dots) with the presence records of the taxon colour-coded by biome. One plot per pair of variables.', width=27),collapse='\n'), adj=c(0,0.5), cex=2)
      }


      missing_photographes=c()
      { # The atlas itself
          for(pol in 1:length(POLTYPES)){
          #for(pol in 1:3){
              print(names(POLTYPES)[pol])

              # Defining layout page 1
              layout(matrix(c(c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21),
                              c(21,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,21),
                              c(21,19,19,19,19,19,4,4,4,4,4,5,5,5,5,5,21),
                              c(21,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,21),
                              c(21,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,21),
                              c(21,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,21),
                              c(20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21)),
                            byrow=TRUE,ncol=17),widths=c(1.5,rep(18/15,15),1.5),heights=c(1.5, c(2,6,6,8,8)/30*26.7, 1.5))

              # Top banner
              par(mar=rep(0,4))
              plot(0,0,type='n',frame=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
              if(substr(toupper(names(POLTYPES)[pol]), nchar(toupper(names(POLTYPES)[pol]))-4, nchar(toupper(names(POLTYPES)[pol]))) == '-TYPE') {
                  if(toupper(names(POLTYPES)[pol]) == 'SCROPHULARIACEAE-TYPE'){
                      text(0,0.9,toupper(names(POLTYPES)[pol]),adj=c(0,1),cex=3,font=2)
                  } else {
                      text(0,0.9,strsplit(toupper(names(POLTYPES)[pol]), '-TYPE')[[1]][1],adj=c(0,1),cex=3,font=4)
                      text(strwidth(strsplit(toupper(names(POLTYPES)[pol]), '-TYPE')[[1]][1], cex=3, font=4), 0.9, '-TYPE',adj=c(0,1),cex=3,font=2)
                  }
              } else if (substr(toupper(names(POLTYPES)[pol]), nchar(toupper(names(POLTYPES)[pol]))-4, nchar(toupper(names(POLTYPES)[pol]))) == 'ACEAE') {
                  text(0,0.9,toupper(names(POLTYPES)[pol]),adj=c(0,1),cex=3,font=2)
              } else {
                  text(0,0.9,toupper(names(POLTYPES)[pol]),adj=c(0,1),cex=3,font=4)
              }
              text(0, 0.55, paste0('(',length(unique(POLTYPES[[pol]][,3])),' species)'),adj=c(0,1),cex=2.5,font=1)

              plot(0,0,type='n',frame=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
              plot(0,0,type='n',frame=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")


              # Climate maps with distributions
              for(v in variables)  plotmap(v, pol)

              # Histograms
              for(v in variables)  makeHistograms(v, pol)

              # Plot pdfs
              for(v in variables)  plotPDFS(v, pol)

              # Plot statistical indices and species distributions
              par(mar=rep(0,4))
              plot(0,0,type='n', xlim=c(0,10.1), ylim=c(-0.5,11), axes=FALSE, frame=FALSE)

              { ## Distribution sizes
                  segments(-0.1,7.2,9.5,7.2, lwd=0.5)
                  segments(-0.1,8.8,9.5,8.8, lwd=0.5)
                  segments(2.5, 7.2,2.5,8.8, lwd=0.5)
                  segments(2.5,8.8, 3.5 ,11.2, lwd=0.5)
                  segments(3.5,11.2,11,11.2, lwd=0.5)
                  segments(-0.1,7.2,-0.1,8.8, lwd=0.5)

                  text(0.2, 8,'Number of\nspecies', adj=c(0,0.5), cex=8/5)
                  text(0.2, 10,'Size of the\ndistribution\n(# grid cells)', adj=c(0,0.5), cex=8/5)
                  LABELS=c("1-10", "10-20", "20-50", "50-100", "100-200", "200-500", ">500")

                  freqpol=as.vector(tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],length))
                  h=hist(freqpol, breaks=c(1,10,20,50,100,200,500,1000), plot=FALSE)$counts
                  for(i in 1:7){
                      text(2.5+(i-0.5), 8, h[i],  adj=c(0.5,0.5), cex=8/5)
                      segments(2.5+i, 7.2,2.5+i,8.8, lwd=0.5)
                      segments(2.5+i,8.8, 2.5+(i+1) ,11.2, lwd=0.5)
                      text(3+(i-0.8),9,LABELS[i], adj=c(0,0), srt=60, cex=8/5)
                  }
              }

              { ## Statistical indices
                  segments(2.5,6,10,6, lwd=0.5)
                  segments(-0.1,5,10,5, lwd=0.5)
                  text(3.25,5.5,'Opt.', adj=c(0.5,0.5), cex=8/5)
                  text(4.75,5.5,'Mean', adj=c(0.5,0.5), cex=8/5)
                  text(6.25,5.5,'St.Dev.', adj=c(0.5,0.5), cex=8/5)
                  text(7.75,5.5,'Skew.', adj=c(0.5,0.5), cex=8/5)
                  text(9.25,5.5,'Range', adj=c(0.5,0.5), cex=8/5)
                  segments(-0.1,5,10,5, lwd=0.5)
                  segments(2.5,6,10,6, lwd=0.5)
                  segments(-0.1,5,-0.1,0, lwd=0.5)
                  for(i in c(2.5,4,5.5,7,8.5,10)) segments(i,6,i,0, lwd=0.5)
                  for(v in 1:length(variables)){
                      text(0.2,5.5-v,VARIABLES_NAMES[[variables[v]]][4], adj=c(0,0.5), cex=8/5)
                      text(3.25,5.5-v,formatC(STATS[[names(POLTYPES)[pol]]][[variables[v]]][1], format='f', digits=ROUNDS[[variables[v]]]), adj=c(0.5,0.5), cex=8/5)
                      text(4.75,5.5-v,formatC(STATS[[names(POLTYPES)[pol]]][[variables[v]]][2], format='f', digits=ROUNDS[[variables[v]]]), adj=c(0.5,0.5), cex=8/5)
                      text(6.25,5.5-v,formatC(STATS[[names(POLTYPES)[pol]]][[variables[v]]][3], format='f', digits=ROUNDS[[variables[v]]]), adj=c(0.5,0.5), cex=8/5)
                      text(7.75,5.5-v,ifelse(STATS[[names(POLTYPES)[pol]]][[variables[v]]][4]>=0,"+","-"), adj=c(0.5,0.5), cex=8/5)
                      text(9.25,5.5-v,formatC(STATS[[names(POLTYPES)[pol]]][[variables[v]]][5], format='f', digits=ROUNDS[[variables[v]]]), adj=c(0.5,0.5), cex=8/5)
                      segments(-0.1,5-v,10,5-v, lwd=0.5)
                  }
              }

              PAGE_NB = addPageNumber(PAGE_NB)

              # Defining layout page 1
              pictures <- list.files('photographs', full.names=TRUE, pattern='*.png')
              pictures <- pictures[grep(paste0(names(POLTYPES)[pol],'-'), pictures)]
              #if(length(pictures) == 1) {
                  layout(matrix(c(c(18,18,18,18,18,18,18,18,18,18,18,18,18,18),
                                  c(18,11,11,11,12,12,12,13,13,13,14,14,14,18),
                                  c(18,16,16,16,16,16,16,16,16,19,19,19,19,18),
                                  c(18,16,16,16,16,16,16,16,16,15,15,15,15,18),
                                  c(18,16,16,16,16,16,16,16,15,15,15,15,15,18),
                                  c(18,16,16,16,16,16,16,16,15,15,15,15,15,18),
                                  c(18,16,16,16,16,16,16,16,15,15,15,15,15,18),
                                  c(18,1,1,1,2,2,2,19,15,15,15,15,15,18),
                                  c(18,1,1,1,2,2,2,19,15,15,15,15,15,18),
                                  c(18,3,3,3,4,4,4,5,5,5,6,6,6,18),
                                  c(18,7,7,7,8,8,8,9,9,9,10,10,10,18),
                                  c(18,18,18,18,18,18,18,18,18,18,18,18,18,17)),
                                ncol=14,byrow=TRUE),width=c(1.5,rep(18/12,12),1.5),heights=c(1.5, c(2,1.15,5.37/2,5.37/2,5.37/2,5.37/2,5.37/2,5.37/2,5.37,5.37)/30*27.7, 1.5))
              #}
              if(length(pictures) == 0) {
                  print(c('no photograph', names(POLTYPES)[pol]))
                  missing_photographes = c(missing_photographes, names(POLTYPES)[pol])
              }

              # Scatterplots
              scatterplot("Prec_Cold_Q","Prec_Warm_Q",pol)
              scatterplot("Prec_Cold_Q","Aridity",pol,pos=2)
              scatterplot("Prec_Cold_Q","Tmean_Wet_Q",pol)
              scatterplot("Prec_Cold_Q","Tmin_Cold_M",pol)
              scatterplot("Prec_Warm_Q","Aridity",pol,pos=2)
              scatterplot("Prec_Warm_Q","Tmean_Wet_Q",pol)
              scatterplot("Prec_Warm_Q","Tmin_Cold_M",pol)
              scatterplot("Aridity","Tmean_Wet_Q",pol)
              scatterplot("Aridity","Tmin_Cold_M",pol)
              scatterplot("Tmean_Wet_Q","Tmin_Cold_M",pol)

              # Top banner
              par(mar=c(0,0,0,0))
              plot(0,0,type='n', axes=FALSE, frame=FALSE,xlim=c(0.02,0.98), ylim=c(-0.3,1))  ;  {
                  biome <- 'Deserts and xeric shrublands'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.24,0.65, 'Deserts\nXeric', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.24,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                        }else{
                            berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=NA, border='grey40')
                            text(0.24,0.65, 'Deserts\nXeric', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                            text(0.24,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )
                       }
                    }
                  biome <- 'Mediterranean Forests, woodlands and scrubs'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.76,0.65, 'Fynbos', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.76,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=NA, border='grey40')
                          text(0.76,0.65, 'Fynbos', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.76,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )                  }
                  }
              }
              plot(0,0,type='n', axes=FALSE, frame=FALSE,xlim=c(0.02,0.98), ylim=c(-0.3,1))  ;  {
                  biome <- 'Montane grasslands and shrublands'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.24,0.65, 'Montane\nGrasslands', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.24,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=NA, border='grey40')
                          text(0.24,0.65, 'Montane\nGrasslands', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.24,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )
                     }
                  }
                  biome <- 'Flooded grasslands and savannas'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.76,0.65, 'Flooded\nGrasslands', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.76,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=NA, border='grey40')
                          text(0.76,0.65, 'Flooded\nGrasslands', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.76,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )                  }
                  }
              }
              plot(0,0,type='n', axes=FALSE, frame=FALSE,xlim=c(0.02,0.98), ylim=c(-0.3,1))  ;  {
                  biome <- 'Tropical and subtropical grasslands, savannas and shrublands'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.24,0.65, 'Savannas', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.24,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=NA, border='grey40')
                          text(0.24,0.65, 'Savannas', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.24,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )
                     }
                  }
                  biome <- 'Tropical and subtropical moist broadleaf forests'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.76,0.65, 'Moist\nforests', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.76,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=NA, border='grey40')
                          text(0.76,0.65, 'Moist\nforests', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.76,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )
                      }
                  }
              }
              plot(0,0,type='n', axes=FALSE, frame=FALSE,xlim=c(0.02,0.98), ylim=c(-0.3,1))  ;  {
                  biome <- 'Tropical and subtropical dry broadleaf forests'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.24,0.65, 'Dry\nforests', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.24,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0,0,0.48,1,rounding=0.5, col=NA, border='grey40')
                          text(0.24,0.65, 'Dry\nforests', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.24,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )
                      }
                  }
                  biome <- 'Mangroves'  ;  {
                      if(( BIOMES_UNIQUE[[pol]][biome] >= 1 ) & ( ! is.na(BIOMES_UNIQUE[[pol]][biome]) )){
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=COL.BIOMES[[biome]], border=COL.BIOMES[[biome]])
                          text(0.76,0.65, 'Mangroves', cex=TEXT_SIZE, adj=c(0.5,0.5),col='white', font=2)
                          text(0.76,0.15, paste0(BIOMES_UNIQUE[[pol]][biome],'%'), cex=TEXT_SIZE, adj=c(0.5,0),col='white' )
                      }else{
                          berryFunctions::roundedRect(0.52,0,1,1,rounding=0.5, col=NA, border='grey40')
                          text(0.76,0.65, 'Mangroves', cex=TEXT_SIZE, adj=c(0.5,0.5),col='grey40', font=2)
                          text(0.76,0.15, '0%', cex=TEXT_SIZE, adj=c(0.5,0),col='grey40' )
                      }
                  }
              }

              # Adding pollen photographes
              plot(0,0,type='n', xlim=c(0,1), ylim=c(0,1), axes=FALSE, frame=FALSE)
              for(pict in pictures){
                  png = readPNG(pict)
                  addImg(png, x = 0.5, y = 0.45, width = 1)
              }

              # Distribution and Biomes
              plot(0,0,type="n",xlim=extent(M1)[1:2]+c(.7,-.7),ylim=extent(M1)[3:4]+c(-1,0),asp=1,frame=FALSE,axes=FALSE,xlab="",ylab="")
              plot(M1,add=TRUE, col='black', border='floralwhite')
              points(POLTYPES_UNIQUE[[pol]][,1:2],pch=20,cex=0.8,col=apply(POLTYPES_UNIQUE[[pol]],1,function(x) return(COL.BIOMES[[x[8]]])))
              PAGE_NB = addPageNumber(PAGE_NB)
        }
    }


    { # Summary taxa ranking
        PAGE_NB = addEmptyPage(PAGE_NB)

        layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
        PAGE_NB = addPageNumber(PAGE_NB)
        plot(0,type='n',xlim=c(0,18),ylim=c(0,26.7),frame=FALSE,axes=FALSE,xlab="",ylab="")
        text(18/2,19.5,"POLLEN TAXA CLIMATE AFFINITIES",adj=c(0.5,0.5),font=2, cex=4)
        txt=paste0( 'The five climatic optima of each taxon as defined by the pdfs were extracted for each variable.\n\n',
                    '           ->  The pdfs of each pollen type are plotted for each variable. The pdfs are colour\n\n',
                    '                 coded based on their optimum.\n\n',
                    '           ->  The climate optima were sorted for each variable and plotted non-linearly. The\n\n',
                    '                 taxa with the lowest/highest optima are at the bottom/top and colour coded\n\n',
                    '                 based on the colour scale used in the previous climate maps.\n\n',
                    '           ->  The climate optima for the five variables were analysed using a Principal Component\n\n',
                    '                 Analysis (PCA) to highlight the taxa with similar climate niches. We represent the\n\n',
                    '                 position of the taxa on the first 3 components (>90% variance).\n\n',
                    '\n\n'
                  )
        text(2,11.5,txt,adj=c(0,0.5),font=1, cex=2)
        text(18/2,6,'Note: Some of these these plots do not include the climate tolerance or\nthe possible multimodality of the responses, and are not absolute.',adj=c(0.5,0.5),font=2, cex=2)
        rect(0,0,18,26.7,cex=2)
        rect(0.2,0.2,17.8,26.5,cex=2)
    }


    { # All the pdfs
        for(v in variables) {
            xx=seq(XRANGE[[v]][1],XRANGE[[v]][2],length.out=500)

            PDFPOLS = list()
            for(pol in names(POLTYPES)) {
                species=tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],length)
                pdf=tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],function(x) return(pdfsp(x,CLASS_WIDTH[[v]],VARIABLES_NAMES[[v]][3],xx)))
                pdfpol=rep(0,500) ; pdfpol.w=0
                for(sp in 1:length(species)){
                    if(species[sp]>=20){
                        pdfpol=pdfpol+sqrt(species[sp])*pdf[[names(species)[sp]]]
                        pdfpol.w=pdfpol.w+sqrt(species[sp])
                    }
                }
                PDFPOLS[[pol]] = pdfpol/pdfpol.w
            }

            ymx = quantile(unlist(lapply(PDFPOLS, function(x) return(max(x)))), probs=ifelse(v==variables[4], 0.9, 1))

            for(page in 1:4) {
                m1 <- matrix(c(1, 1, 1, 14, 15, 49, 16:48), ncol=3, byrow=TRUE)
                m1 <- cbind(c(1, 2:13), m1)
                m1 <- rbind(rep(54, 4), m1, c(54, 50:52), rep(54, 4))
                if(PAGE_NB %% 2 == 0) {
                    m1 <- cbind(c(rep(54, 15), 53), m1, rep(54, 16))
                } else {
                    m1 <- cbind(rep(54, 16), m1, c(rep(54, 15), 53))
                }

                layout(m1, height=c(1.5, 0.5, rep(25.7/12, 12), 0.5, 1.5), width=c(1.5,c(0.7, 5, 5, 5)/16*18,1.5))
                par(mar=c(0,0,0,0), ps=7)

                plot(NA, NA, type='n', xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', frame=FALSE, axes=FALSE)
                text(0.5, 1, paste0(VARIABLES_NAMES[[v]][1], ' ', VARIABLES_NAMES[[v]][5], ' [page ', page, ']'), font=2, cex=3, adj=c(0.5, 1))

                for(i in 1:12) {
                    plot(0,0,type='n',xlim=c(0,1),ylim=c(0,ymx),
                             frame=FALSE, axes=FALSE, xlab="",ylab="", xaxs='i', yaxs='i')
                    segments(1, 0, 1, ymx, lwd=0.5)
                    for(j in axTicks(2)) {
                        if (j <= max(ymx) & j > 0) {
                            segments(1, j, 0.9, j, lwd=0.5)
                            text( 0.85, j, j, cex=1.3, adj=c(1, 0.45))
                        }
                    }
                    hbars=axTicks(2)
                }

                for(pol in names(POLTYPES)[(page-1)*35+1:35]) {
                    plot(0,0,type='n',xlim=range(xx)+c(0,0.05)*diff(range(xx)),ylim=c(0,ymx),
                             frame=FALSE, axes=FALSE, xlab="",ylab="", xaxs='i', yaxs='i')

                    for(h in hbars) segments(xx[1], h, xx[500], h, lwd=0.5, col='grey80')
                    polygon(xx[c(1,1:500, 500)],c(0, PDFPOLS[[pol]], 0), lwd=1.5, col=VARIABLE.COL1000[[v]][which.max(PDFPOLS[[pol]])])
                    segments(xx[which.max(PDFPOLS[[pol]])],0,xx[which.max(PDFPOLS[[pol]])],max(PDFPOLS[[pol]]),lwd=1.2,lty=4)

                    rect(xx[495], ymx*0.9, xx[495]-strwidth(pol, cex=1.5), ymx*0.9-strheight(pol, cex=1.5)*1.5, col="white", border=NA)

                    if(substr(pol, nchar(pol)-4, nchar(pol)) == '-type') {
                        if(toupper(pol) == 'SCROPHULARIACEAE-TYPE'){
                            text(xx[495], ymx*0.9, pol, adj=c(1,1), cex=1.5)
                        } else {
                            text(xx[495], ymx*0.9, '-type', adj=c(1,1), cex=1.5)
                            text(xx[495]-strwidth('-type', cex=1.5), ymx*0.9, strsplit(pol, '-type')[[1]][1], adj=c(1,1), cex=1.5, font=3)
                        }
                    } else if (substr(pol, nchar(pol)-4, nchar(pol)) == 'aceae') {
                        text(xx[495], ymx*0.9, pol, adj=c(1,1), cex=1.5)
                    } else {
                        text(xx[495], ymx*0.9, pol, adj=c(1,1), cex=1.5, font=3)
                    }
                }


                plot(NA, NA, type='n', xlim=c(-30,530), ylim=c(-0.1,1.1), xaxs='i', yaxs='i', frame=FALSE, axes=FALSE)
                for (i in 1:500) {
                    rect(i-1, 0.30, i, 0.7, border=VARIABLE.COL1000[[v]][i], col=VARIABLE.COL1000[[v]][i], lwd=0.2)
                }
                rect(0, 0.3, 500, 0.7, lwd=0.5)
                for(i in seq(XRANGE[[v]][1], XRANGE[[v]][2], CLASS_WIDTH[[v]]*2)) {
                    j = min(which(xx >= i))
                    segments(j, 0.7, j, 0.75, lwd=0.5)
                    text(j, 0.8, i, adj=c(0.5, 0), cex=1.3)
                }
                for(i in seq(XRANGE[[v]][1]+CLASS_WIDTH[[v]], XRANGE[[v]][2], CLASS_WIDTH[[v]]*2)) {
                    j = min(which(xx >= i))
                    segments(j, 0.3, j, 0.25, lwd=0.5)
                    text(j, 0.2, i, adj=c(0.5, 1), cex=1.3)
                }

                for(i in 1:3) {
                    plot(0,0,type='n',xlim=range(xx)+c(0,0.05)*diff(range(xx)),ylim=c(0,1),
                             frame=FALSE, axes=FALSE, xlab="",ylab="", xaxs='i', yaxs='i')
                    segments(xx[1], 1, xx[500], 1, lwd=0.5)
                    for(j in axTicks(1)) {
                        if (j <= max(xx)) {
                            segments(j, 1, j, 0.85, lwd=0.5)
                            text(j, 0.75, j, cex=1.3, adj=c(0.5, 1))
                        }
                    }
                }



                PAGE_NB = addPageNumber(PAGE_NB)


            }
        }
    }


    { # Taxa relative ranking
        clim_opt=list()
        for(v in variables) clim_opt[[v]] = unlist(sapply(names(POLTYPES), function(x) return(STATS[[x]][[v]][1])))
        clim_opt <- as.data.frame(matrix(unlist(clim_opt), ncol=5, byrow=FALSE))
        colnames(clim_opt) <- variables
        rownames(clim_opt) <- names(POLTYPES)

        layout(matrix(c( c(12,12,12,12),
                         c(12, 1, 2,12),
                         c(12, 3, 4,12),
                         c(12, 5, 6,12),
                         c(12, 7, 8,12),
                         c(12, 9,10,12),
                         c(11,12,12,12)),
                        ncol=4, byrow=TRUE), height=c(1.5, rep(26.7/5, 5), 1.5), width=c(1.5,c(1,19)/20*18,1.5))

        par(mar=c(0,0,0,0))

        for(v in variables) {
            w <- order(clim_opt[, v])[1:(nrow(clim_opt)/2+1)]
            val=table(VARIABLES[[v]][,3]%/%CLASS_WIDTH[[v]])
            plot(0, 0, xlim=c(0,1), ylim=c(0,1), type='n', axes=FALSE, frame=FALSE, xlab='', ylab='', main='')
            text(0.5, 0.5, VARIABLES_NAMES[[v]][6], cex=2, font=2, srt=-90)

            plot(0, 0, xlim=c(2.5, nrow(clim_opt)/2+1), ylim=c(0,1), type='n', axes=FALSE, frame=FALSE, xlab='', ylab='', main='')
            idx=1
            for(pol in rownames(clim_opt)[w]){
                text(idx, 0.55, pol, adj=c(1, 0.5), srt=-90, cex=1.4)
                text(idx, 0.15, formatC(clim_opt[pol, v], format='f', digits=ROUNDS[[v]]), adj=c(0,0.5), srt=-90, cex=1.4)
                rect(idx-0.5, 0.2, idx+0.5, 0.5, border=VARIABLE.COL[[v]][which(names(val) == clim_opt[pol, v]%/%CLASS_WIDTH[[v]])], col=VARIABLE.COL[[v]][which(names(val) == clim_opt[pol, v]%/%CLASS_WIDTH[[v]])])
                idx = idx + 1
            }
            rect(0.5, 0.2, idx-0.5, 0.5, lwd=0.5, col=NA, border='black')
        }
        PAGE_NB = addPageNumber(PAGE_NB)

        layout(matrix(c( c(12,12,12,12),
                         c(12, 1, 2,12),
                         c(12, 3, 4,12),
                         c(12, 5, 6,12),
                         c(12, 7, 8,12),
                         c(12, 9,10,12),
                         c(12,12,12,11)),
                        ncol=4, byrow=TRUE), height=c(1.5, rep(26.7/5, 5), 1.5), width=c(1.5,c(19,1)/20*18,1.5))

        par(mar=c(0,0,0,0))
        for(v in variables) {
            w <- order(clim_opt[, v])[(nrow(clim_opt)/2+1.5):nrow(clim_opt)]
            val=table(VARIABLES[[v]][,3]%/%CLASS_WIDTH[[v]])

            plot(0, 0, xlim=c(1, nrow(clim_opt)/2-2), ylim=c(0,1), type='n', axes=FALSE, frame=FALSE, xlab='', ylab='', main='')
            idx=1
            for(pol in rownames(clim_opt)[w]){
              text(idx, 0.55, pol, adj=c(1, 0.5), srt=-90, cex=1.5)
              text(idx, 0.15, formatC(clim_opt[pol, v], format='f', digits=ROUNDS[[v]]), adj=c(0,0.5), srt=-90, cex=1.5)
              rect(idx-0.5, 0.2, idx+0.5, 0.5, border=VARIABLE.COL[[v]][which(names(val) == clim_opt[pol, v]%/%CLASS_WIDTH[[v]])], col=VARIABLE.COL[[v]][which(names(val) == clim_opt[pol, v]%/%CLASS_WIDTH[[v]])])
              idx = idx + 1
            }
            rect(0.5, 0.2, idx-0.5, 0.5, lwd=0.5, col=NA, border='black')

            plot(0, 0, xlim=c(0,1), ylim=c(0,1), type='n', axes=FALSE, frame=FALSE, xlab='', ylab='', main='')
            text(0.5, 0.5, VARIABLES_NAMES[[v]][6], cex=2, font=2, srt=-90)
        }
        PAGE_NB = addPageNumber(PAGE_NB)
    } # End Taxa relative ranking


    { # Plotting PCAs of climate optima
        colnames(clim_opt) = c('TWetQ', 'TColdM', 'PWarmQ', 'PColdQ', 'Aridity')
        acp=ade4::dudi.pca(clim_opt, center=TRUE, scale=TRUE, scannf=FALSE, nf=5)

        myplot12 = factoextra::fviz_pca_biplot(acp, axes=c(2,1),
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,     # Avoid text overlapping
                                               labelsize = 3,
                                               title = 'PC2 (x-axis) v. PC1 (y-axis)')
        myplot12 = myplot12 + ggplot2::theme(text = ggplot2::element_text(size = 5),
                                             plot.margin =ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
                                             plot.title = ggplot2::element_text(size=10, face=2),
                                             axis.title = ggplot2::element_text(size = 7),
                                             axis.text = ggplot2::element_text(size = 7))

        myplot13 = factoextra::fviz_pca_biplot(acp, axes=c(3,1),
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,     # Avoid text overlapping
                                               labelsize = 2,
                                               title = 'PC3 (x-axis) v. PC1 (y-axis)')
        myplot13 = myplot13 + ggplot2::theme(text = ggplot2::element_text(size = 5),
                                             plot.margin =ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
                                             plot.title = ggplot2::element_text(size=10, face=2),
                                             axis.title = ggplot2::element_text(size = 7),
                                             axis.text = ggplot2::element_text(size = 7))


        myplot23 = factoextra::fviz_pca_biplot(acp, axes=c(3, 2),
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,     # Avoid text overlapping
                                               labelsize = 2,
                                               title = 'PC3 (x-axis) v. PC2 (y-axis)')
        myplot23 = myplot23 + ggplot2::theme(text = ggplot2::element_text(size = 5),
                                             plot.margin =ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
                                             plot.title = ggplot2::element_text(size=10, face=2),
                                             axis.title = ggplot2::element_text(size = 7),
                                             axis.text = ggplot2::element_text(size = 7))

        page1 = ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0.92,0.95,label=PAGE_NB), size=2.5) +
                ggplot2::xlab(NULL) +
                ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
                ggplot2::geom_text(hjust=1, vjust=1)

        page2 = ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0.1,0.95,label=PAGE_NB), size=2.5) +
                ggplot2::xlab(NULL) +
                ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
                ggplot2::geom_text(hjust=0, vjust=1)


        layout_matrix1 = matrix(c( c(3,3,3),
                                   c(3,1,3),
                                   c(2,3,3)),
                                ncol=3, byrow=TRUE)
        grid.arrange(myplot12, page1, layout_matrix = layout_matrix1, widths=c(1.5,18,1.5), heights=c(1.5,26.7,1.5))

        layout_matrix2 = matrix(c( c(4,4,4),
                                   c(4,1,4),
                                   c(4,2,4),
                                   c(4,4,3)),
                                ncol=3, byrow=TRUE)
        grid.arrange(myplot13, myplot23, page2, layout_matrix = layout_matrix2, widths=c(1.5,18,1.5), heights=c(1.5,rep(26.7/2, 2),1.5))

        PAGE_NB = PAGE_NB + 2
    } # End plot PCAs


    { # Taxonomy
        PAGE_NB = addEmptyPage(PAGE_NB)

        layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
        PAGE_NB = addPageNumber(PAGE_NB)
        plot(0,type='n',xlim=c(0,18),ylim=c(0,26.7),frame=FALSE,axes=FALSE,xlab="",ylab="")
        text(18/2,19.5,"POLLEN TAXONOMY",adj=c(0.5,0.5),font=2, cex=4)
        rect(0,0,18,26.7,cex=2)
        rect(0.2,0.2,17.8,26.5,cex=2)
    }


    { # Printing list of taxa
        max_lines = 60
        START = TRUE
        for(pol in names(POLTYPES)){
            if(START) {
                line_count = new_page(max_lines, max_lines)
                PAGE_NB = line_count[2]  ;  line_count = line_count[1]
                START = FALSE
            }
            if(line_count > 1) {
                line_count = new_page(line_count, max_lines) ## Adding an empty line before new taxon names, if not at the top of the page.
                PAGE_NB = line_count[2]  ;  line_count = line_count[1]
            }
            if(line_count >= max_lines - 1) {
                line_count = new_page(max_lines, max_lines) ## Going to the next page if the last item is a new taxon name
                PAGE_NB = line_count[2]  ;  line_count = line_count[1]
            }

            if(substr(pol, nchar(pol)-4, nchar(pol)) == '-type') {
                if(toupper(pol) == 'SCROPHULARIACEAE-TYPE'){
                    text(0, line_count, pol, font=2, cex=2, adj=c(0,0.5))
                } else {
                    text(0, line_count, strsplit(pol, '-type')[[1]][1], font=4, cex=2, adj=c(0,0.5))
                    text(strwidth(strsplit(pol, '-type')[[1]][1], font=4, cex=2), line_count, '-type', font=2, cex=2, adj=c(0,0.5))
                }
            } else if (substr(pol, nchar(pol)-4, nchar(pol)) == 'aceae') {
                text(0, line_count, pol, font=2, cex=2, adj=c(0,0.5))
            } else {
                text(0, line_count, pol, font=4, cex=2, adj=c(0,0.5))
            }

            #text(0, line_count, pol, font=2, cex=2, adj=c(0,0.5))
            line_count = new_page(line_count, max_lines)
            PAGE_NB = line_count[2]  ;  line_count = line_count[1]
            genus =  unique(apply(TAXONOMY[[pol]], 1, function(x) return(paste(x[1:6], collapse=' > '))))
            for(g in genus) {
                text(0.05, line_count, g, cex=1.5, adj=c(0,0.5))
                line_count = new_page(line_count, max_lines, step=0.7)
                PAGE_NB = line_count[2]  ;  line_count = line_count[1]
                gg = strsplit(g, '> ')[[1]][6]
                sp = TAXONOMY[[pol]][TAXONOMY[[pol]][,6] == gg, 7]
                sp = unlist(lapply(sapply(sp, function(x) strsplit(x, paste0(gg, ' '))), function(y) return(paste(substr(gg,1,1), y[2], sep='. '))))
                sp = strwrap(paste(sp, collapse=', '),width=170)
                for(j in sp) {
                    text(0.1, line_count, j, cex=1.2, adj=c(0,0.5), font=3)
                    line_count = new_page(line_count, max_lines, step=0.7)
                    PAGE_NB = line_count[2]  ;  line_count = line_count[1]
                }
                line_count = new_page(line_count, max_lines, step=0.3)
                PAGE_NB = line_count[2]  ;  line_count = line_count[1]
            }
        }
    } # End list of taxa
}  ;  dev.off()







#-;
