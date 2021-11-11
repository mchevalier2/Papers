
{
    ## Fit the species pdfs.
    pdfsp=function(clim,width,type,xx,npix=20){
        if(length(clim)<npix){return(NA)}
        xbar=mean(clim,na.rm=TRUE)
        s2=var(clim,na.rm=TRUE)
        if(type=="normal"){
            return(1/sqrt(2*pi*s2)*exp(-((xx-xbar)**2)/(2*s2)))
        }else{
            if(xx[1]<=0) xx[1]=1e-12
            mu=log(xbar)-log(1+s2/xbar**2)/2
            sigma2=log(1+s2/xbar**2)
            return(1/sqrt(2*pi*sigma2*xx**2)*exp(-((log(xx)-mu)**2)/(2*sigma2)))
        }
    }


    ## Make the scatterplot of Variable 1 v. variable 2
    scatterplot=function(v1,v2,pol,pos=1){
        par(mar=c(2,1.8,0.3,0.2),mgp=c(3,0.4,0))
        x1=XRANGE[[v1]] ; x1[1]=x1[1]-diff(x1)/15
        x2=XRANGE[[v2]] ; x2[1]=x2[1]-diff(x2)/15
        plot(0,0,type='n',xlim=x1,ylim=x2,frame=FALSE,axes=FALSE,xlab="",ylab="")
        rect(x1[1], x2[1], x1[2], x2[2], lwd=0.5)

        points(VARIABLES[[v1]][,3],VARIABLES[[v2]][,3],col="black", bg=NA,pch=21,cex=0.6, lwd=0.3)
        points(POLTYPES_UNIQUE[[pol]][,v1],POLTYPES_UNIQUE[[pol]][,v2],col=unlist(apply(POLTYPES_UNIQUE[[pol]],1,function(x) return(COL.BIOMES[[x[8]]]))),pch=16,cex=0.8)

        #w=which(!is.na(VARIABLES[["Aridity"]][,3]))
        #if(pos==1){
        #    text(XRANGE[[v1]][2]-diff(XRANGE[[v1]])/7,XRANGE[[v2]][2]-diff(XRANGE[[v2]])/10,eval(substitute(expression(rho == b*a ), list(a = "%",b = round(cor(VARIABLES[[v1]][w,3],VARIABLES[[v2]][w,3]),3)*100))),col="grey70",adj=c(1,0),cex=8/5)
        #    text(XRANGE[[v1]][2]-diff(XRANGE[[v1]])/7,XRANGE[[v2]][2]-diff(XRANGE[[v2]])/5,eval(substitute(expression(rho == b*a ), list(a = "%",b = round(cor(POLTYPES_UNIQUE[[pol]][,v1],POLTYPES_UNIQUE[[pol]][,v2]),3)*100))),col="black",adj=c(1,0),cex=8/5)
        #}else{
        #    text(XRANGE[[v1]][2]-diff(XRANGE[[v1]])/7,XRANGE[[v2]][1]+diff(XRANGE[[v2]])/7,eval(substitute(expression(rho == b*a ), list(a = "%",b = round(cor(VARIABLES[[v1]][w,3],VARIABLES[[v2]][w,3]),3)*100))),col="grey70",adj=c(1,1),cex=8/5)
        #    text(XRANGE[[v1]][2]-diff(XRANGE[[v1]])/7,XRANGE[[v2]][1]+diff(XRANGE[[v2]])/20,eval(substitute(expression(rho == b*a ), list(a = "%",b = round(cor(POLTYPES_UNIQUE[[pol]][,v1],POLTYPES_UNIQUE[[pol]][,v2]),3)*100))),col="black",adj=c(1,1),cex=8/5)
        #}

        segments(XRANGE[[v1]][1]-diff(XRANGE[[v1]])/30,quantile(VARIABLES[[v2]][,3],0.1,na.rm=TRUE),XRANGE[[v1]][1]-diff(XRANGE[[v1]])/30,quantile(VARIABLES[[v2]][,3],0.9,na.rm=TRUE),col="grey70",lwd=0.7)
        segments(XRANGE[[v1]][1]-diff(XRANGE[[v1]])/20,quantile(POLTYPES_UNIQUE[[pol]][,v2],0.1,na.rm=TRUE),XRANGE[[v1]][1]-diff(XRANGE[[v1]])/20,quantile(POLTYPES_UNIQUE[[pol]][,v2],0.9,na.rm=TRUE),col="black",lwd=0.7)

        segments(quantile(VARIABLES[[v1]][,3],0.1,na.rm=TRUE),XRANGE[[v2]][1]-diff(XRANGE[[v2]])/30,quantile(VARIABLES[[v1]][,3],0.9,na.rm=TRUE),XRANGE[[v2]][1]-diff(XRANGE[[v2]])/30,col="grey70",lwd=0.7)
        segments(quantile(POLTYPES_UNIQUE[[pol]][,v1],0.1,na.rm=TRUE),XRANGE[[v2]][1]-diff(XRANGE[[v2]])/20,quantile(POLTYPES_UNIQUE[[pol]][,v1],0.9,na.rm=TRUE),XRANGE[[v2]][1]-diff(XRANGE[[v2]])/20,col="black",lwd=0.7)

        axis(2,lwd.ticks=0.5,lwd=0,pos=x1[1],tck=-0.03)
        axis(1,lwd.ticks=0.5,lwd=0,pos=x2[1],tck=-0.03)
        mtext(paste(VARIABLES_NAMES[[v1]][4],VARIABLES_NAMES[[v1]][5],sep=" "),side=1,line=0.8,adj=0.5667)
        mtext(paste(VARIABLES_NAMES[[v2]][4],VARIABLES_NAMES[[v2]][5],sep=" "),side=2,line=0.9,adj=0.5667)
    }


    ## Plot the map of variable v, with pollen type 'pol' in option.
    plotmap=function(v, pol=NA, rel=1){
        par(mar=rep(0,4))
        plot(0,0,type="n",xlim=extent(M1)[1:2]+c(-rel + 1.7, -0.7),ylim=extent(M1)[3:4]+c(0,1.5),asp=1,frame=FALSE,axes=FALSE,xlab="",ylab="")
        image(RASTERS[[v]], col=VARIABLE.COL[[v]], add=TRUE, interpolate=FALSE)
        val=table(VARIABLES[[v]][,3]%/%CLASS_WIDTH[[v]])
        contour(RASTERS[[v]],levels=seq(as.numeric(names(val))[1],as.numeric(names(val))[length(val)],2),drawlabels=FALSE,add=TRUE,col="grey50",lwd=0.5)
        plot(M1,add=TRUE, lwd=0.5)
        text(mean(extent(M1)[1:2]),extent(M1)[4]+1.8,VARIABLES_NAMES[[v]][1],adj=c(0.5,1),font=2,cex=rel*8/5)
        if(!is.na(pol))   points(POLTYPES_UNIQUE[[pol]][,1:2],pch=20,cex=0.3)
    }


    ## Plot the histograms of variable v, with pollen type 'pol' in option.
    makeHistograms=function(v, pol=NA, rel=1, title=TRUE){
        par(mar=c(2,1.8+15*(rel-1),1,0.2)/rel**2,mgp=c(3,0.4*rel,0))
        val=table(VARIABLES[[v]][,3]%/%CLASS_WIDTH[[v]])
        plot(0,0,type='n',xlim=XRANGE[[v]],ylim=c(0,YMAX[[v]]),frame=FALSE,axes=FALSE,xlab="",ylab="")
        if(title) title(paste(VARIABLES_NAMES[[v]][4],sep=""),font=2,cex.main=TEXT_SIZE,line=-0.25)
        mtext("Number of quarter-degree grid cells",side=2,line=0.9*rel, cex=rel)
        mtext(paste(VARIABLES_NAMES[[v]][5],sep=""),side=1,line=0.8*rel, cex=rel)
        for(i in 1:length(val)){
            rect((i-1+as.numeric(names(val)[1]))*CLASS_WIDTH[[v]],0,(i+as.numeric(names(val)[1]))*CLASS_WIDTH[[v]],val[i],col=VARIABLE.COL[[v]][i])
        }
        if(!is.na(pol)) {
            val=table(POLTYPES_UNIQUE[[pol]][,v]%/%CLASS_WIDTH[[v]])
            for(i in 1:length(val)){
                rect((i-1+as.numeric(names(val)[1]))*CLASS_WIDTH[[v]],0,(i+as.numeric(names(val)[1]))*CLASS_WIDTH[[v]],val[i],col="black")
            }
        }
        axis(2,at=c(0,max(YMAX[[v]],axTicks(2))),labels=c("",""),tck=0,lwd=0.5,pos=XRANGE[[v]][1])
        axis(2,lwd.ticks=0.5,lwd=0,pos=XRANGE[[v]][1],tck=-0.03/rel,cex.axis=TEXT_SIZE)
        axis(1,at=XRANGE[[v]],labels=c("",""),tck=0,lwd=0.5,pos=0)
        axis(1,lwd.ticks=0.5,lwd=0,pos=0,tck=-0.03/rel,cex.axis=TEXT_SIZE)
    }


    ## Fit the pdf of the pollen pol for variable v.
    ## The first 5 elements are statistics of the pdf, not the pdf itself.
    fitPDFS=function(v, pol){
        res=rep(NA,5)
        xx=seq(XRANGE[[v]][1],XRANGE[[v]][2],length.out=500)

        species=tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],length)
        pdf=tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],function(x) return(pdfsp(x,CLASS_WIDTH[[v]],VARIABLES_NAMES[[v]][3],xx)))

        pdfpol=rep(0,500) ; pdfpol.w=0
        for(sp in 1:length(species)){
            if(species[sp]>=20){
                pdfpol=pdfpol+sqrt(species[sp])*pdf[[names(species)[sp]]]
                pdfpol.w=pdfpol.w+sqrt(species[sp])
            }
        }
        pdfpol=pdfpol/pdfpol.w


        res[1] = xx[which.max(pdfpol)]
        res[2] = stats::weighted.mean(xx, pdfpol)
        res[3] = radiant.data::weighted.sd(xx, pdfpol, na.rm=TRUE)
        res[4] = res[2]-res[1]
        res[5] = diff(quantile(POLTYPES[[pol]][,v], c(0, 1), na.rm=TRUE))
        return(c(res, pdfpol))
    }


    ## Plot the pdf of the pollen pol for variable v.
    plotPDFS=function(v, pol){
        xx=seq(XRANGE[[v]][1],XRANGE[[v]][2],length.out=500)

        species=tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],length)
        pdf=tapply(POLTYPES[[pol]][,v],POLTYPES[[pol]][,"species"],function(x) return(pdfsp(x,CLASS_WIDTH[[v]],VARIABLES_NAMES[[v]][3],xx)))

        par(mar=c(2,1.8,1,0.2),mgp=c(3,0.4,0))
        plot(0,0,type='n',xlim=range(xx),ylim=c(0,max(unlist(lapply(pdf,function(x) return(max(c(0,x),na.rm=TRUE)))),na.rm=TRUE)),frame=FALSE,axes=FALSE,xlab="",ylab="")
        title(main=paste(VARIABLES_NAMES[[v]][4],sep=""),cex.main=TEXT_SIZE,family="sans",line=-0.25)

        pdfpol=rep(0,500) ; pdfpol.w=0
        for(sp in 1:length(species)){
            if(species[sp]>=20){
                points(xx,pdf[[names(species)[sp]]],type='l',lwd=1,col="grey60")
                pdfpol=pdfpol+sqrt(species[sp])*pdf[[names(species)[sp]]]
                pdfpol.w=pdfpol.w+sqrt(species[sp])
            }
        }
        pdfpol=pdfpol/pdfpol.w

        for(i in seq(1,498,2)) rect(xx[i],0,xx[i+2],pdfpol[i+1],col=VARIABLE.COL1000[[v]][i],border=VARIABLE.COL1000[[v]][i])
        segments(xx[which.max(pdfpol)],0,xx[which.max(pdfpol)],max(pdfpol),lwd=1.2,lty=4)
        points(xx,pdfpol,type='l',lwd=2)

        axis(2,at=c(0,max(max(unlist(lapply(pdf,function(x) return(max(c(0,x),na.rm=TRUE)))),na.rm=TRUE),axTicks(2))),labels=c("",""),tck=0,lwd=0.5,pos=XRANGE[[v]][1])
        axis(2,lwd.ticks=0.5,lwd=0,pos=XRANGE[[v]][1],tck=-0.03,cex.axis=TEXT_SIZE)
        axis(1,at=XRANGE[[v]],labels=c("",""),tck=0,lwd=0.5,pos=0)
        axis(1,lwd.ticks=0.5,lwd=0,pos=0,tck=-0.03,cex.axis=TEXT_SIZE)

        mtext("Density of probability",side=2,line=0.9,adj=0.5667)
        mtext(paste(VARIABLES_NAMES[[v]][5],sep=""),side=1,line=0.8,adj=0.5667)
    }


    addImg <- function( obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
                        x = NULL, # mid x coordinate for image
                        y = NULL, # mid y coordinate for image
                        width = NULL, # width of image (in x coordinate units)
                        interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing.
                       ){
        if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
        USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
        PIN <- par()$pin # The current plot dimensions, (width, height), in inches
        DIM <- dim(obj) # number of x-y pixels for the image
        ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
        WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
        HEIi <- WIDi * ARp # height in inches
        HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
        rasterImage(image = obj,
            xleft = x-(width/2), xright = x+(width/2),
            ybottom = y-(HEIu/2), ytop = y+(HEIu/2),
            interpolate = interpolate)
      }


      addPageNumber <- function(PAGE_NB) {
          par(mar=c(0,0,0,0))
          plot(0,0,type='n', xlim=c(0,1), ylim=c(0, 1), frame=FALSE, axes=FALSE)
          if(gtools::even(PAGE_NB)) {
              text(1,1, PAGE_NB, adj=c(1,1), cex=1.5)
          } else {
              text(0,1, PAGE_NB, adj=c(0,1), cex=1.5)
          }
          return(PAGE_NB+1)
      }


      new_page <- function(line_count, max_lines, page_nb=PAGE_NB, step=1) {
          if (line_count >= max_lines) {
              layout(matrix(c(3,3,3,3,2,3,1+2*(PAGE_NB%%2),3,3-2*(PAGE_NB%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
              par(mar=c(0,0,0,0))
              page_nb = addPageNumber(page_nb)
              plot(0,0, type='n', xlim=c(0.01,1), ylim=c(max_lines,1), axes=FALSE, frame=FALSE)
              return(c(0, page_nb))
          }
          return(c(line_count + step, page_nb))
      }

      addEmptyPage <- function(page_nb) {
          layout(matrix(c(3,3,3,3,2,3,1+2*(page_nb%%2),3,3-2*(page_nb%%2)), ncol=3, byrow=TRUE), height=c(1.5, 26.7, 1.5), width=c(1.5,18,1.5))
          page_nb = addPageNumber(page_nb)
          plot(0,type='n',xlim=c(0,1),ylim=c(0,1),frame=FALSE,axes=FALSE,xlab="",ylab="")
          return(page_nb)
      }

}
