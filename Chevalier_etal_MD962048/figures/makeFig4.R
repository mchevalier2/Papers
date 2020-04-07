## Figure 3: Presentation of the reconstruction
##

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


}
