## Figure 3: Presentation of the reconstruction
##
## Loading necessary data
load(url('https://github.com/mchevalier2/Papers/raw/master/Chevalier_etal_MD962048/figures/Fig3.RData'))

OUTPUT_FOLDER=getwd()
s <- readline(prompt=paste0("Where should the figure be saved?\nDefault is current workin directory (",OUTPUT_FOLDER,"): "))
if(s != '') OUTPUT_FOLDER <- s

pkg2install=c()
if (! ("sp" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'sp')
if (! ("raster" %in% rownames(installed.packages()))) pkg2install=c(pkg2install, 'raster')


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
