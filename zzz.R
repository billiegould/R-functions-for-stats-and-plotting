.onLoad <- function(libname, pkgname) {
    
    ## find and set the DEEPSEA_HOME
    deepseahome <- Sys.getenv('DEEPSEA_HOME')
    if (deepseahome == '') {
        deepseahome <- dirname(Sys.which('lbwf'))
    }
    if (!grepl('bin$', deepseahome) && deepseahome != '') 
        deepseahome <- file.path(deepseahome, 'bin')
    
    ## find and set the LBWFDATA (pipeline related data)
    lbwfdata <- Sys.getenv("LBWFDATA")
    if (!file.exists(lbwfdata)) {
        lbwfdata <- '/prednet/data03/deepseadata/data_1.7.0'
        if (!file.exists(lbwfdata)) {
            lbwfdata <- '/opt/ngsTools/share'
        }
    }
    if (!file.exists(lbwfdata)) stop('lbwfdata is not defined!')
    
    targetInfoDir <- file.path(lbwfdata, "panels")
    mrdPanelDir <- file.path(lbwfdata, "MRD_panels")
    options(deepseahome=deepseahome)
    options(lbwfdata=lbwfdata)
    options(targetInfoDir=targetInfoDir)
    options(mrdPanelDir=mrdPanelDir)
    options(outputBySampleDir='/prednet/data02/OutputBySample')
    options(outputByProjectDir='/prednet/data02/OutputByProject')
    options(runRootDir='/prednet/data10/OutputByRun10')
    options(CLIA_results='/prednet/CLIA_results')
    
    ## VEP host IP:
    vephost <- '10.10.0.63'
    options(vephost=vephost)
    
    invisible()
}
