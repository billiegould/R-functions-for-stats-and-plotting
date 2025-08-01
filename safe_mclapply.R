## hack for getting mclapply to return all warning messages encounterd
## from https://stackoverflow.com/questions/21486658/warnings-suppressed-with-mclapply-in-r 7/24/2023

safe_mclapply <- function(X, FUN, mc.cores, stop.on.error=T, ...){
  fun <- function(x){
    res_inner <- tryCatch({
      withCallingHandlers(
        expr = {
          FUN(x, ...)
        }, 
        warning = function(e) {
          message_parallel(trimws(paste0("WARNING [element ", x,"]: ", e)))
          # this line is required to continue FUN execution after the warning
          invokeRestart("muffleWarning")
        },
        error = function(e) {
          message_parallel(trimws(paste0("ERROR [element ", x,"]: ", e)))
        }
      )},
      error = function(e){
        # error is returned gracefully; other results of this core won't be affected
        return(e)
      }
    )
    return(res_inner)
  }
  
  res <- mclapply(X, fun, mc.cores=mc.cores)
  failed <- sapply(res, inherits, what = "error")
  if (any(failed == T)){
    error_indices <- paste0(which(failed == T), collapse=", ")
    error_traces <- paste0(lapply(res[which(failed == T)], function(x) x$message), collapse="\n\n")
    error_message <- sprintf("Elements with following indices failed with an error: %s. Error messages: \n\n%s", 
                             error_indices,
                             error_traces)
    if (stop.on.error)
      stop(error_message)
    else
      warning(error_message, "\n\n### Errors will be ignored ###")
  }
  return(res[!failed])
}

#' Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}