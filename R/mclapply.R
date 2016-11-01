mclapply.or.stop <- function
### Run mclapply but stop with an error instead of silently returning
### a try-error.
(...
### passed to parallel::mclapply
){
  result.list <- mclapply(...)
  is.error <- sapply(result.list, inherits, "try-error")
  if(any(is.error)){
    print(result.list[is.error])
    stop("errors in mclapply")
  }
  result.list
### List of values, same as lapply.
}

### Set mc.cores option from an environment variable.
PPN.cores <- function
(variable="PBS_NUM_PPN"
### The PBS_NUM_PPN variable is defined by PBS,
### for example it will be 4 when qsub -l nodes=1:ppn=4
){
  stopifnot(is.character(variable))
  stopifnot(length(variable) == 1)
  ppn.txt <- Sys.getenv(variable)
  ppn <- as.integer(ppn.txt)
  if(is.finite(ppn)){
    options(mc.cores=ppn)
  }
  ppn
### The new value of options(mc.cores).
}

