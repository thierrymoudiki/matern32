
# global reference to packages (will be initialized in .onLoad)
jax <- NULL
jaxlib <- NULL

install_miniconda_ <- function(silent = TRUE)
{
  try(reticulate::install_miniconda(),
      silent = silent)
}


install_packages <- function(pip=TRUE) {
    reticulate::py_install("jax", pip = pip)
    reticulate::py_install("jaxlib", pip = pip)
}


.onLoad <- function(libname, pkgname) {
  
  if (.Platform$OS.type == "unix")
  {
    do.call("install_miniconda_", list(silent=TRUE))
    do.call("install_packages", list(pip=TRUE))
    
    jax <<- reticulate::import("jax", delay_load = TRUE)
    jaxlib <<- reticulate::import("jaxlib", delay_load = TRUE)
  }
  
}
