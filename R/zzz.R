.onAttach <- function(libname, pkgname) {
  required_packages <- c("cpc", "geigen", "ggplot2", "mvtnorm", "numDeriv", "pracma", "RRembo", "tidyr", "dplyr", "alphahull")

  missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
  if(length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  }
}
