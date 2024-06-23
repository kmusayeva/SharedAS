#Install required packages
required_packages <- c("geigen", "ggplot2", "mvtnorm", "numDeriv", "pracma", "RRembo", "tidyr", "dplyr", "alphahull")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
print(missing_packages)
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Install the package cpc
install.packages("remotes")
remotes::install_github("tpepler/cpc")
remotes::install_github("mbinois/RRembo")

# Load devtools
library(remotes)

# Path to the local package directory
package_path <- "./"

# Install the package from the local directory
install_local(package_path)

