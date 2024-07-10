#!/usr/bin/env Rscript

# Function to check and install packages
install_if_missing = function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of required packages
required_packages = c("data.table", "stringr", "dplyr", "tidyr",
                      "ggplot2", "vcfR")

# Check and install missing packages
install_if_missing(required_packages)
