#### ====================================== ####
####     Install/load required packages     ####
#### ====================================== ####

#  Check required packages  ...
#  .. and install the required/not installed ones. 

list.of.packages <- 
  c("ggplot2", "data.table", "optparse", "dplyr", "BiocManager")

new.packages <- 
  list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

if (!("GenomicAlignments" %in% installed.packages()[,"Package"])){
  BiocManager::install("GenomicAlignments")
}

if (!("Rsamtools" %in% installed.packages()[,"Package"])){
  BiocManager::install("Rsamtools")
}

if (!("Biobase" %in% installed.packages()[,"Package"])){
  BiocManager::install("Biobase")
}

#  Quietly Load the required packages
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(GenomicAlignments)))
suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(Biobase)))
