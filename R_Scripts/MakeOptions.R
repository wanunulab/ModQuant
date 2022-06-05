#### ====================================== ####
####    Make/define the main code options   ####
#### ====================================== ####
wdir = getwd( ) 

args = commandArgs(trailingOnly = FALSE)
file.path = strsplit(args[4], "=")
#print(file.path)
#print(file.path[[1]][2])
file.path = strsplit(file.path[[1]][2], split="FeatureExtract.R")

if (file.path[[1]] == ""){
  source("R_Scripts/HandlePackages.R")
} else{
  source(paste(file.path[[1]], "/R_Scripts/HandlePackages.R", sep=""))
}

#source(paste(wdir,"/R/HandlePackages.R",sep = ""))
#source("HandlePackages.R")
#source(paste(file.path[[1]], "/R_Scripts/HandlePackages.R", sep=""))


option_list = list(
  make_option(c("-b", "--bam_file"), type="character", default="", 
              help="path to bam file with aligned reads (native or direct), Ex: ...", metavar="character"),
  make_option(c("-s", "--summary_file"), type="character", default="", 
              help="nanopolish summary file, Ex: ...", metavar="character"),
  make_option(c("-e", "--eventalign_file"), type="character", default="", 
              help="nanopolish event align file with raw signal data, Ex: ...", metavar="character"),
  
  make_option(c("-m", "--minimum_mapping_score"), type="character", default=NULL, 
              help="minimum mapping quality score for each read, Ex: ...", metavar="character"),
  make_option(c("-l", "--minimum_read_length"), type="character", default=NULL, 
              help="minimum length (nt) for each read, Ex: ...", metavar="character"),
  
  make_option(c("-t", "--target_position"), type="character", default="", 
              help="target base of interest, Ex: ...", metavar="character"),
  make_option(c("-r", "--target_region"), type="character", default="15", 
              help="region spanning target base, Ex: ...", metavar="character"),
  make_option(c("-p", "--preserved_correct_calls"), type="character", default="12", 
              help="minimum correct basecalls required inorder to parse features, Ex: ...", metavar="character"),
  make_option(c("-n", "--number_of_bases_and_5mers"), type="character", default="5", 
              help="number of basecalls and 5mer signal frames to extract, Ex: ...", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="Mod_Features.csv", 
              help="output file name of ML training dataframe, Ex: ...", metavar="character"))







#features.df = feature_compile_func(bam.input = filtered.bam.df, summary.file = summary, event.align.file = raw.signals.df, 
#                                   position.of.interest = 511, num.region.of.interest = 15, region.error.limit = 12,
#                                   num.bases.and.5mers = 5)
    
  #make_option(c("-p", "--position"), type="character", default=NULL, 
  #            help="position", metavar="character"),
  #make_option(c("-c", "--chromosome"), type="character", default=NULL, 
  #            help="chromosome", metavar="character"),
  #make_option(c("-o", "--out"), type="character", default="out.pdf", 
  #            help="output file name of signal plot", metavar="character"))