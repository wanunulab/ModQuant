#### ========================================== ####
####           Main: Feature Extraction Prep    ####
#### ========================================== ####
#!/usr/bin/env Rscript

start.time = Sys.time()

wdir = getwd( ) 
print(wdir)

#input = file('stdin', 'r')
#row = readLines(input, n=1)
#print(row)##

#options(echo=TRUE)
args = commandArgs(trailingOnly = FALSE)
file.path = strsplit(args[4], "=")
#print(file.path)
#print(file.path[[1]][2])
file.path = strsplit(file.path[[1]][2], split="FeatureExtract.R")
#print(file.path)

cat("Initializing feature extraction package...")

if (file.path[[1]] == ""){
  source("R_Scripts/MakeOptions.R")
} else{
  source(paste(file.path[[1]], "/R_Scripts/MakeOptions.R", sep=""))
}
#source(paste(wdir,"/R/MakeOptions.R",sep = ""))
#source(paste(file.path[[1]], "/R_Scripts/MakeOptions.R", sep=""))
#source("R_Scripts/MakeOptions.R")

if (file.path[[1]] == ""){
  source("R_Scripts/Feature_Preprocessing_Functions.R")
} else{
  source(paste(file.path[[1]], "/R_Scripts/Feature_Preprocessing_Functions.R", sep=""))
}

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



#source(paste(wdir,"/R/Feature_Preprocessing_Functions.R",sep = ""))
#source(psate(file.path[[1]], "/R_Scripts/Feature_Preprocessing_Functions.R"))


#########Loading bam file########
#D1b = as.character(opt$bam_file) #example path
#allD1 = Rsamtools::scanBam(file = D1b);
#bam.df = as.data.frame(allD1)
#################################

cat("\nLoading Data...")

#########Loading bam file########
D1b = as.character(opt$bam_file) #example path
allD1 = Rsamtools::scanBam(file = D1b);
bam.df = as.data.frame(allD1)
#################################

#######Loading eventalign summary file########
summary = read.table(as.character(opt$summary_file), header = F, skip = 1)
##############################################

#######Loading eventalign data file########
raw.signals.df = read.table(as.character(opt$eventalign_file), header = T)
###########################################

###########Loading read input filters###############
minimum.read.length = as.numeric(opt$minimum_read_length)
minimum.mapping.score = as.numeric(opt$minimum_mapping_score)
#cat(sprintf("\nMinimum read length filter: %0.0f nts", minimum.read.length))
#cat(sprintf("\nMinimum mapping score: %0.0f ", minimum.mapping.score))

if (is.na(minimum.read.length)){
  minimum.read.length = NULL
}

if (is.na(minimum.mapping.score)){
  minimum.mapping.score = NULL
}
######################################################

#########Loading feature extraction conditions###############
position.interest = as.numeric(opt$target_position)
num.region.interest = as.numeric(opt$target_region)
region.limit = as.numeric(opt$preserved_correct_calls)
num.bases.5mers = as.numeric(opt$number_of_bases_and_5mers)
#############################################################

filtered.bam.df = read_filter_func(whole.bam.df = bam.df, 
                                   min.mapping.qual.score = minimum.mapping.score, 
                                   min.read.length = minimum.read.length) #50, 565

cat("\nExtracting Features...")

features.df = feature_compile_func(bam.input = filtered.bam.df, 
                                   summary.file = summary, 
                                   event.align.file = raw.signals.df, 
                                   position.of.interest = position.interest, 
                                   num.region.of.interest = num.region.interest, 
                                   region.error.limit = region.limit,
                                   num.bases.and.5mers = num.bases.5mers)


fwrite(features.df, opt$out, na = 'NA')
cat("\nFeature Extraction is done! Output file is ready!")

end.time = Sys.time()
processing.time = difftime(end.time, start.time, units = 'secs')
cat(sprintf("\nTotal run time: %0.3f seconds", processing.time))


