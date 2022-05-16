# install libraries
if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicAlignments")
install.packages("ggplot2")
install.packages("data.table")
#install.packages("Rsamtools")

# Load Necessary Libraries
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)
library(data.table)

# the function
Decipher_cigar_func <- function(cigar, start.pos, sequence, base.quality = ''){
  cig.df  = cigarToRleList(cigar)[[1]] 
  seq.temp =strsplit(as.character(sequence),"")[[1]]
  phred.q = base.quality 
  phred.q <- PhredQuality(phred.q)
  phred.q <- as(phred.q, "IntegerList")  # quality scores
  phred.q.lst = phred.q[[1]]
  
  values = cig.df@values
  cig.length = cig.df@lengths
  
  positions.list = c() 
  nucleutide.list = c()
  base.qual.list = c()
  
  counter.on.sequence = 1
  last.avail.pos = start.pos - 1
  # loop over all the cigar values, decifer the read
  for (j in c(1:length(values))){
    value = values[j]
    cig.length.tmp = cig.length[j]
    # #M	Match 		Exact match of x positions
    if (value == "M"){
      positions.list = c(positions.list,c((last.avail.pos+1):(last.avail.pos+cig.length.tmp)))
      nucleutide.list = c(nucleutide.list,seq.temp[c(counter.on.sequence:(counter.on.sequence+cig.length.tmp-1))])
      base.qual.list = c(base.qual.list,phred.q.lst[c(counter.on.sequence:(counter.on.sequence+cig.length.tmp-1))])
      
      counter.on.sequence = counter.on.sequence + cig.length.tmp
      last.avail.pos = last.avail.pos + cig.length.tmp
    }
    #D	Deletion 	Next x positions on ref don't match
    if (value == "D"){
      positions.list = c(positions.list,c((last.avail.pos+1):(last.avail.pos+cig.length.tmp)))
      nucleutide.list = c(nucleutide.list,rep("-",cig.length.tmp))
      base.qual.list = c(base.qual.list,rep(-1,cig.length.tmp))
      
      last.avail.pos = last.avail.pos + cig.length.tmp
    }
    #I	Insertion 	Next x positions on query don't match
    if (value == "I"){
      positions.list = c(positions.list,rep((last.avail.pos),cig.length.tmp))
      nucleutide.list = c(nucleutide.list,rep("+",cig.length.tmp))
      base.qual.list = c(base.qual.list,rep(-1,cig.length.tmp))
      
      counter.on.sequence = counter.on.sequence + cig.length.tmp
      #last.avail.pos = last.avail.pos + cig.length.tmp
    }
    #N	Alignment gap 	Next x positions on ref don't match
    if (value == "N"){
      last.avail.pos = last.avail.pos + cig.length.tmp
    }
    #S	Alignment gap 	Next x positions on ref don't match
    if (value == "S"){
      counter.on.sequence = counter.on.sequence + cig.length.tmp
    }
  } # end of looping over the cigar values
  pos.base.df = data.frame("pos"=positions.list,"base" = nucleutide.list, "base_quality"=base.qual.list)
  #############
  return(pos.base.df)
}


#### read bam file
# Path to the bam file
#D1b = "~/repositories/2020/nanopore/DataSet/Hela/bamFiles/Hela_Dir2_sorted.bam"

####100 percent modified data experiment 1####
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment1\\Analyze_BAM\\MCM5.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment1\\Analyze_BAM\\MRPS14.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment1\\Analyze_BAM\\PRPSAP1.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment1\\Analyze_BAM\\PSMB2.bam"


####100 percent modified data experiment 2####
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment2\\Analyze_BAM\\MCM5_updated_sorted.bam" #This is the synthetic transcript bam file
D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment2\\Analyze_BAM\\MRPS14_updated_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment2\\Analyze_BAM\\PRPSAP1_updated_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment2\\Analyze_BAM\\PSMB2_updated_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment2\\Analyze_BAM\\PTTG1IP_subset_updated_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\100_perc_mod\\syn_100_experiment2\\Analyze_BAM\\RNF7_subset_updated_sorted.bam"


####0 percent Unmodified data experiment 2####
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\0_perc_mod\\syn_0_experiment2\\Analyze_BAM\\MCM5_updated_small_sorted.bam" #This is the synthetic transcript bam file
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\0_perc_mod\\syn_0_experiment2\\Analyze_BAM\\MRPS14_updated_small_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\0_perc_mod\\syn_0_experiment2\\Analyze_BAM\\PRPSAP1_updated_small_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\0_perc_mod\\syn_0_experiment2\\Analyze_BAM\\PSMB2_updated_small_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\0_perc_mod\\syn_0_experiment2\\Analyze_BAM\\PTTG1IP_updated_small_sorted.bam"
#D1b = "C:\\Users\\amrma\\OneDrive\\Desktop\\Psi_ML\\0_perc_mod\\syn_0_experiment2\\Analyze_BAM\\RNF7_updated_small_sorted.bam"


####0 percent modified data####
#D1b = "Z:/MCM5_small_sorted.bam"
#D1b = "Z:/PRPSAP1_small_sorted.bam"
#D1b = "Z:/PSMB2_small_sorted.bam"
#D1b = "Z:/MRPS14_small_sorted.bam"

# read bam file
allD1 = Rsamtools::scanBam(file = D1b);

# convert bam file to a data frame
bam.df = as.data.frame(allD1)

# load the summary file
#summary = read.table("Z:/MRPS14_summary",header = F,skip = 1)
#summary = read.table("Z:/PSMB2_summary",header = F,skip = 1)
#summary = read.table("Z:/PRPSAP1_summary",header = F,skip = 1)
#summary = read.table("Z:/MCM5_summary",header = F,skip = 1)

#summary = read.table("MCM5_summary",header = F,skip = 1)
#summary = read.table("MRPS14_summary",header = F,skip = 1)
#summary = read.table("PRPSAP1_summary",header = F,skip = 1)
#summary = read.table("PSMB2_summary",header = F,skip = 1)

###Summary for 100% synthetic experiment 1###
#summary = read.table("Summary_Files\\MCM5_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\MRPS14_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\PRPSAP1_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\PSMB2_summary",header = F,skip = 1)

###Summary for 100% synthetic experiment 2###
#summary = read.table("Summary_Files\\MCM5_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\MRPS14_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\PRPSAP1_summary",header = F,skip = 1, nrows=10000000)
#summary = read.table("Summary_Files\\PSMB2_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\PTTG1IP_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\RNF7_summary",header = F,skip = 1)

###Summary for 0% synthetic experiment 2###
#summary = read.table("Summary_Files\\MCM5_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\MRPS14_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\PRPSAP1_summary",header = F,skip = 1)
summary = read.table("Summary_Files\\PSMB2_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\PTTG1IP_summary",header = F,skip = 1)
#summary = read.table("Summary_Files\\RNF7_summary",header = F,skip = 1)

print(nrow(summary))
View(summary)
#print(summary[summary$V2 == "81a08372-83ee-44ca-aa0d-370acf8ec8a9", ])


###load the txt file that contains the raw signals####
#raw.signals.df = read.table("Event_Aligned_Raw\\MCM5.txt",header = T)
#raw.signals.df = read.table("Event_Aligned_Raw\\MRPS14.txt",header = T)
#raw.signals.df = read.table("Event_Aligned_Raw\\PRPSAP1.txt",header = T, nrows=15000000)
raw.signals.df = read.table("Event_Aligned_Raw\\PSMB2.txt",header = T, nrows=15000000)

#raw.signals.df = read.table("Event_Aligned_Raw\\PTTG1IP.txt",header = T, nrows=15000000)
#raw.signals.df = read.table("Event_Aligned_Raw\\RNF7.txt",header = T, nrows=15000000)

#raw.signals.df = read.table("Z:/MCM5.txt",header = T)
#raw.signals.df = read.table("MRPS14.txt",header = T, nrows=8000000)
#raw.signals.df = read.table("PRPSAP1.txt",header = T, nrows=8000000)
#raw.signals.df = read.table("PSMB2.txt",header = T, nrows=8000000)
head(raw.signals.df)
print(nrow(raw.signals.df))
#######Looping through bam file and extracting metadata######
print(length(bam.df))
print(nrow(bam.df))
r = 0
dr = 0
for(row in 1:nrow(bam.df)){
  r = r+1
  if (nchar(bam.df$seq[row])>700){
    dr = dr+1
  }
}
print(r)
print(dr)
######

#Filtering bam.df to only contain reads of interest
nrow(bam.df[nchar(bam.df$seq)>700,])
bam.df = bam.df[nchar(bam.df$seq)>800,]
nrow(bam.df)
head(bam.df)
bam.df = bam.df[nchar(bam.df$seq)>=545 & bam.df$mapq>50,] #condition with length and mapping quality
nrow(bam.df)
print(bam.df$mapq)


# look at a specific read (alignment)
row = 103 #406 is a good row
print(bam.df[row, ])
bam.df$qname[row]
bam.df$rname[row]
bam.df$strand[row]
bam.df$mapq[row]
bam.df$cigar[row]
bam.df$seq[row]
bam.df$qual[row]
nchar(bam.df$seq[row])
print(bam.df$seq[row])
bam.df$pos[row] #There seems to be an offset in the function, missing first and last nucleotides (row=406)  

####This code checks to make sure the meta dataframe contains the correct cigar info
x = bam.df[bam.df$qname == "b1327ab8-b5a1-4304-8044-abe259fe168b", ]
print(x$qname)
print(x)
print(bam.df$qname[row])
####################
print(which(bam.df$qname == "d980357e-b00f-4030-92fe-8fc9df1301b9"))
row = 24
bam.df$qname[row]
# how to use the decode function:
seq.df = Decipher_cigar_func(cigar = bam.df$cigar[row],start.pos = bam.df$pos[row],
                             sequence = bam.df$seq[row], base.quality = bam.df$qual[row])

View(seq.df)
print(nrow(seq.df))

##########################################
#Looping through 15mer to make sure most of the base calls are preserved
#############
head(seq.df)
target.positions = seq(504, 518, 1)
test.df = seq.df[seq.df$pos %in% target.positions,]

target.5mer = seq(509, 513, 1)
test.df = seq.df[seq.df$pos %in% target.5mer,]
test.df = test.df[test.df$base != "+" & test.df$base != "-",]
print(test.df)
#############


previous.pos = 0
counter.bases = 0
for (j in 1:nrow(test.df)){
  #print(test.df[j, ])
  current.base = test.df$base[j]
  current.pos = test.df$pos[j]
  if (current.pos != previous.pos){
    previous.pos = test.df$pos[j]
    if (current.base == "A" | current.base == "C" | current.base == "G" | current.base == "T"){
      counter.bases = counter.bases + 1
    }
  }
}
###############################################

##This function filters reads on how many bases they have in the 15mer region##
synthetic_filter_func <- function(current.read.df, target.positions){
  sliced.read.df = current.read.df[current.read.df$pos %in% target.positions,]
  previous.pos = 0
  counter.bases = 0
  for (j in 1:nrow(sliced.read.df)){
    current.base = sliced.read.df$base[j]
    current.pos = sliced.read.df$pos[j]
    if (current.pos != previous.pos){
      previous.pos = sliced.read.df$pos[j]
      if (current.base == "A" | current.base == "C" | current.base == "G" | current.base == "T"){
        counter.bases = counter.bases + 1
      }
    }
  }
  return(counter.bases)
}




##This function filters reads on how many bases they have in the 15mer region##
synthetic_filter_func <- function(current.read.df, target.positions){
  sliced.read.df = current.read.df[current.read.df$pos %in% target.positions,]
  #print(sliced.read.df)
  #sprintf("Sliced out targert length: %s", nrow(sliced.read.df))
  #if (is.na(sliced.read.df)){
  if (nrow(sliced.read.df) == 0){
    #print("Read does not contain target region")
    counter.bases = 0
  }
  else{
    #print(sliced.read.df)
    previous.pos = 0
    counter.bases = 0
    for (j in 1:nrow(sliced.read.df)){
      current.base = sliced.read.df$base[j]
      current.pos = sliced.read.df$pos[j]
      if (current.pos != previous.pos){
        previous.pos = sliced.read.df$pos[j]
        if (current.base == "A" | current.base == "C" | current.base == "G" | current.base == "T"){
          counter.bases = counter.bases + 1
        }
      }
    }
  }
  return(counter.bases)
}
############################################
##This function extracts and extends the current mean, standard deviation, and raw sample points for a 5mer from nanopolish##
raw_trace_meta_extract_func <- function(read.raw.signal.df, kmer.pos){
  read.raw.signal.df.pos = read.raw.signal.df[which(read.raw.signal.df$position == kmer.pos),]
  if(nrow(read.raw.signal.df.pos) >= 1){
    #print('HIT')
    samples_string = c()
    samples = c()
    for (i in c(1:nrow(read.raw.signal.df.pos))){
      # splitting the string of floats
      #print(nchar(read.raw.signal.df.pos$samples[i]))
      splitted.string = strsplit(read.raw.signal.df.pos$samples[i],",")[[1]]
      # convert splitted strings to numeric values
      samples.tmp = as.numeric(splitted.string)
      # append the current row's signals to the main list
      samples = c(samples, samples.tmp)
      samples_string = c(samples_string, read.raw.signal.df.pos$samples[i])
    }
    samples_string = paste(samples_string, collapse = ",")
    signal.mean = mean(samples) #full 5mer sample mean
    signal.sd   = sd(samples) #full 5mer sample std
  } else {
    #print('No raw current')
    samples_string = NA
    signal.mean = NA #full 5mer sample mean
    signal.sd   = NA #full 5mer sample std
  }
  #signal.meta = c(signal.mean, signal.sd, samples_string)
  signal.meta.df = data.frame("kmer_mean"=signal.mean,"kmer_std"=signal.sd,"raw_current"=samples_string)
  #############
  return(signal.meta.df)
}

#############################################
row = 600
nrow(bam.df)
bam.df = bam.df[nchar(bam.df$seq)>565 & bam.df$mapq>50,] #565 condition with length and mapping quality
nrow(bam.df)
bam.df = bam.df[bam.df$qname != "81a08372-83ee-44ca-aa0d-370acf8ec8a9", ]
print(bam.df$mapq)
bam.df$qname[row]

#####################################
target.positions = seq(504, 518, 1)
target.5mer = seq(509, 513, 1)
fifteenMer.error.limit = 12

#read index and read ID in summary file
read.IDs = c()
read.indices = c()
#Base Calls at 509, 510, 511(Psi), 512, 513
bn2 = c()
bn1 = c()
b0 = c()
bp1 = c()
bp2 = c()
#Quality Score at 509, 510, 511(Psi), 512, 513
qn2 = c()
qn1 = c()
q0 = c()
qp1 = c()
qp2 = c()
######################################################################
#Time-domain features of target kmer with psi at position n2
signal.means.n2 = c()
signal.stds.n2 = c()
#raw.current.samples.n2 = c()
string.raw.current.samples.n2 = c()

#Time-domain features of target kmer with psi at position n1
signal.means.n1 = c()
signal.stds.n1 = c()
#raw.current.samples.n1 = c()
string.raw.current.samples.n1 = c()

#Time-domain features of target kmer with psi at position 0
signal.means.0 = c()
signal.stds.0 = c()
#raw.current.samples.0 = c()
string.raw.current.samples.0 = c()

#Time-domain features of target kmer with psi at position p1
signal.means.p1 = c()
signal.stds.p1 = c()
#raw.current.samples.p1 = c()
string.raw.current.samples.p1 = c()

#Time-domain features of target kmer with psi at position p2
signal.means.p2 = c()
signal.stds.p2 = c()
#raw.current.samples.p2 = c()
string.raw.current.samples.p2 = c()
######################################################################

######################################################################
#Time-domain features of kmer 12 bases upstream with psi at position n2
upstream.signal.means.n2 = c()
upstream.signal.stds.n2 = c()
#upstream.raw.current.samples.n2 = c()
upstream.string.raw.current.samples.n2 = c()

#Time-domain features of kmer 12 bases upstream with psi at position n1
upstream.signal.means.n1 = c()
upstream.signal.stds.n1 = c()
#upstream.raw.current.samples.n1 = c()
upstream.string.raw.current.samples.n1 = c()

#Time-domain features of kmer 12 bases upstream with psi at position 0
upstream.signal.means.0 = c()
upstream.signal.stds.0 = c()
#upstream.raw.current.samples.0 = c()
upstream.string.raw.current.samples.0 = c()

#Time-domain features of kmer 12 bases upstream with psi at position p1
upstream.signal.means.p1 = c()
upstream.signal.stds.p1 = c()
#upstream.raw.current.samples.p1 = c()
upstream.string.raw.current.samples.p1 = c()

#Time-domain features of kmer 12 bases upstream with psi at position p2
upstream.signal.means.p2 = c()
upstream.signal.stds.p2 = c()
#upstream.raw.current.samples.p2 = c()
upstream.string.raw.current.samples.p2 = c()
################################################################################
##################1021 position is 3' and 0 position is 5'######################
pos.target.n2 = 510 #kmer with psi in -2 position in nanopolish 
pos.target.n1 = 509 #kmer with psi in -1 position in nanopolish 
pos.target.0 = 508 #kmer with psi in 0 position in nanopolish 
pos.target.p1 = 507 #kmer with psi in +1 position in nanopolish 
pos.target.p2 = 506 #kmer with psi in +2 position in nanopolish 

pos.upstream.n2 = 522 #kmer 12 bases upstream from psi in position n2 in nanopolish
pos.upstream.n1 = 521 #kmer 12 bases upstream from psi in position n1 in nanopolish
pos.upstream.0 = 520 #kmer 12 bases upstream from psi in position 0 in nanopolish
pos.upstream.p1 = 519 #kmer 12 bases upstream from psi in position p1 in nanopolish
pos.upstream.p2 = 518 #kmer 12 bases upstream from psi in position p2 in nanopolish

######################################################################

######################################################################
n=8000
#for (j in 1:nrow(bam.df)){
for (j in 1:n){
  row = j
  seq.df = Decipher_cigar_func(cigar = bam.df$cigar[row],start.pos = bam.df$pos[row],
                               sequence = bam.df$seq[row], base.quality = bam.df$qual[row])
  preserved.15mer = synthetic_filter_func(current.read.df = seq.df, 
                                          target.positions = target.positions)
  if (preserved.15mer >= fifteenMer.error.limit){
    read.name = bam.df$qname[row]
    # find the read index 
    read.indx = summary$V1[which(summary$V2 == read.name)]
    print(read.indx)
    read.indx = unique(read.indx)
    # all the rows of the nanopolish output that correspond to our targeted read
    target.rows = which(raw.signals.df$read_index == read.indx)###issue
    #head(target.rows)
    # make a new data frame that has only our read of interest
    read.raw.signal.df = raw.signals.df[target.rows,]
################################################################################################  
    meta.current.target.n2 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                         kmer.pos = pos.target.n2)
    samples_string = paste(meta.current.target.n2$raw_current, 
                           collapse = ",")
    string.raw.current.samples.n2 = c(string.raw.current.samples.n2, 
                                      samples_string)
    signal.means.n2 = c(signal.means.n2, meta.current.target.n2$kmer_mean) #append mean to full list
    signal.stds.n2 = c(signal.stds.n2, meta.current.target.n2$kmer_std) #append std to full list
    
    
    
    meta.current.target.n1 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                        kmer.pos = pos.target.n1)
    samples_string = paste(meta.current.target.n1$raw_current, 
                           collapse = ",")
    string.raw.current.samples.n1 = c(string.raw.current.samples.n1, 
                                     samples_string)
    signal.means.n1 = c(signal.means.n1, meta.current.target.n1$kmer_mean) #append mean to full list
    signal.stds.n1 = c(signal.stds.n1, meta.current.target.n1$kmer_std) #append std to full list
    
    
    meta.current.target.0 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                    kmer.pos = pos.target.0)
    samples_string = paste(meta.current.target.0$raw_current, 
                           collapse = ",")
    string.raw.current.samples.0 = c(string.raw.current.samples.0, 
                                   samples_string)
    signal.means.0 = c(signal.means.0, meta.current.target.0$kmer_mean) #append mean to full list
    signal.stds.0 = c(signal.stds.0, meta.current.target.0$kmer_std) #append std to full list
    
    
    meta.current.target.p1 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                         kmer.pos = pos.target.p1)
    samples_string = paste(meta.current.target.p1$raw_current, 
                           collapse = ",")
    string.raw.current.samples.p1 = c(string.raw.current.samples.p1, 
                                      samples_string)
    signal.means.p1 = c(signal.means.p1, meta.current.target.p1$kmer_mean) #append mean to full list
    signal.stds.p1 = c(signal.stds.p1, meta.current.target.p1$kmer_std) #append std to full list
    
    
    
    meta.current.target.p2 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                         kmer.pos = pos.target.p2)
    samples_string = paste(meta.current.target.p2$raw_current, 
                           collapse = ",")
    string.raw.current.samples.p2 = c(string.raw.current.samples.p2, 
                                      samples_string)
    signal.means.p2 = c(signal.means.p2, meta.current.target.p2$kmer_mean) #append mean to full list
    signal.stds.p2 = c(signal.stds.p2, meta.current.target.p2$kmer_std) #append std to full list 
###############################################################################################
    upstream_meta.current.target.n2 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                         kmer.pos = pos.upstream.n2)
    samples_string = paste(upstream_meta.current.target.n2$raw_current, 
                           collapse = ",")
    upstream.string.raw.current.samples.n2 = c(upstream.string.raw.current.samples.n2, 
                                      samples_string)
    upstream.signal.means.n2 = c(upstream.signal.means.n2, upstream_meta.current.target.n2$kmer_mean) #append mean to full list
    upstream.signal.stds.n2 = c(upstream.signal.stds.n2, upstream_meta.current.target.n2$kmer_std) #append std to full list
    
    
    upstream_meta.current.target.n1 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                                  kmer.pos = pos.upstream.n1)
    samples_string = paste(upstream_meta.current.target.n1$raw_current, 
                           collapse = ",")
    upstream.string.raw.current.samples.n1 = c(upstream.string.raw.current.samples.n1, 
                                               samples_string)
    upstream.signal.means.n1 = c(upstream.signal.means.n1, upstream_meta.current.target.n1$kmer_mean) #append mean to full list
    upstream.signal.stds.n1 = c(upstream.signal.stds.n1, upstream_meta.current.target.n1$kmer_std) #append std to full list
    
    
    upstream_meta.current.target.0 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                                  kmer.pos = pos.upstream.0)
    samples_string = paste(upstream_meta.current.target.0$raw_current, 
                           collapse = ",")
    upstream.string.raw.current.samples.0 = c(upstream.string.raw.current.samples.0, 
                                               samples_string)
    upstream.signal.means.0 = c(upstream.signal.means.0, upstream_meta.current.target.0$kmer_mean) #append mean to full list
    upstream.signal.stds.0 = c(upstream.signal.stds.0, upstream_meta.current.target.0$kmer_std) #append std to full list
    
    
    upstream_meta.current.target.p1 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                                  kmer.pos = pos.upstream.p1)
    samples_string = paste(upstream_meta.current.target.p1$raw_current, 
                           collapse = ",")
    upstream.string.raw.current.samples.p1 = c(upstream.string.raw.current.samples.p1, 
                                               samples_string)
    upstream.signal.means.p1 = c(upstream.signal.means.p1, upstream_meta.current.target.p1$kmer_mean) #append mean to full list
    upstream.signal.stds.p1 = c(upstream.signal.stds.p1, upstream_meta.current.target.p1$kmer_std) #append std to full list
    
    
    upstream_meta.current.target.p2 = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                                  kmer.pos = pos.upstream.p2)
    samples_string = paste(upstream_meta.current.target.p2$raw_current, 
                           collapse = ",")
    upstream.string.raw.current.samples.p2 = c(upstream.string.raw.current.samples.p2, 
                                               samples_string)
    upstream.signal.means.p2 = c(upstream.signal.means.p2, upstream_meta.current.target.p2$kmer_mean) #append mean to full list
    upstream.signal.stds.p2 = c(upstream.signal.stds.p2, upstream_meta.current.target.p2$kmer_std) #append std to full list
    
    ##############################################################################################
    read.IDs = c(read.IDs, read.name)
    read.indices = c(read.indices, read.indx) #append read index to full index list
    #####################################################################
      
    #Filter out 5mer of interest in the full read, will contain insertions and deletions
    current.read.5mer = seq.df[seq.df$pos %in% target.5mer,]
    #Filter out only the 5mer with base calls and deletions because position for additions could be multiple 
    current.read.5mer = current.read.5mer[current.read.5mer$base != "+", ]#& current.read.5mer$base != "-",] #Deletions can occur instead of U to C
      
    bn2.tmp = current.read.5mer$base[1]
    bn2 = c(bn2, bn2.tmp)
      
    bn1.tmp = current.read.5mer$base[2]
    bn1 = c(bn1, bn1.tmp)
      
    b0.tmp = current.read.5mer$base[3]
    b0 = c(b0, b0.tmp)
      
    bp1.tmp = current.read.5mer$base[4]
    bp1 = c(bp1, bp1.tmp)
      
    bp2.tmp = current.read.5mer$base[5]
    bp2 = c(bp2, bp2.tmp)
      
    qn2.tmp = current.read.5mer$base_quality[1]
    qn2 = c(qn2, qn2.tmp)
      
    qn1.tmp = current.read.5mer$base_quality[2]
    qn1 = c(qn1, qn1.tmp)
      
    q0.tmp = current.read.5mer$base_quality[3]
    q0 = c(q0, q0.tmp)
      
    qp1.tmp = current.read.5mer$base_quality[4]
    qp1 = c(qp1, qp1.tmp)
      
    qp2.tmp = current.read.5mer$base_quality[5]
    qp2 = c(qp2, qp2.tmp)
      
  }
}

##########
#Checking results individually#
head(raw.current.samples)
View(read.indices[1:1640])
View(b0)
View(bn2)
print(length(read.indices))
tail(read.indices)
tail(raw.signals.df)
print(length(read.indices))
print(length(read.indices[1:1640]))
print(length(read.IDs))
print(length(bn2))
print(length(bn1))
print(length(b0))
print(length(bp1))
print(length(bp2))
print(length(qp2))
print(length(signal.means.n2))
print(length(signal.stds))
print(length(raw.current.samples))
print(length(string.raw.current.samples))
head(raw.current.samples, n = 2L)
head(string.raw.current.samples, n = 2L)
######################################

read.indices = read.indices[1:1640] #This was to deal with issue with 100% PSMB2 prep

upstream.signal.means.n2 = c()
upstream.signal.stds.n2 = c()
#upstream.raw.current.samples.n2 = c()
upstream.string.raw.current.samples.n2 = c()

n=3
read.IDs[1:1640]
#"read_index"=read.indices,
tmp.meta.df = data.frame("read_ID"=read.IDs, 
                         "Bn2"=bn2, "Bn1"=bn1, "B0"=b0, "Bp1"=bp1, "Bp2"=bp2, 
                         "Qn2"=qn2, "Qn1"=qn1, "Q0"=q0,"Qp1"=qp1, "Qp2"=qp2, 
                         
                         "current_mean_n2"=signal.means.n2,
                         "current_mean_n1"=signal.means.n1,
                         "current_mean_0"=signal.means.0,
                         "current_mean_p1"=signal.means.p1,
                         "current_mean_p2"=signal.means.p2,
                         
                         "current_std_n2"=signal.stds.n2,
                         "current_std_n1"=signal.stds.n1,
                         "current_std_0"=signal.stds.0,
                         "current_std_p1"=signal.stds.p1,
                         "current_std_p2"=signal.stds.p2,
                         
                         "raw_current_n2"=string.raw.current.samples.n2,
                         "raw_current_n1"=string.raw.current.samples.n1,
                         "raw_current_0"=string.raw.current.samples.0,
                         "raw_current_p1"=string.raw.current.samples.p1,
                         "raw_current_p2"=string.raw.current.samples.p2,
                         
                         "upstream_current_mean_n2"=upstream.signal.means.n2,
                         "upstream_current_mean_n1"=upstream.signal.means.n1,
                         "upstream_current_mean_0"=upstream.signal.means.0,
                         "upstream_current_mean_p1"=upstream.signal.means.p1,
                         "upstream_current_mean_p2"=upstream.signal.means.p2,
                         
                         "upstream_current_std_n2"=upstream.signal.stds.n2,
                         "upstream_current_std_n1"=upstream.signal.stds.n1,
                         "upstream_current_std_0"=upstream.signal.stds.0,
                         "upstream_current_std_p1"=upstream.signal.stds.p1,
                         "upstream_current_std_p2"=upstream.signal.stds.p2,
                         
                         "upstream_raw_current_n2"=upstream.string.raw.current.samples.n2,
                         "upstream_raw_current_n1"=upstream.string.raw.current.samples.n1,
                         "upstream_raw_current_0"=upstream.string.raw.current.samples.0,
                         "upstream_raw_current_p1"=upstream.string.raw.current.samples.p1,
                         "upstream_raw_current_p2"=upstream.string.raw.current.samples.p2
                         
                         
                         )
####Correct csv prep
View(tmp.meta.df)

Filtered.meta.df = na.omit(tmp.meta.df)
View(Filtered.meta.df)

print(nrow(tmp.meta.df))
print(nrow(Filtered.meta.df))

fwrite(tmp.meta.df, "PSMB2_0perc_unmod_pre_finalfilter_experiment2_correct.csv")
fwrite(Filtered.meta.df, "PSMB2_0perc_unmod_post_finalfilter_experiment2_correct.csv")

test.tpm.meta.df = tmp.meta.df
View(test.tpm.meta.df)
fwrite(test.tpm.meta.df, "PSMB2_0perc_unmod_pre_finalfilter_experiment2_TEST_correct.csv", na = 'NA')


#####

test = as.data.frame(cbind(read.indices, read.IDs, bn2, bn1, b0, bp1, bp2,
                           qn2, qn1, q0, qp1, qp2, signal.means, signal.stds,
                           string.raw.current.samples))

View(test)



Filtered.meta.df = na.omit(test)
View(Filtered.meta.df)

#sliced.read.df = current.read.df[current.read.df$pos %in% target.positions,]
#pos.base.df = data.frame("pos"=positions.list,"base" = nucleutide.list, 
#                         "base_quality"=base.qual.list)

#head(seq.df)
#print(preserved.15mer)
#################################

####This code works to combine filtered reads with eventalign results###





# how to convert binary base quality to list of integers
base.quality = bam.df$qual[row]
phred.q = base.quality
phred.q <- PhredQuality(phred.q)
phred.q <- as(phred.q, "IntegerList")  # quality scores
phred.q.lst = phred.q[[1]]

print(phred.q.lst)

###################
processTime<-c(.05,.06)
deadtime<-c(.38,.36)

cycles<-list(list(.08,.10,.07), list(.07,.09,.38))
print(cycles)
cycles<-list(c(.08,.10,.07), c(.07,.09,.38))
print(cycles)
print(cycles[[1]][[2]])

test<-as.data.frame(cbind(processTime,cycles,deadtime))
View(test)
test$cycles[[1]][[2]]

Cycles = c()
print(Cycles)
Cycles = c(Cycles, list(L1))
Cycles = c(Cycles, list(L2))
Cycles = c(Cycles, list(L3))
L1 = c(.08,.10,.07)
L2 = c(.07,.09,.38)
L3 = c(4,4,4,5,5,5,5,6)
Combined = list(L1, L2)
print(Combined)
print(Cycles)

print(summary$V1[which(summary$V2 == "81a08372-83ee-44ca-aa0d-370acf8ec8a9")])

