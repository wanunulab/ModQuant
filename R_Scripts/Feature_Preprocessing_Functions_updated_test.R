########This function the cigar string of the read to output the########## 
########base position, basecall, and quality score for every nt in a read#
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
      ##################################################################################################################################
      #nucleutide.list = c(nucleutide.list,rep("+",cig.length.tmp)) #Uncomment to get original output
      inserted.nucleotides = seq.temp[c(counter.on.sequence:(counter.on.sequence+cig.length.tmp-1))]
      inserted.nucleotides = gsub("[[:<:]]", "ins", inserted.nucleotides, perl=TRUE)
      nucleutide.list = c(nucleutide.list, inserted.nucleotides)
      ##################################################################################################################################
      #base.qual.list = c(base.qual.list,rep(-2,cig.length.tmp))#switched -1 to -2 for insertion, uncomment to place -2 QS for insertions
      base.qual.list = c(base.qual.list,phred.q.lst[c(counter.on.sequence:(counter.on.sequence+cig.length.tmp-1))])
      
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


############################################
########This function extracts and extends the current mean,######## 
########standard deviation, dwell time, and raw sample points####### 
########for 5mer events from nanopolish for a single read###########
raw_trace_meta_extract_func <- function(read.raw.signal.df, kmer.pos){
  read.raw.signal.df.pos = read.raw.signal.df[which(read.raw.signal.df$position == kmer.pos),]
  if(nrow(read.raw.signal.df.pos) >= 1){
    #print('HIT')
    samples_string = c()
    samples = c()
    dwell_time = 0 #initialize dwell time for 5mer
    
    for (i in c(1:nrow(read.raw.signal.df.pos))){
      # splitting the string of floats
      splitted.string = strsplit(read.raw.signal.df.pos$samples[i],",")[[1]]
      dwell.tmp = as.numeric(read.raw.signal.df.pos$event_length[i])#[[1]] ###
      #print(dwell.tmp) ###
      # convert splitted strings to numeric values
      samples.tmp = as.numeric(splitted.string)
      # append the current row's signals to the main list
      samples = c(samples, samples.tmp)
      samples_string = c(samples_string, read.raw.signal.df.pos$samples[i])
      dwell_time = dwell_time + dwell.tmp ###
      #print(dwell_time)###
    }
    samples_string = paste(samples_string, collapse = ",") #full 5mer current datapoints
    signal.mean = mean(samples) #full 5mer sample mean
    signal.sd = sd(samples) #full 5mer sample std
  } else {
    #print('No raw current')
    samples_string = NA
    signal.mean = NA #full 5mer sample mean
    signal.sd = NA #full 5mer sample std
    dwell_time = NA #full 5mer sample std
  }
  signal.meta.df = data.frame("kmer_mean"=signal.mean,"kmer_std"=signal.sd, "kmer_dt"=dwell_time,"raw_current"=samples_string)
  #############
  return(signal.meta.df)
}
#############################################

##This function filters reads on how many bases they have in the 15mer region##
synthetic_filter_func <- function(current.read.df, target.positions){
  sliced.read.df = current.read.df[current.read.df$pos %in% target.positions,]
  if (nrow(sliced.read.df) == 0){
    counter.bases = 0
  }
  else{
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
        ######################################
        #if (current.base == "A" | current.base == "C" | current.base == "G" | current.base == "T"){
        #  counter.bases = counter.bases + 1
        #}
        ######################################
      }
    }
  }
  return(counter.bases)
}


###Function for raw signal prep from eventalign###
#Check if aligned reads are on plus or minus strand
strand_check_func <- function(whole.bam.df){
  strand.pos = head(whole.bam.df$strand[1])
  return(strand.pos)
}

#This function allows for optional filtering of reads based on mapq and length
read_filter_func <- function(whole.bam.df, min.mapping.qual.score = NULL, min.read.length = NULL ){
  
  if (!is.null(min.read.length)){
    min.read.length = as.integer(min.read.length)
  }
  
  if (!is.null(min.mapping.qual.score)){
    min.mapping.qual.score = as.integer(min.mapping.qual.score)
  }
  
  if (!is.null(min.read.length) && is.null(min.mapping.qual.score)){
    cat("\nOnly provided minimum read length filter")
    if (min.read.length <= 0){
      stop("\nError: Minimum read length should be a value greater than zero!")
    }
    else{
      updated.bam.df = whole.bam.df[nchar(whole.bam.df$seq) >= min.read.length, ]
    }
  }
  
  else if (!is.null(min.mapping.qual.score) && is.null(min.read.length)){
    cat("\nOnly provided minimum mapping quality filter")
    if (min.mapping.qual.score < 0 || min.mapping.qual.score > 60){
      stop("\nError: Minimum mapping quality score needs to be between zero and sixty (you should know this)!")
    }
    else{
      updated.bam.df = whole.bam.df[whole.bam.df$mapq >= min.mapping.qual.score, ]
    }
  }
  
  else if (!is.null(min.mapping.qual.score) && !is.null(min.read.length)){
    cat("\nProvided both a minimum mapping quality and read length filter")
    if (min.read.length <= 0 || min.mapping.qual.score < 0 || min.mapping.qual.score > 60){
      stop("\nError: Mapping quality score needs to be between zero and sixty and read length should be a value greater than zero!")
    }
    else{
      updated.bam.df = whole.bam.df[nchar(whole.bam.df$seq) >= min.read.length & whole.bam.df$mapq >= min.mapping.qual.score, ]
    }
  }
  
  else {
    updated.bam.df = whole.bam.df
    cat("\nProvided neither a minimum mapping quality nor read filter condition")
  }
  
  return(updated.bam.df)
}



#####Function for compiling all features into one dataframe######
###Add a condition called con.insertions###
feature_compile_func <- function(bam.input, summary.file, event.align.file, 
                                 position.of.interest, num.region.of.interest, region.error.limit,
                                 num.bases.and.5mers, consider.insertions){
  
  # strand check 
  rna.strand = strand_check_func(bam.input)
  
  # number of bases in target region to extract must be odd
  if ((num.region.of.interest %% 2) == 0){
    num.region.of.interest = num.region.of.interest + 1
  }
  
  # number of bases and 5mers to extract must be odd
  if ((num.bases.and.5mers %% 2) == 0){
    num.bases.and.5mers = num.bases.and.5mers + 1
  }
  
  # deal with offset between position of base in BAM and 5mer frame in eventalign
  pos.target.5mer.0 = position.of.interest - 3 
  
  region.of.interest = seq(position.of.interest-floor(num.region.of.interest/2), position.of.interest+floor(num.region.of.interest/2), 1)
  target.bases = seq(position.of.interest-floor(num.bases.and.5mers/2), position.of.interest+floor(num.bases.and.5mers/2), 1) #(509, 510 ,511, 512, 513)
  
  
  
  # set up a vector with the positions of 5mer frames to extract from every read
  pos.all.5mers = seq(pos.target.5mer.0-floor(num.bases.and.5mers/2), pos.target.5mer.0+floor(num.bases.and.5mers/2), 1) #(506 ,507, 508, 509, 510)
  
  # denotes variable and list names of 5mer frames in the selected region
  basecalls.and.5mers.positions = seq(0-floor(num.bases.and.5mers/2), 0+floor(num.bases.and.5mers/2), 1)#(-2, -1, 0, +1, +2)
  
  # vector containing frames of 5mers that are of interest and order depends on rna strand
  all.5mer.frames = c()
  idx = 1
  for (i in basecalls.and.5mers.positions){
    if (rna.strand == "+"){
      flipped.pos.all.5mers = rev(pos.all.5mers)
      all.5mer.frames = c(all.5mer.frames, flipped.pos.all.5mers[idx])
    }
    if (rna.strand == "-"){
      ##assign(name.5mer.positions, pos.all.5mers[idx])##removed
      all.5mer.frames = c(all.5mer.frames, pos.all.5mers[idx])
    }
    idx = idx + 1
  }
  
  read.IDs = c()
  Insertions.Freq.A = c()
  Insertions.Freq.C = c()
  Insertions.Freq.G = c()
  Insertions.Freq.T = c()
  Insertions.Total = c()
  for (i in basecalls.and.5mers.positions) {
    
    if (i < 0){
      name.basecalls = paste("B.n", abs(i), sep = "")
      name.basequals = paste("Q.n", abs(i), sep = "")
      
      name.samples = paste("raw.current.samples.n", abs(i), sep = "")
      name.means = paste("signal.mean.n", abs(i), sep = "")
      name.stds = paste("signal.std.n", abs(i), sep = "")
      name.dts = paste("signal.dt.n", abs(i), sep = "")
      
      assign(name.basecalls, c())
      assign(name.basequals, c())
      
      assign(name.samples, c())
      assign(name.means, c())
      assign(name.stds, c())
      assign(name.dts, c())
    }
    
    if (i == 0){
      name.basecalls = paste("B.", i, sep = "")
      name.basequals = paste("Q.", i, sep = "")
      
      name.samples = paste("raw.current.samples.", i, sep = "")
      name.means = paste("signal.mean.", i, sep = "")
      name.stds = paste("signal.std.", i, sep = "")
      name.dts = paste("signal.dt.", i, sep = "")
      
      assign(name.basecalls, c())
      assign(name.basequals, c())
      
      assign(name.samples, c())
      assign(name.means, c())
      assign(name.stds, c())
      assign(name.dts, c())
    }
    
    if (i > 0){
      name.basecalls = paste("B.p", i, sep = "")
      name.basequals = paste("Q.p", i, sep = "")
      
      name.samples = paste("raw.current.samples.p", i, sep = "")
      name.means = paste("signal.mean.p", i, sep = "")
      name.stds = paste("signal.std.p", i, sep = "")
      name.dts = paste("signal.dt.p", i, sep = "")
      
      assign(name.basecalls, c())
      assign(name.basequals, c())
      
      assign(name.samples, c())
      assign(name.means, c())
      assign(name.stds, c())
      assign(name.dts, c())
    }
  }
  
  #target.positions = seq(504, 518, 1)
  #target.5mer = seq(509, 513, 1)
  #fifteenMer.error.limit = 12
  # main loop
  
  ##options(width = 80)
  ##extra <- nchar('||100%')
  ##width <- options()$width
  ##step <- round(j / nrow(bam.input) * (width - extra))
  ##text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
  ##                strrep(' ', width - step - extra), round(j / nrow(bam.input) * 100))
  ##cat(text)
  ##Sys.sleep(0.05)
  ##cat(if (j == nrow(bam.input)) '\n' else '\014')
  
  #pb = txtProgressBar(min = 0, max = nrow(bam.input), style = 3, width = 80, char = "\U03A8")
  pb = txtProgressBar(min = 0, max = nrow(bam.input), style = 3, width = 50, char = "=")
  
  #print("Reads in Bam file: ")
  #print(nrow(bam.input))
  preserved.inc = 0
  for (j in 1:nrow(bam.input)){
    
    close(pb)
    setTxtProgressBar(pb, j)
    ########################################################################################
    row = j
    seq.df = Decipher_cigar_func(cigar = bam.input$cigar[row],start.pos = bam.input$pos[row],
                                 sequence = bam.input$seq[row], base.quality = bam.input$qual[row])
    
    preserved.region = synthetic_filter_func(current.read.df = seq.df, 
                                             target.positions = region.of.interest)
    if (preserved.region >= region.error.limit){
      preserved.inc = preserved.inc + 1
      read.name = bam.input$qname[row]
      
      #####################
      ##cat("\n")
      ##cat(read.name)
      ##cat("\n")
      ##cat(preserved.region)
      ##cat("\n")
      #####################
      
      # find the read index 
      read.indx = summary.file$V1[which(summary.file$V2 == read.name)]
      read.indx = unique(read.indx)
      # all the rows of the nanopolish output that correspond to our targeted read
      target.rows = which(event.align.file$read_index == read.indx)###issue
      # make a new data frame that has only our read of interest
      read.raw.signal.df = event.align.file[target.rows,]
      ###############################################################################################
      read.IDs = c(read.IDs, read.name)
      ###############################################################################################
      #Parse out bases of interest in the full read, including insertions
      current.read.target.bases = seq.df[seq.df$pos %in% target.bases,] #positions 509, 510, 511, 512, 513
      
      ##########Here we set a condition to consider reads without insertions in the target region##########
      if (consider.insertions == "no"){
        #current.read.target.bases = current.read.target.bases[current.read.target.bases$base != "+", ]#Deletions can be included, but not insertions
        current.read.target.bases = current.read.target.bases[!grepl("ins", current.read.target.bases$base), ]
      }
      
        ########################################################################################################################################################
        ##############This algorithm adds the highest frequency insertion basecall and the number of insertions in the window of kmers as features##############
        # if (nrow(highest.frequency.base) > 1){
        #   #check for quality score tie breaker
        #   prev.cum.sum.quality.score = 0
        #   for(j in 1:nrow(highest.frequency.base)){
        #     insertion.base.call = highest.frequency.base$insertion_base[j]
        #     print("current insertion base call: ")
        #     print(insertion.base.call)
        #     cum.sum.quality.score = sum(insertions.current.read.target.bases[insertions.current.read.target.bases$base == insertion.base.call, ]$base_quality)
        #     print("current insertion base call cum sum qs: ")
        #     print(cum.sum.quality.score)
        #     if (cum.sum.quality.score > prev.cum.sum.quality.score){
        #       base.largest.cum.sum.quality.score = insertion.base.call
        #     } else if (cum.sum.quality.score == prev.cum.sum.quality.score){
        #       first.index.current.insertion.base = which(current.read.target.bases$base == insertion.base.call)[1]
        #       first.index.prev.insertion.base = which(current.read.target.bases$base == highest.frequency.base$insertion_base[j-1])[1]
        #       print("current index")
        #       print(first.index.current.insertion.base)
        #       print("previous index")
        #       print(first.index.prev.insertion.base)
        #       
        #       if(first.index.current.insertion.base < first.index.prev.insertion.base){
        #         base.largest.cum.sum.quality.score = insertion.base.call
        #       }
        #       else{
        #         base.largest.cum.sum.quality.score = highest.frequency.base$insertion_base[j-1]
        #       }
        #     }
        #     prev.cum.sum.quality.score = cum.sum.quality.score
        #   }
        #   highest.frequency.base.insertion = base.largest.cum.sum.quality.score
        # }
        ########################################################################################################################################################
      #else{
      #  highest.frequency.base.insertion = 0
      #  max.number.insertions = 0
      #}
      
      #print("...")
      #print(current.read.target.bases)
      
      #print(insertions.n)
      #print(insertion.base.frequencies)
      
      #print(highest.frequency.base)
      #print(nrow(highest.frequency.base))
      #print("Results: ")
      #print(highest.frequency.base.insertion)
      #print(max.number.insertions)

      #print(which.max(insertion.base.frequencies))
      #print(current.read.target.bases[!grepl("ins", current.read.target.bases$base), ])
      #print(current.read.target.bases[grepl("ins", current.read.target.bases$base), ])
      #print(preserved.region)
      ######################################################
      ###########Take reverse if strand is "-"##############
      ##############Then take the complement################
      if (rna.strand == "-"){
        current.read.target.bases = apply(current.read.target.bases, 2, rev) #reverse the dataframe containing positions, basecalls, and quality scores
        current.read.target.bases = as.data.frame(current.read.target.bases)
        #current.read.target.bases$base <- gsub('A', 'T',
        #           gsub('T', 'A',
        #           gsub('G', 'C', 
        #           gsub('C', 'G', current.read.target.bases$base)))) #take the complement of the basecalls in the dataframe
        
        current.read.target.bases$base = chartr("ACGT", "TGCA", current.read.target.bases$base)
      }
      
      
      insertions.n = nrow(current.read.target.bases[grepl("ins", current.read.target.bases$base), ])#sample size of insertions in the window of interest
      #insertion.A.frequency = 0
      #insertion.C.frequency = 0
      #insertion.G.frequency = 0
      #insertion.T.frequency = 0
      #insertion.base.frequencies = 0
      #max.number.insertions = 0
      #highest.frequency.base = 0
      #highest.frequency.base.insertion = 0
      insertions.current.read.target.bases = 0
      
      if (insertions.n > 0){
        insertions.current.read.target.bases = current.read.target.bases[grepl("ins", current.read.target.bases$base), ]#dataframe only containing insertions
        insertion.base.frequencies = table(insertions.current.read.target.bases$base) #Frequency of insertions for each base 
        insertion.base.frequencies = as.data.frame(insertion.base.frequencies) #convert table to dataframe
        names(insertion.base.frequencies)[names(insertion.base.frequencies) == 'Var1'] <- 'insertion_base' #rename column
        
        insertion.A.frequency = insertion.base.frequencies[insertion.base.frequencies$insertion_base == "insA", ]$Freq#/insertions.n
        if (length(insertion.A.frequency) == 0){
          insertion.A.frequency = 0
        }
        else{
          insertion.A.frequency = round(insertion.A.frequency/insertions.n, digits=3)
        }

        
        insertion.C.frequency = insertion.base.frequencies[insertion.base.frequencies$insertion_base == "insC", ]$Freq#/insertions.n
        if (length(insertion.C.frequency) == 0){
          insertion.C.frequency = 0
        }
        else{
          insertion.C.frequency = round(insertion.C.frequency/insertions.n, digits=3)
        }
        
        insertion.G.frequency = insertion.base.frequencies[insertion.base.frequencies$insertion_base == "insG", ]$Freq#/insertions.n
        if (length(insertion.G.frequency) == 0){
          insertion.G.frequency = 0
        }
        else{
          insertion.G.frequency = round(insertion.G.frequency/insertions.n, digits=3)
        }
        
        insertion.T.frequency = insertion.base.frequencies[insertion.base.frequencies$insertion_base == "insT", ]$Freq#/insertions.n
        if (length(insertion.T.frequency) == 0){
          insertion.T.frequency = 0
        }
        else{
          insertion.T.frequency = round(insertion.T.frequency/insertions.n, digits=3)
        }
      } 
      
      else{
        insertions.n = 0
        insertion.A.frequency = 0
        insertion.C.frequency = 0
        insertion.G.frequency = 0
        insertion.T.frequency = 0
      }
      
      ##during alignment the insertions on the antisense strand undergo reverse complement##
      if (rna.strand == "-"){
        insertion.C.frequency.comp = insertion.G.frequency
        insertion.G.frequency.comp = insertion.C.frequency
        insertion.A.frequency.comp = insertion.T.frequency
        insertion.T.frequency.comp = insertion.A.frequency
        
        insertion.C.frequency = insertion.C.frequency.comp
        insertion.G.frequency = insertion.G.frequency.comp
        insertion.A.frequency = insertion.A.frequency.comp
        insertion.T.frequency = insertion.T.frequency.comp
        
        
      } 
        
        #max.number.insertions = max(insertion.base.frequencies)
        #insertion.base.frequencies = as.data.frame(insertion.base.frequencies)
        #names(insertion.base.frequencies)[names(insertion.base.frequencies) == 'Var1'] <- 'insertion_base'
        #highest.frequency.base = insertion.base.frequencies[insertion.base.frequencies$Freq == max.number.insertions, ]
        
        #highest.frequency.base.insertion = highest.frequency.base$insertion_base#insertion.base.frequencies$insertion_base[insertion.base.frequencies$Freq == max.number.insertions, ]
        
      #######################################################
      #print("Total number of insertions: ")
      #print(insertions.current.read.target.bases)
      
      #print("Number of insertion calls for A: ") 
      #print(insertion.A.frequency)
      
      #print("Number of insertion calls for C: ") 
      #print(insertion.C.frequency)
      
      #print("Number of insertion calls for G: ") 
      #print(insertion.G.frequency)
      
      #print("Number of insertion calls for T: ") 
      #print(insertion.T.frequency)
      
      #print("Total Number of Insertions in read: ") 
      #print(insertions.n)
      
      for (i in 1:length(all.5mer.frames)){ #1, 2, 3, 4, 5
        current.frame.pos = all.5mer.frames[i] #(510, 509, 508, 507, 506)
        current.local.pos = basecalls.and.5mers.positions[i] #(-2, -1, 0, +1, +2)
        
        ######################################################################################################
        #####Here we need to take into account the new insertion features#####
        current.read.target.bases = current.read.target.bases[!grepl("ins", current.read.target.bases$base), ]
        bc.tmp = current.read.target.bases$base[i] #basecall at a particular target position 
        qs.tmp = current.read.target.bases$base_quality[i] #quality score at a particular target position
        ######################################################################################################        
        meta.current.target = raw_trace_meta_extract_func(read.raw.signal.df = read.raw.signal.df, 
                                                          kmer.pos = current.frame.pos)
        
        samples_string = paste(meta.current.target$raw_current, collapse = ",")
        
        if (current.local.pos < 0){
          
          assign(paste("B.n", abs(current.local.pos), sep=""), 
                 c(get(paste("B.n", abs(current.local.pos), sep="")), bc.tmp))
          
          assign(paste("Q.n", abs(current.local.pos), sep=""), 
                 c(get(paste("Q.n", abs(current.local.pos), sep="")), qs.tmp))
          
          assign(paste("raw.current.samples.n", abs(current.local.pos), sep=""), 
                 c(get(paste("raw.current.samples.n", abs(current.local.pos), sep="")), samples_string ))
          
          assign(paste("signal.mean.n", abs(current.local.pos), sep=""), 
                 c(get(paste("signal.mean.n", abs(current.local.pos), sep="")), meta.current.target$kmer_mean))
          
          assign(paste("signal.std.n", abs(current.local.pos), sep=""), 
                 c(get(paste("signal.std.n", abs(current.local.pos), sep="")), meta.current.target$kmer_std))
          
          assign(paste("signal.dt.n", abs(current.local.pos), sep=""), 
                 c(get(paste("signal.dt.n", abs(current.local.pos), sep="")), meta.current.target$kmer_dt))
          
        }
        
        if (current.local.pos == 0){
          
          assign(paste("B.", current.local.pos, sep=""), 
                 c(get(paste("B.", current.local.pos, sep="")), bc.tmp))
          
          assign(paste("Q.", current.local.pos, sep=""), 
                 c(get(paste("Q.", current.local.pos, sep="")), qs.tmp))
          
          assign(paste("raw.current.samples.", current.local.pos, sep=""), 
                 c(get(paste("raw.current.samples.", current.local.pos, sep="")), samples_string ))
          
          assign(paste("signal.mean.", current.local.pos, sep=""), 
                 c(get(paste("signal.mean.", current.local.pos, sep="")), meta.current.target$kmer_mean))
          
          assign(paste("signal.std.", current.local.pos, sep=""), 
                 c(get(paste("signal.std.", current.local.pos, sep="")), meta.current.target$kmer_std))
          
          assign(paste("signal.dt.", current.local.pos, sep=""), 
                 c(get(paste("signal.dt.", current.local.pos, sep="")), meta.current.target$kmer_dt))
          
        }
        
        if (current.local.pos > 0){
          
          assign(paste("B.p", current.local.pos, sep=""), 
                 c(get(paste("B.p", current.local.pos, sep="")), bc.tmp))
          
          assign(paste("Q.p", current.local.pos, sep=""), 
                 c(get(paste("Q.p", current.local.pos, sep="")), qs.tmp))
          
          assign(paste("raw.current.samples.p", current.local.pos, sep=""), 
                 c(get(paste("raw.current.samples.p", current.local.pos, sep="")), samples_string ))
          
          assign(paste("signal.mean.p", current.local.pos, sep=""), 
                 c(get(paste("signal.mean.p", current.local.pos, sep="")), meta.current.target$kmer_mean))
          
          assign(paste("signal.std.p", current.local.pos, sep=""), 
                 c(get(paste("signal.std.p", current.local.pos, sep="")), meta.current.target$kmer_std))
          
          assign(paste("signal.dt.p", current.local.pos, sep=""), 
                 c(get(paste("signal.dt.p", current.local.pos, sep="")), meta.current.target$kmer_dt))
          
        }
        
      }
      #read.IDs = c(read.IDs, read.name)
      Insertions.Freq.A = c(Insertions.Freq.A, insertion.A.frequency)
      Insertions.Freq.C = c(Insertions.Freq.C, insertion.C.frequency)
      Insertions.Freq.G = c(Insertions.Freq.G, insertion.G.frequency)
      Insertions.Freq.T = c(Insertions.Freq.T, insertion.T.frequency)
      Insertions.Total = c(Insertions.Total, insertions.n)
      
    }
    
  }
  
  #read.IDs = c(1, 2, 3, 4)#
  read.ids.df = data.frame("read_ID"=read.IDs)
  basecalls.df = data.frame("read_ID"=read.IDs)
  quality.scores.df = data.frame("read_ID"=read.IDs)
  
  signal.current.samples.df = data.frame("read_ID"=read.IDs)
  signal.means.df = data.frame("read_ID"=read.IDs)
  signal.stds.df = data.frame("read_ID"=read.IDs)
  signal.dts.df = data.frame("read_ID"=read.IDs)
  
  insertions.a.df = data.frame("Insertions_Freq_A"=Insertions.Freq.A)
  insertions.c.df = data.frame("Insertions_Freq_C"=Insertions.Freq.C)
  insertions.g.df = data.frame("Insertions_Freq_G"=Insertions.Freq.G)
  insertions.t.df = data.frame("Insertions_Freq_T"=Insertions.Freq.T)
  total.insertions.df = data.frame("Total_Insertions"=Insertions.Total)
  
  for (i in 1:length(all.5mer.frames)){  #(1, 2, 3, 4, 5)
    current.local.pos = basecalls.and.5mers.positions[i]
    
    if (current.local.pos < 0){
      
      col_name = paste("B.n", abs(current.local.pos), sep="")
      basecalls.df[col_name] = get(paste("B.n", abs(current.local.pos), sep=""))
      
      col_name = paste("Q.n", abs(current.local.pos), sep="")
      quality.scores.df[col_name] = get(paste("Q.n", abs(current.local.pos), sep=""))
      
      col_name = paste("raw.current.samples.n", abs(current.local.pos), sep="")
      signal.current.samples.df[col_name] = get(paste("raw.current.samples.n", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.mean.n", abs(current.local.pos), sep="")
      signal.means.df[col_name] = get(paste("signal.mean.n", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.std.n", abs(current.local.pos), sep="")
      signal.stds.df[col_name] = get(paste("signal.std.n", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.dt.n", abs(current.local.pos), sep="")
      signal.dts.df[col_name] = get(paste("signal.dt.n", abs(current.local.pos), sep=""))
      
    }
    
    if (current.local.pos == 0){
      
      col_name = paste("B.", abs(current.local.pos), sep="")
      basecalls.df[col_name] = get(paste("B.", abs(current.local.pos), sep=""))
      
      col_name = paste("Q.", abs(current.local.pos), sep="")
      quality.scores.df[col_name] = get(paste("Q.", abs(current.local.pos), sep=""))
      
      col_name = paste("raw.current.samples.", abs(current.local.pos), sep="")
      signal.current.samples.df[col_name] = get(paste("raw.current.samples.", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.mean.", abs(current.local.pos), sep="")
      signal.means.df[col_name] = get(paste("signal.mean.", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.std.", abs(current.local.pos), sep="")
      signal.stds.df[col_name] = get(paste("signal.std.", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.dt.", abs(current.local.pos), sep="")
      signal.dts.df[col_name] = get(paste("signal.dt.", abs(current.local.pos), sep=""))
      
    }
    
    if (current.local.pos > 0){
      
      col_name = paste("B.p", abs(current.local.pos), sep="")
      basecalls.df[col_name] = get(paste("B.p", abs(current.local.pos), sep=""))
      
      col_name = paste("Q.p", abs(current.local.pos), sep="")
      quality.scores.df[col_name] = get(paste("Q.p", abs(current.local.pos), sep=""))
      
      col_name = paste("raw.current.samples.p", abs(current.local.pos), sep="")
      signal.current.samples.df[col_name] = get(paste("raw.current.samples.p", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.mean.p", abs(current.local.pos), sep="")
      signal.means.df[col_name] = get(paste("signal.mean.p", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.std.p", abs(current.local.pos), sep="")
      signal.stds.df[col_name] = get(paste("signal.std.p", abs(current.local.pos), sep=""))
      
      col_name = paste("signal.dt.p", abs(current.local.pos), sep="")
      signal.dts.df[col_name] = get(paste("signal.dt.p", abs(current.local.pos), sep=""))
      
      
    }
    
  }
  # Drop read ID columns from all the feature columns
  basecalls.df =  basecalls.df[ , !(names(basecalls.df) %in% "read_ID")]
  quality.scores.df =  quality.scores.df[ , !(names(quality.scores.df) %in% "read_ID")]
  signal.means.df = signal.means.df[ , !(names(signal.means.df) %in% "read_ID")]
  signal.stds.df = signal.stds.df[ , !(names(signal.stds.df) %in% "read_ID")]
  signal.dts.df = signal.dts.df[ , !(names(signal.dts.df) %in% "read_ID")]
  signal.current.samples.df = signal.current.samples.df[ , !(names(signal.current.samples.df) %in% "read_ID")]
  
  #print("Number of times preserved condition was met")
  #print(preserved.inc)
  #print(nrow(read.ids.df))
  #print(nrow(insertions.a.df))
  #print(nrow(insertions.c.df))
  #print(nrow(insertions.g.df))
  #print(nrow(insertions.t.df))
  #print(nrow(total.insertions.df))
  #print(nrow(basecalls.df))
  
  tmp.meta.df = bind_cols(read.ids.df, basecalls.df, quality.scores.df, 
                          insertions.a.df, insertions.c.df, insertions.g.df,
                          insertions.t.df, total.insertions.df,
                          signal.means.df, signal.stds.df, signal.dts.df,
                          signal.current.samples.df)
  
  
  return(tmp.meta.df)
  
}

















