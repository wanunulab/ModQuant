# ModQuant
This package can be used to quantify site-specific RNA modifications in direct RNA nanopore sequenced (DRS) reads through supervised machine learning (ML) models that are trained with synthetic transcripts that replicate the sequence of the native reads and contain the modified or unmodified base (control) of interest. The labeled dataset for ML training comprises basecalling features and corresponding signal features from the ionic current of the aligned reads. The main input files for this tool are the aligned reads (bam file) and the resquiggled eventalign signal data generated through **Nanopolish**, which can be downloaded from [here](https://github.com/jts/nanopolish). Our tool was developed for training ML models for quantifying the occupancy of pseudouridine sites (ψ) in native human mRNAs. For each site, a model was trained with a corresponding DRS dataset of ψ-modified and unmodifed synthetic mRNA standards. Please read the step by step guide on how to generate and prepare your DRS data for ML training and testing at these specific site.  
## System Requirement

### Hardware Requirements
Package requires only a standard computer with enough RAM (>=16GB recommended) and storage to support the in-memory operations.

### Software Requirements
#### OS Requirements 
This package is supported for Windows and macOS has been tested on the following systems:

  + Windows10 (19044.1706) 
  + mcOS:Catalina(10.15.6)

#### R Dependecies
##### R 
  + R (required | Version>/4.0.1 recommended)
  + Rstudio (recommended)
##### R Packages
  + ggplot2(3.3.5), data.table(1.14.2), optparse(1.7.1), dplyr(1.0.8), BiocManager(1.30.16), Biobase(>=2.50.0), Rsamtools(>=2.6.0), GenomicAlignments(1.26.0) **No action is required on your end. These Packages will be installed automatically on your system the first time you run the code (if not already installed)**

#### Python Dependecies
##### Python
  + Python (required | Version>/3.8 recommended)
  + Anaconda (recommended)
##### Python Packages
  + matplotlib(>=3.2.0), numpy(>=1.18.1), pandas(>=1.0.2), scikit-learn(>=1.0.2), scipy(>=1.3.1), seaborn(>=0.11.2) (**No action is required on your end. These Packages will be installed automatically on your system once you run "pip install -r requirements.txt" and it's recommended that you create a conda environment first).**


## Installation Guide
### Installation with Anaconda
If you do not already have anaconda (latest version works) on your system, download it from [here](https://www.anaconda.com/python-r-distribution?utm_campaign=python&utm_medium=online-advertising&utm_source=google&utm_content=anaconda-download&gclid=Cj0KCQjwqPGUBhDwARIsANNwjV5PuBoSMd9M6wqi0PNKTpaCYmY7G9iIUtehV9XzetwajllP-sFybKcaAtIeEALw_wcB)

*We recommend anaconda because it is a convenient and flexible platform for building python environments, and it provides the Jupyter Lab/Notebook IDE, which is great for building machine learning models with scikit-learn if you wish to try ML algorithms that are different from the 5 models used in this package*

Next, if you do not already have RStudio (>=4.0.1) on your system, download it from [here](https://www.rstudio.com/products/rstudio/download/)

There are two ways to download the package: 

1) Click on the green button on top left of this page that says 'code', then click on 'Download ZIP'. Then unzip the compressed file. Move the package into the directory of your choice.

**OR**

2) Open anaconda prompt and install git with ```conda install git``` and then go to your directory of choice and clone the package with ```git clone https://github.com/wanunulab/ModQuant.git```

Navigate to the directory where the package was downloaded.

```cd <\path\to\ModQuant>```

Next, create an environement (python version 3.8) that has pip as dependency (required for downloading the other python dependencies) ```conda create -n <name_of_env> python=3.8 pip``` (takes around a minute) and then activate the environment with ```conda activate <name_of_env>```. Install all the required dependencies for the package with ```pip install -r requirements.txt``` (takes around a minute). 

Next, you will need to download all the required R dependencies, which can be done in the anaconda command prompt. First, you will need to find the location of Rscript.exe if you did not set it up as an environment variable. From our system, the directory path is "C:\Program Files\R\R-4.0.1\bin\Rscript.exe". Next, run ```"<path_to_Rscript.exe>" FeatureExtract.R -h``` and the required R libraries will be installed. *Note that this might take some time (couple of minutes for our computer) and make sure to have the to your Rscript.exe in quotations*. The command prompt should print the input options. Now your environment is ready!

*Note that you can run these commands from a different directory, just make sure to provide the correct directory path to the modules in the command line*

## Instructions for Basecalling, Sequence Alignment, and Resquiggling with 
### Basecalling ONT Sequencing Data
You will need to an ONT basecaller (Guppy or Dorado) on your system, which can be downloaded [here](https://community.nanoporetech.com/downloads) if you have an ONT account. Run the following command to basecall your DRS data with GPU-enabled Guppy:    

```guppy/bin/guppy_basecaller -i <path\to\fast5\folder> -s <path\to\output\fastq\folder> -c guppy\data\rna_r9.4.1_70bps_hac.cfg -x auto -r --u_substitution```

### Aligning DRS Data with Minimap2
For this work, the sequence alignment tool we used was Minimap2, which can be downloaded [here](https://github.com/lh3/minimap2). The following commands should be run to prepare a sequence alignment file for each synthetic or native mRNA dataset containing the ψ-modified site of interest:

1) Merge fastq sequencing files into one file: 
```cat path\to\fastqs\*.fastq > <filname>.fastq```
2) Align merged fastq file with the corresponding reference file: 
```minimap2 -ax splice -uf -k14 <filename>.fastq > <filename>.sam```
3) Convert .sam to .bam: 
```samtools view -h -Sb <filename>.sam > <filename>.bam```
4) Sort the bam file: 
```samtools sort <filename>.bam -o <filename>_sorted.bam```
5) Index the bam file to make a **.bai** file: 
```samtools index <filename>_sorted.bam``` 

### Aligning Ionic Current Data with Sequence Aligned Data using Nanopolish 
Once you have downloaded the [Nanopolish package](https://github.com/jts/nanopolish), run the following commands to generate your resquiggled DRS data: 

1) Index fast5 files with the merged fastq file: 
```./nanopolish index -d path\to\fast5\folder path\to\merged.fastq``` 
This command may take some time and will output the following files: **.index**, **.index.fai**, **.index.gzi**, and **.index.readdb**.
2) Run the *eventalign* resquiggle tool to create *eventalign* **summary** and *eventalign* **.txt** files:  
```./nanopolish eventalign --read path\to\merged.fastq --bam path\to\gene.bam --genome path\to\reference --scale-events > gene.txt --summary=gene_summary --samples --signal-index``` 


## ModQuant modules
### Feature Extraction
Parses basecalling and signal features from a user-specified location of interest on aligned RNA reads. Required input files are the sorted and indexed .bam file, the nanopolish eventalign .txt file that contains the resquiggled ionic current data from the **Fast5** files, and the nanopolish summary file. The output file is a .csv file with the user-specified features (columns) for every read (rows) that has passed the filtering conditions. To see the description of all the input options for feature extraction, run ```"<path_to_Rscript.exe>" FeatureExtract.R -h``` in the command line. 

To run feature extraction on our synthetic RNA standards, you will need to download the **.bam**, **.bam.bai**, *eventalign* **summary**, and *eventalign* **.txt** files from [here](https://discover.pennsieve.io/datasets/340). All four of these files are required for features to be extracted and compiled for each training set (and native  dataset). Next, move them into a folder that indicates these files carry information about synthetic reads that replicate the mRNA sequence of mRNA transcripts bearing a site-specific ψ modification. To extract the basecalls, quality scores, and k-mer signal information (mean, standard deviation, samples, and dwell-time of the ionic current) surrounding the ψ-modified site **target base position (-t)**, and the **2 neighboring bases in the 5' and 3' direction for a total of 5 (-n)** on synthetic or native reads with a **minimum read length (-l)** and **minimum mapping quality (-m)** of 550nt and 50, respectively, we can run the following in the command line:  

```<path_to_Rscript.exe> FeatureExtract.R -b <path\to\gene_sorted.bam> -s <path\to\gene_summary> -e <path\to\gene.txt> -m 50 -l 550 -t 511 -r 15 -p 10 -n 5 -i yes -o <path\to\gene_100_perc_modified_features.csv>```

### Feature Preparation
Prior to training ML models, the data containing the extracted features from the unmodified and modified reads needs to be labeled, combined, and prepared into one dataframe. The user has the option of adding the Fourier coefficients of the k-mer current signals to the feature space. By default, the 1st quartile (25%), 2nd quartile (50%), 3rd quartile (75%) of each k-mer signal of interest are added to the feature space. We can easily do this with the following command:

```python PrepFeatures.py -cf path\to\gene_0_perc_unmodified_features.csv> -mf <path\to\gene_100_perc_modified_features.csv> -fc 3 -o gene_training_data```
The output will be a .pkl file stored in a folder called *prepared_training_data*. 

### Training, Testing, and Applying ML models for site-specific ψ quantification
Once the feature space is prepared with both modifed and unmodified reads, we can train supervised machine learning models and assess their performance on ψ classification. We can also save these models and apply them on these exact sites found in native human mRNA to quantify their occupancy. Please go to the [tutorial](https://github.com/wanunulab/ModQuant/blob/main/Pseudouridine_ML_Quantification_Tutorial.ipynb) to generate, test, and apply these models.
