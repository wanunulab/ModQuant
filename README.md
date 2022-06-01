# ModQuant
This package can be used to quantify site-specific RNA modifications in nanopore sequenced reads through supervised machine learning (ML) models that are trained with synthetic transcripts that replicate the sequence of the native reads and contain the modified or unmodified base (control) of interest. The labeled dataset for ML training comprises basecalling features and corresponding signal features from the ionic current of the aligned reads. The main input files for this tool are the aligned reads (bam file) and the resquiggled eventalign signal data generated through nanopolish (eventalign .txt and summary files). Our tool was developed for training ML models for quantifying the occupancy of ψ putative sites in native HeLa mRNA. For each site, a model was trained with a corresponding synthetic modified and unmodifed dataset. Please read the step by step guide on how to generate and prepare your data for ML training and testing.  
## System Requirement

### Hardware Requirements
Package requires only a standard computer with enough RAM (>=16GB recommended) and storage to support the in-memory operations.

### Software Requirements
#### OS Requirements 
This package is supported for Windows and macOS has been tested on the following systems:

Windows10:19044.1706 
mcOS:Catalina(10.15.6)

#### R Dependecies
##### R 
  + R (required | Version>/4.0.1 recommended)
  + Rstudio (recommended)
##### R Packages
  + ggplot2(3.3.5), data.table(1.14.2), optparse(1.7.1), dplyr(1.0.8), BiocManager(1.30.16), Biobase(>=2.50.0), Rsamtools(>=2.6.0), GenomicAlignments(1.26.0) (**No action is required on your end. These Packages will be installed automatically on your system the first time you run the code (if not already installed).**)
#### Python Dependecies
##### Python
  + Python (required | Version>/3.8 recommended
  + Anaconda (recommended)
##### Python Packages
  + matplotlib(>=3.2.0), numpy(>=1.18.1), pandas(>=1.0.2), scikit-learn(>=1.0.2), scipy(>=1.3.1), seaborn(>=0.11.2) (**No action is required on your end. These Packages will be installed automatically on your system once you run "pip install -r requirements.txt" (recommended that you create a conda environment first).**)


## Installation Guide
