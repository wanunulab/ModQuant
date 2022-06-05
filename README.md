# ModQuant
This package can be used to quantify site-specific RNA modifications in nanopore sequenced reads through supervised machine learning (ML) models that are trained with synthetic transcripts that replicate the sequence of the native reads and contain the modified or unmodified base (control) of interest. The labeled dataset for ML training comprises basecalling features and corresponding signal features from the ionic current of the aligned reads. The main input files for this tool are the aligned reads (bam file) and the resquiggled eventalign signal data generated through nanopolish (eventalign .txt and summary files). Our tool was developed for training ML models for quantifying the occupancy of ψ putative sites in native HeLa mRNA. For each site, a model was trained with a corresponding synthetic modified and unmodifed dataset. Please read the step by step guide on how to generate and prepare your data for ML training and testing.  
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
  + matplotlib(>=3.2.0), numpy(>=1.18.1), pandas(>=1.0.2), scikit-learn(>=1.0.2), scipy(>=1.3.1), seaborn(>=0.11.2) (**No action is required on your end. These Packages will be installed automatically on your system once you run "pip install -r requirements.txt" (recommended that you create a conda environment first).**


## Installation Guide
### Windows installation with Anaconda
If you do not already have anaconda (latest version works) on your system, download it from [here](https://www.anaconda.com/python-r-distribution?utm_campaign=python&utm_medium=online-advertising&utm_source=google&utm_content=anaconda-download&gclid=Cj0KCQjwqPGUBhDwARIsANNwjV5PuBoSMd9M6wqi0PNKTpaCYmY7G9iIUtehV9XzetwajllP-sFybKcaAtIeEALw_wcB)

*We recommend anaconda because it is a convenient and flexible platform for building python environments, and it provides the Jupyter Lab/Notebook IDE, which is great for building machine learning models with scikit-learn if you wish to try ML algorithms that are different from the 5 models used in this package*

Next, if you do not already have RStudio (>=4.0.1) on your system, download it from [here](https://www.rstudio.com/products/rstudio/download/)

There are two ways to download the package: 

#. Click on the green button on top left of this page that says 'code', then click on 'Download ZIP'. Then unzip the compressed file. Move the package into the directory of your choice.
#. Open anaconda prompt and install git with ```conda install git``` and then go to your directory of choice and clone the package with ```git clone https://github.com/wanunulab/ModQuant.git```

Navigate to the directory where the package was downloaded.

```cd <\path\to\ModQuant>```

Next, create an environement (python version 3.8) that has pip as dependency (required for downloading the other python dependencies) ```conda create -n <name_of_env> python=3.8 pip``` (takes around a minute) and then activate the environment with ```conda activate <name_of_env>```. Install all the required dependencies for the package with ```pip install -r requirements.txt``` (takes around a minute). You should now be able to print out all the versions of scikit-learn, numpy, and pandas downloaded along with the input options with ```python modML.py -h```  

Next, you will need to download all the required R dependencies, which can be done in the anaconda command prompt. First, you will need to find the location of Rscript.exe if you did not set it up as an environment variable. From our system, the directory path is "C:\Program Files\R\R-4.0.1\bin\Rscript.exe". Next, run ```"<path_to_Rscript.exe>" FeatureExtract.R -h``` and the required R libraries will be installed. *Note that this might take some time (couple of minutes for our computer) and make sure to have the to your Rscript.exe in quotations*. The command prompt should print the input options. Now your environment is ready!

*Note that you can run these commands from a different directory, just make sure to provide the correct directory path to modML.py and FeatureExtract.R in the command line*



