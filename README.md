# Palm Functional Diversity Research Project

R code for data analysis in my MSc research project on functional diversity in palms. 
MSc Biological sciences, track Ecology & Evolution
University of Amsterdam

## Getting Started: cloning the repository

First install Git on your system if you haven't already. For detailed instructions see chapter 1.5 of the Git book: https://git-scm.com/book/en/v2

Second, create a directory for the project on your file system.
Then you can clone the repository:
```
cd <directory>
git clone https://github.com/DeltaIsaire/palmFD.git
```

Third, you have to manually add the 'data/' directory, including the datafiles, into the root directory of the project (i.e. the same place this readme is located). Without the data, the scripts have nothing to work with. The data files are not part of the repository because they are not publicly available at this time.

## Prerequisites

The code of this project is written for R statistical software version 3.4.4 (R Core Team, 2018).
See https://www.R-project.org/ for more information and installation instructions. 

In addition, the following R packages are required:
```
parallel - should come pre-installed with your R installation
plyr
magrittr
BHPMF - install with devtools::install_github("fisw10/BHPMF")
reshape2
FD
```

### Compatibility with Windows

The code has been developed and tested on Linux Mint 18. I have tried to make the code fully compatible with Windows, but compatibility is not guaranteed.

In particular, the null models for functional diversity are generated with parallel processing using 'parallel::mclapply'. This R function uses forking, which is not supported on Windows.


## Running the code

There are two series of scripts.
The first is the test series, whose filenames begin with 'test\_'. This script series is used for code development. There is a parent script ('00\_run\_all\_tests\_.R') that will automatically execute all test scripts.
The second series of scripts are the main scripts, which are like the no-nonsense final version of the code. There is also a parent script for these ('00\_run\_all\_main.R').

Each script contains a description of what it does, at the top of the file.

To run any of the scripts, launch R and set the working directory to the root folder of the project. For example:
```
setwd("/home/delta/R_Projects/palm_FD")
```

Then you can run any of the scripts with source(). Keep in mind that each script requires the output of the previous script in the sequence, so for a first-time run you should execute the parent scripts.
```
source(file = "00_run_all.R")
source(file = "00_run_all_tests.R")
```

