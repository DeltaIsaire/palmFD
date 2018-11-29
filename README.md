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
plyr
magrittr
BHPMF
```

## Running the code

The main file (00_run_all.R) automatically executes all scripts in sensible order,
thus providing a one-click way to re-generate all output. To run it, launch R and set the working directory to the root folder of the project. For example:
```
setwd("/home/delta/R_Projects/palm_FD")
```

Then you can run the script to execute all the code at once:
```
source(file="00_run_all.R")
```

