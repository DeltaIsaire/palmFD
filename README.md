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

The code of this project is written for R statistical software version 3.6.0 (R Core Team, 2019).
See https://www.R-project.org/ for more information and installation instructions. 

In addition, the following R packages (and their prerequisites) are required:
```
car
corrplot
cowplot
FD
ggplot2
grid
gridExtra
leaps
magrittr
MuMIn
ncf
parallel
plyr
reshape2
scales
sf
spatialreg
spdep
viridis
```

### Compatibility with Windows

The code has been developed and tested on Linux Mint 18 MATE (64 bit). I have tried to make the code fully compatible with Windows, but compatibility is not guaranteed.

## Running the code

The complete analysis is performed in 10 scripts, each of which depends on the output of one or more previous scripts. The filenames have number prefixes to indicate the proper order.

There is a parent script (00_run_all.R) that will execute all scripts in their proper order. Running this parent script is the easiest way to re-generate all output at once.

Each script contains a description of what it does, at the top of the file.

To run any of the scripts, launch R and set the working directory to the root folder of the project. For example:
```
setwd("/home/delta/R_Projects/palm_FD")
```

Then you can run any of the scripts with source(). Keep in mind that each script requires the output of the previous script in the sequence, so for a first-time run you should probably execute the parent script.
```
source(file = "00_run_all.R")
```

Note that the total runtime of the scripts is considerable, in particular the null model procedures. Calculation of multi-dimensional functional diversity indices is computationally expensive, and the null models require it to be done thousands of times. It can take a few weeks on a dedicated server.
