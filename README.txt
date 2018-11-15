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

Third, you have to manually add the following subdirectories in the root directory of the project (i.e. the same place this readme is located):
```
data/
output/
graphs/
```
The scripts will only work if these directories are present.
These folders (and their contents) are not part of the repository because they contain the data (raw or processed) used for this project, which are not publicly available at this time.

Finally, you have to make sure the data/ directory contains the required datasets.

## Prerequisites

The code of this project is written for R statistical software (R Core Team, 2018).
See https://www.R-project.org/ for more information and installation instructions. This project was written for R version 3.4.4.

In addition, the following R packages are required:
```
plyr
ape
```

## Running the code

The main file (00_main.R) automatically executes all scripts in sensible order,
thus providing a one-click way to re-generate all output. To run it, launch R and set the working directory to the root folder of the project. For example:
```
setwd("/home/delta/R_Projects/palm_FD")
```

Then you can run the script 00_main.R to execute all the code at once:
```
source(file="00_main.R")
```

