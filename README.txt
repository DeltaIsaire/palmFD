========================================================
R code for "palm functional diversity" research  project
MSc Biological sciences, track Ecology & Evolution
University of Amsterdam
========================================================

File paths in all scripts are relative to the /palm_FD directory. Code depends
on the following subdirectories being present:
data, functions, output, graphs
Only the functions directory is available from the repository. The other directories
contain datasets I shouldn't publish, else my supervisors would kill me :/

The main file (00_main.R) automatically executes all scripts in sensible order,
thus providing a one-click way to re-generate all output, e.g.
> setwd("/home/delta/R_Projects/palm_FD")
> source(file="00_main.R")

List of required R packages:
plyr
ape

