# Author's Log
This is my daily log of activities and notes related to this project.


## Monday november 12, 2018

Start of Author's Log.
Jun suggested to keep a daily log of all my activities on the palm FD project, as a way to keep track of all the files and code I will be producing. That's a great idea, so let's do this.

Met up with Jun this morning to discuss method details. We talked mostly about gap-filling the palm trait data. Goal for this week should be getting trait-filling code to work, at least on a basic level.
Monday-morning meetup with Jun will likely become a weekly thing (Mondays 10am).

Re-organized my code and data files. 
All input and output data files are kept in my R projects folder, ~/R_Projects/Palm_FD.
All code is kept inside one folder within the palm FD folder:~/R_Projects/Palm_FD/code.
Your standard R working directory should therefore be
setwd("/home/delta/R_Projects/Palm_FD")

I have not set up github yet, but I can see that becoming a thing.

Structure of code filing:
For every project step, such as trait-filling, I will create a parent codefile in /R_code (e.g. trait_filling.R). This file contains the executive code for this step.
Any functions created for the process will have their code in a seperate file under R_code/functions. These functions are called in the parent codefile, which should also load them with the source() function.
see also https://nicercode.github.io/guides/functions/

So today I'm learning about functions.
In fact, I'm going to take a few hours today reading some of the guides on https://nicercode.github.io/, because it's a great introduction to getting your hands dirty in R, and I need that introduction.
And also https://arxiv.org/pdf/1210.0530.pdf (downloaded, see code/supporting_documents).
I also vow to abide to Google's R Style Guide: https://google.github.io/styleguide/Rguide.xml (downloaded, see code/supporting_documents).

I should totally set up git / github, for version control and code communication.


## Tuesday November 13, 2018
This morning I studied how Git works, and set up a repository for the project:
https://github.com/DeltaIsaire/palmFD
Now we can code properly and share the work.

Next I begin coding (a.k.a. puzzling) for trait-filling. Git will keep track of what I am doing.
I begin by cross-referencing the species in the distribution, traits and tree datafiles to see if they agree on which palm species exist, and what the discrepancies are. A good exercise in using my new coding methods...


## Wednesday November 14, 2018
A day spent coding.
I now have three complete functions: isOneDimensional, crossReference and multiCrossRef.

Subsequently I used the cross-referencing functions to compare the species listed in the palm distribution, palm traits and palm phylogenetic tree, and produced a list of the species names present in all three datasets.

Starting tomorrow I should subset the palm distribution and palm traits data to the agreed species list, and see if I can calculate genus-level means for palm traits. Removing species from the phylogenetic tree will be more difficult...


## Monday November 19, 2018
Met Jun again this morning. We discussed how to do BHPMF, in particular the use of a hierarchical matrix rather than a phylogenetic tree. 

Last week I finished the code for gap-filling with genus means.
Today I will try to get BHPMF working properly. 

Planning update:
The official starting data is 22 October 2018. I actually started working on the project one week earlier, so the unofficial starting date is 15 October 2018.
That makes today the beginning of week 6 (or officially, week 5). 
According to the planning we made, the proposal should be finished by week 6, and week 7-10 are reserved for trait filling. I spent last week on trait filling already, so there is a bit of overlap.
In conclusion, there's 2 weeks left for the proposal (re-structuring the introduction and finishing the method section), and 3 weeks for trait-filling.
Also, in week 7 or 8 there should be an interim assessment. You should plan a date for that this week. 


## Tuesday November 20, 2018
This morning I read Noble 2009: A quick guide to organizing computational biology projects. There are some useful tips in there, some of them new to me. Definitely recommended material for students like me, together with Wilson et al 2013: best practices for scientific computing. Both PDFs are stored in palm_FD/reference_literature; together with Google's R style guide. While I'm at it, I should also mention the Nice R Code blog, which is a great reference for beginning R coders: https://nicercode.github.io/

Subsequently I reorganized some files (see git log). I also split some of the non-essential testing code from the trait filling script into a seperate file (again see git log). When developing code, you will often be running all kinds of checks, and it is good to store those in a seperate file. That way the main script remains clean, but you still have some reference material available for bug-testing and the like. The general pattern is to name these reference files scriptname_test.R. So essentially, when making a script you work on two files concurrently: the script and the script_test.

The goal for today is comparing the gap-filling results obtained with genus-means and BHPMF. This involves making plots, for which I will be writing some custom plotting functions. These functions combine my best knowledge of plotting with advice from the Nice R Code blog, and will help with plotting now and in the future. First custom function is Scatterplot, cause I will be making scatterplots today.

What I ended up doing:
- wrote functions Scatterplot and GraphPDF, see palm_FD/functions/plotting.functions.R
- Made the first Scatterplot, comparing Genus mean and BHPMF estimates of stem height.

A lot of time today was spent working out making quality plotting functions. That time investment will pay off later. 


## Wednesday November 20, 2018
Today I should finish comparing the results of genus-mean gap-filling with BHPMF gap-filling.

I added the MultiScatter function, which provides a quick way to create decent quality scatterplots. Subsequently I used the function for comparing genus-mean trait predictions with BHPMF trait predictions, for our three traits stem height, blade length and fruit length.

The next task is comparing predicted trait values with boxplots. hile preparing data for making these boxplots, I found out that the BHPMF results aren't quite what they should be. The BHPMF output appears to also estimate trait values that were already known. Case in point is the palm species Eremospatha_tessmanniana.

The original data:
```
> traits[which(traits$species == "Eremospatha_tessmanniana"), ]
                      species       genus stem.height blade.length fruit.length
1413 Eremospatha_tessmanniana Eremospatha         150          0.8           NA
```
Data filled with BHPMF:
```
> traits.BHPMF[which(traits.BHPMF$species == "Eremospatha_tessmanniana"), ]
                      species       genus stem.height blade.length fruit.length
1209 Eremospatha_tessmanniana Eremospatha     47.3923        0.803       2.4764
```

In hindsight, I could have known. The description of the BHPMF::GapFilling function says "Will save the mean and standard deviations (uncertainties) of the predictions for both missing and observed values into the given file paths."

That's all very well, but it might not be what we need. We should probably assume that the observed values are 'correct' and do not need to be predicted. Then instead of taking the BHPMF output as-is, we should extract only the predicted values for trait values which were missing. In theory that's easy enough, relatively speaking.

I have implemented this. In practice it wasn't entirely straightforward. It involved updating my custom GapFill function to handle the case where the 'by' column in the two provided dataframes is a factor, with differing levels between the datasets. And also, bugfix the gapfilling routine: I replaced an unlist() with as.numeric(). And then of course I had to check if everything now works as intended, which thankfully it does. 


## Friday November 23, 2018
Yesterday I emailed my daily supervisor with some questions about BHPMF and the inputs/outputs used in the process. Now it's time to get our hands dirty and see if we can make this work. 

Yesterday I already started playing around with the inputs of the BHPMF GapFilling function. There are a few issues to work out. Specifically:
- BHPMF got stuck in an infinite loop if some observed trait values are exactly 0
	Tested solution: set these values to a small non-zero value.
- BHPMF doesn't work if a species has NA for all included traits, even if only one column of trait values (i.e. only a single trait) is provided as input.
	Tested solutions: remove these species, include a dummy trait with all values 1.
All the test code is in '1_trait_filling_test.R'

Today I continue testing, and the goal is to report outcomes as well as comparisons of outputs between different BHPMF attemps and with Genus-Mean gap-filling.

First I run a series of BHPMF Gap-filling with different inputs:
```
From 1_trait_filling.test.R:
We begin with the base unimputed trait matrix, which has a single 'oberved'
trait value for each species, with gaps (NAs).
The main script generates:
1. Matrix gap-filled using genus means.
2. Matrix gap-filled using BHPMF with our three traits of interest, extracting
   from the BHPMF output only the estimates for the gaps in the original matrix.
Here, we generate:
3. Matrix gap-filled using BHPMF with our three traits of interest,
   using the full output of BHPMF (including estimates of 'observed' trait data)
4. Same as #2, but including a fourth dummy trait, with all values = 1
5. Same as #2, but including the following additional traits:
   MaxStemDia_cm, MaxLeafNumber, Max_Rachis_Length_m, Max_Petiole_length_m,
   AverageFruitWidth_cm
```
The next thing to do is comparing the outputs of these five gap-filling attempts. The code for that will be in 2_trait_filling_comparison_test.R, and I'm making a document trait_filling_comparison.docx as an update. To import graphs into a writer document it is best to export as SVG rather than PDF, so I made a custom GraphSVG function first. UPDATE: I might end up using SVGs from now on, cause they're way more awesome than PDFs.

At the end of the day I have some neat combined boxplots, comparing for each
trait the estimated values for each of the five methods. 

## Monday November 26, 2018
This morning I quickly finished making comparison graphs, mostly scatterplots comparing gap-filling method two (standard BHPMF) with the four other methods, for each trait. Subsequently I met up with my daily supervisor again to discuss gap-filling methods. 

Some observations:
For our purposes, boxplots really aren't that informative. Maybe histograms would be better. Then again, we should just move on.
For the comparisons, it looks like there is no strong bias for any of the five methods, i.e. the results make sense, variation is symmetric, et cetera.

There was one strange outlier value in method three. Let's investigate:
```
# Original data:
> traits[which(traits$species == "Nypa_fruticans"), ]
            species genus stem.height blade.length fruit.length
1954 Nypa_fruticans  Nypa           0         10.2         11.5
# Input to method three:
> trait.matrix[which(rownames(trait.matrix) == "Nypa_fruticans"), ]
 stem.height blade.length fruit.length 
      0.0001      10.2000      11.5000
> hierarchy.matrix[which(hierarchy.matrix[, 2] == "Nypa"), ]
         species            genus            tribe        subfamily 
"Nypa_fruticans"           "Nypa"      "Nypoideae"      "Nypoideae"
> hierarchy.matrix[which(hierarchy.matrix[, 3] == "Nypoideae"), ]
         species            genus            tribe        subfamily 
"Nypa_fruticans"           "Nypa"      "Nypoideae"      "Nypoideae" 
> hierarchy.matrix[which(hierarchy.matrix[, 4] == "Nypoideae"), ]
         species            genus            tribe        subfamily 
"Nypa_fruticans"           "Nypa"      "Nypoideae"      "Nypoideae" 
# Result of method three:
> filled.three[which(filled.three$species == "Nypa_fruticans"), ]
            species genus stem.height blade.length fruit.length
1697 Nypa_fruticans  Nypa     -0.0017     -47.1934      11.5002
> test.std[which(test.std$species == "Nypa_fruticans"), ]
            species genus stem.height blade.length fruit.length
1697 Nypa_fruticans  Nypa      0.0955      19.3328       0.0966
```
That estimated blade length value is weird. The uncertainty (stdev) in that blade length estimate is also huge. From the input, we can see it is the only species in its genus, tribe and subfamily, interestingly enough.



There are three things I should look into now:
1) Explicitly compare BHPMF estimates for 'known' values with those known values.
2) To get the fullest coverage, we can combine BHPMF with genus-means. There are two ways to do that: taking genus means from the original 'known' values, or taking genus means from the known + BHMPMF estimated values. These two methods will have to be compared.
3) Compare the use of 1 dummy trait with using 5 dummy traits, to assess how using dummies influences the output.
And also, if there is time, you can divide the data into training and test subsets (95% vs 5% of data points), to assess how accurate the gap-filling methods are.

Aside from gap-filling, there is writing to consider. As inspired by a Science mag article my assessor linked us (DOI: 10.1126/science.353.6300.718), I'm going to spend the first hour of my working days on the writing. Starting tomorrow morning.

Update:
I finished 1). Turns out the BHPMF estimates for 'observed' values can differ substantially from those observed values. Moreover, there is an asymmetrical trend: the higher the observed trait value, the lower the estimated trait value tends to be. All three traits show this pattern, but it is most pronounced for stem height.

Update: I finished 3). It looks like adding dummy variables does not introduce bias in the trait value estimates. But does this workaround provide accurate results for those species where BHPMF normally does not work? It seems so. Gap-filling using BHPMF with a dummy trait is a serious option to consider. 

## Wednesday November 28, 2018
It's maintenance day: I am reviewing my code formatting and comment structure and checking if my scripts and stuff need to be reorganized. I have taken some tips from the Tidyverse style guide (https://style.tidyverse.org), and I'm rewriting my code to use piping operators from package magrittr. Once again I wonder why the master's course in data analysis never took a few days to teach the students how to actually code properly in R. Piping makes code way easier to read, and actually makes coding easier and less confusing as well.

Current situation:
Different gap-filling methods are implemented in test_1_trait_filling.R;
and these are compared in test_2_trait_filling_comparison.R.
The main trait filling script is where the best/chosen trait filling method will be implemented. Consider that script a placeholder for the moment.

