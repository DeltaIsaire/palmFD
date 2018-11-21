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
