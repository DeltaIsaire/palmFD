# Palm Project 02: Variable Analysis
# In which the selected palm traits are explored,
# And correlations between GLM climate anomalies from 2 sources are investigated

setwd("/home/delta/R_Projects/palm_FD")

library(plyr)

##############################
# Exploration of palm traits #
##############################

# load selected trait data
palmtraits <- read.csv(file="palmtraits_selected.csv",header=TRUE)

# For each trait we want:
# descriptive statistics: percentage NA, min, max, median, mean
# visualisation (histogram / boxplot)
# Outlier investigation

# dataframe with descriptive statistics
#######################################
# all values except median and max are rounded to 2 decimals for easy reading
traitstats <- data.frame(
	trait=names(palmtraits),
	missing=round(as.numeric(colwise(function(x) sum(is.na(x))/length(x)) (palmtraits)),digits=2),
	min=c(NA,NA,round(as.numeric(numcolwise(function(x) min(x,na.rm=TRUE)) (palmtraits)),digits=2)),
	median=c(NA,NA,as.numeric(numcolwise(function(x) median(x,na.rm=TRUE)) (palmtraits))),
	max=c(NA,NA,as.numeric(numcolwise(function(x) max(x,na.rm=TRUE)) (palmtraits))),
	mean=c(NA,NA,round(as.numeric(numcolwise(function(x) mean(x,na.rm=TRUE)) (palmtraits)),digits=2))
	)

# Pretty Pictures
#################

# Bar percent chart of growthForm
#
# This requires the frequency of growthform, sorted in decreasing order,
# in matrix format.
growthform <- as.matrix(summary(palmtraits$GrowthForm)[order(summary(palmtraits$GrowthForm),decreasing=TRUE)])
# Then we make the barplot
# note use of ylim argument to control height of the bar.
barplot(growthform, width=0.25, horiz=TRUE, col=c("green","red","yellow","black"), xlab="cumulative frequency",ylim=c(0,1))
legend("topright",legend=rownames(growthform),fill=c("green","red","yellow","black"))


# Boxplots of numerical traits
par(mfrow=c(1,4))
#stem height
boxplot(palmtraits$MaxStemHeight_m,notch=FALSE,outline=TRUE,border="black",col="lightgrey",ylab="Maximum stem height (m)",xlab="Maximum stem height",boxwex=1.4)
# Blade length
boxplot(palmtraits$MaxBladeLength_m,notch=FALSE,outline=TRUE,border="black",col="lightgrey",ylab="Maximum blade length (m)",xlab="Maximum blade length",boxwex=1.4)
# Fruit Length
boxplot(palmtraits$AverageFruitLength_cm,notch=FALSE,outline=TRUE,border="black",col="lightgrey",ylab="Average fruit length (cm)",xlab="Average fruit length",boxwex=1.4)
# Fruit Width
boxplot(palmtraits$AverageFruitWidth_cm,notch=FALSE,outline=TRUE,border="black",col="lightgrey",ylab="Average fruit width (cm)",xlab="Average fruit width",boxwex=1.4)
# cancel par
par(mfrow=c(1,1))


# And because it's awesome, boxplots of traits grouped by growthform
# With text labels in top-right of each plot
# four plots grouped 2x2
par(mfrow=c(2,2),mar=c(2,4,1,1)+0.1) # order of mar is bottom, left, top, right
boxplot(MaxStemHeight_m ~ GrowthForm, data=palmtraits,notch=FALSE,outline=TRUE,varwidth=TRUE,border="black", col="lightgrey",ylab="Maximum stem height (m)",boxwex=1)
text(x=3.4,y=150,labels="A",cex=2,pos=1)
boxplot(MaxBladeLength_m ~ GrowthForm, data=palmtraits,notch=FALSE,outline=TRUE,varwidth=TRUE,border="black", col="lightgrey",ylab="Maximum blade length (m)",boxwex=1)
text(x=3.4,y=25,labels="B",cex=2,pos=1)
boxplot(AverageFruitLength_cm ~ GrowthForm, data=palmtraits,notch=FALSE,outline=TRUE,varwidth=TRUE,border="black", col="lightgrey",ylab="Average fruit length (cm)",boxwex=1)
text(x=3.4,y=45,labels="C",cex=2,pos=1)
boxplot(AverageFruitWidth_cm ~ GrowthForm, data=palmtraits,notch=FALSE,outline=TRUE,varwidth=TRUE,border="black", col="lightgrey",ylab="Average fruit width (cm)",boxwex=1)
text(x=3.4,y=20,labels="D",cex=2,pos=1)
par(mfrow=c(1,1))

# Calamus asthonii, the acaulescent palm with stem height of 20 m.
# As far as outliers go, this is the only whose value I doubt.


################################
# Correlation of LGM Anomalies #
################################

# read anomalies file
anomalies <- read.csv(file="tdwg3_climate_anomalies.csv",header=TRUE)

# We have anomalies from two different LGM climate sources:
# CCSM3 and MIROC3.2

# The main question I want to answer:
# Is there a reason not to use the ensemble mean anomalies of both sources?

# Precipitation anomalies
#########################

cor(anomalies[,c("PrecAnom","PrecAnomCCSM","PrecAnomMIROC")])
#               PrecAnom PrecAnomCCSM PrecAnomMIROC
# PrecAnom      1.0000000    0.3173573     0.9227263
# PrecAnomCCSM  0.3173573    1.0000000    -0.0726687
# PrecAnomMIROC 0.9227263   -0.0726687     1.0000000
#
# There are certainly some differences.
# CCSM and MIROC are barely correlated, indicating difference in predicted values.
# In addition, the ensemble mean is related much more closely to MIROC data,
# which suggests that MIROC had a better sample size?

cor(anomalies[,c("TAnom","TAnomCCSM","TAnomMIROC")])
#                TAnom TAnomCCSM TAnomMIROC
# TAnom      1.0000000 0.9671960  0.8319237
# TAnomCCSM  0.9671960 1.0000000  0.6639181
# TAnomMIROC 0.8319237 0.6639181  1.0000000
#
# For Temperature anomalies, the differences are less extreme,
# But CCSM and MIROC still do not completely agree.

# This is exactly the sort of thing to ask my supervisors.

# The data can also be visualized:
par(mfrow=c(1,2))
plot(PrecAnomCCSM~PrecAnomMIROC,data=anomalies,main="Precipitation anomaly")
plot(TAnomCCSM~TAnomMIROC,data=anomalies,main="Temperature anomaly")
par(mfrow=c(1,1))



############################################
# Spatial maps of tdwg3 units (with palms) #
############################################

# Consider the maptools library
# https://www.r-bloggers.com/r-and-gis-working-with-shapefiles/
# library rgdal has a better polygon read function: readOGR()
library(maptools)
library(rgdal)

map <- readOGR(dsn="/home/delta/R_Projects/palm_FD/tdwg", layer="TDWG_level3_Coordinates")
# Do not specify file extension. I am not sure which of the files it actually uses,
# but it works.



