# Palm Project 01: exploring data
# In which the provided datasets are explored and filtered.

setwd("/home/delta/R_Projects/palm_FD")


############################
# 1: Palm distribution data
############################

palmdist <- read.csv(file="palms_in_tdwg3.csv",header=TRUE)
# Data frame with two columns:
#	1. Area_code_L3 - 3-letter code for each botanical country
#	2. SpecName - Species name
#
# This is a list of all species present in each tdwg3 unit.
# There are 198 unique tdwg3 units and 2557 unique species.
# Area code is in alphabetical order.


# One initial objective with this data, is counting species richness for each tdwg3 unit. 
# How many botanical countries have low (<4) species richness?
# Calculation of FD indices requires more species than traits.
# We will use 3 traits, therefore only tdwg3 units with at least 4 species can be used. 

# list richness for each tdwg3 unit:
richness <- data.frame(
	tdwg3=levels(palmdist$Area_code_L3),
	richness=as.numeric(by(palmdist, palmdist$Area_code_L3, nrow))
	)
# to answer the question:
length(which(richness$richness<4))
# 65 (out of 198) botanical countries have richness < 4
# Which leaves 133 usable countries.


# Presence-absence matrix
#########################

# We want to transform the palmdist list into a presence/absence matrix,
# Where rows are tdwg3 units and columns are species.
# Because that is what we need for FD index calculation.

# Don't panic or spass out about needing to 'melt' or 'cast' the data.
# The solution is glaringly simple:
presabs <- as.data.frame.matrix(table(palmdist[1:2]))
# And add twdg3 unit labels:
presabs <- data.frame(
	AreaCodeL3=levels(palmdist$Area_code_L3),
	presabs
	)

# Verification: do the rowsums match the richness we calculated above?
as.numeric(richness$richness - rowSums(presabs[,-1]))
# All values 0, which means we can report success.

# Write to .csv:
write.csv(presabs,file="palm_tdwg3_presabs_matrix.csv",eol="\r\n",row.names=FALSE)





######################
# 2. Palm Traits Data
######################

traitdata <- read.csv(file="PalmTraits_10.csv",header=TRUE)
# Data frame with 31 columns:
#	1. SpecName - binomial species name
#	2. accGenus - only the genus name
#	3. accSpecies - only the species name
#	4. PalmTribe
#	5. PalmSubfamily
#	6-30 - information on 25 traits
#	31. extra reference notes

# Integrity check: are the species in palmdist and traitdata the same?
str(palmdist$SpecName)
str(traitdata$SpecName)
# both are a factor with 2557 levels, which suggests a match.


# Three quantitative traits will be used for this study: stem height, leaf size and fruit size. In addition, growth form is of interest.

# MaxStemHeight_m is the trait for stem height.

# Three traits relate to leaf size: Max_Blade_Length_m, Max_Rachis_Length_m and Max_Petiole_length_m.
# We have to choose one.
# Göldel et al. (2015) used Maximum rachis length.
# Compare correlations:
cor(traitdata[,c("Max_Blade_Length_m","Max_Rachis_Length_m","Max_Petiole_length_m")], use="pairwise.complete.obs")
# Blade length and rachis length are highly correlated (r=0.95),
# But petiole length is only weakly correlated to either (r=0.52 with blade length,
# r=0.51 with rachis length)
# Compare missing data:
length(which(is.na(traitdata$Max_Blade_Length_m)))/length(traitdata$Max_Blade_Length_m)
# 26% missing data
length(which(is.na(traitdata$Max_Rachis_Length_m)))/length(traitdata$Max_Rachis_Length_m)
# 42% missing data
# Given the difference in missing data and the strong correlation,
# the best choice should be maximum blade length.
# So what about the terminology?
# Rachis is the extension of the petiole on a pinnate leaf (Dransfield et al. 2008).
# Does that mean that rachis length does not apply to palmate leaves?
# That would explain the difference in % missing data.
# Turns out this is not true, because costapalmate leaves are common,
# which do have a sort of rachis. 

# Fruitsize is coded by minimum, maximum and average values for length and width.
# Göldel et al. (2015) calculated a fruit volume metric.
# We can look into doing the same at a later stage.
cor(traitdata[,c("AverageFruitLength_cm","AverageFruitWidth_cm")],use="pairwise.complete.obs")
# Length and width are quite strongly correlated (r = 0.88)
# For now I will retain both variables.
# Final check: is there insane intraspecific variation in size variables?
boxplot(traitdata$MaxFruitLength_cm / traitdata$MinFruitLength_cm)
summary(traitdata$MaxFruitLength_cm / traitdata$MinFruitLength_cm)
# range 0.58-4.00
# no extreme outliers. 
# Strange that the max/min ratio for some species is <1,
# meaning the minimum value is larger than the maximum value. 
boxplot(traitdata$MaxFruitWidth_cm / traitdata$MinFruitWidth_cm)
summary(traitdata$MaxFruitWidth_cm / traitdata$MinFruitWidth_cm)
# range 0.68-10.00
# Again for some species the max/min ratio is <1.
# Two outliers on the high end. Mean is 1.475.
quantile((traitdata$MaxFruitWidth_cm / traitdata$MinFruitWidth_cm),.99,na.rm=TRUE)
# 99th percentile is 3.2.
# So actually, the variation is not that extreme.
# Which species are the outliers?
traitdata[which((traitdata$MaxFruitWidth_cm / traitdata$MinFruitWidth_cm)>6),]
# Livistona carinensis, fruit width 0.5-5
# Pelagodoxa henryana, fruit width 2.0-15

# Growth form is coded as three traits: Climbing, Acaulescence and Errect.
# Explore:
summary(as.factor(traitdata$Climbing))
summary(as.factor(traitdata$Acaulescence))
summary(as.factor(traitdata$Errect))
# There are a few NAs in Acaulescence and Errect.
# Also, in addition to binary 0 and 1 values, some entries have value 2.
# The expected pattern is that one of the three traits has value 1,
# and the other two have value 0.
# We can list the exceptions:
traitdata[which(
	(traitdata$Climbing + traitdata$Acaulescence + traitdata$Errect)!=1
	)
,c("SpecName","Climbing","Acaulescence","Errect")]
# The pattern is consistent: these species have a value 2 for two of the three growthform traits.
# Perhaps this indicates the growthform classification was difficult.
# I will have to ask. Until then we'll call these indeterminate.
# One single species (Licuala bruneiana) has two 1 values. This could be a typo
# according to https://plants.jstor.org/stable/10.5555/al.ap.specimen.k000697823
# it is stemless (i.e. acaulescent), so I will classify it as such.
# Next we must investigate the NAs:
traitdata[which(is.na(traitdata$Acaulescence)),c("SpecName","Climbing","Acaulescence","Errect")]
traitdata[which(is.na(traitdata$Errect)),c("SpecName","Climbing","Acaulescence","Errect")]
# These lists overlap almost entirely.
# for traitdata$Errect, the first two species with value NA have value 2 for Climbing + Acaulescence

# I will combine the three growthform variables into one factor.
# For now, species with missing values or with value 2 will be set to NA.
# Licuala bruneiana is set to acaulescent.
GrowthForm<-1:length(traitdata[,1])
GrowthForm[which(traitdata$Climbing==1)] <- "climbing"
GrowthForm[which(traitdata$Acaulescence==1)] <- "acaulescent"
GrowthForm[which(traitdata$Errect==1)] <- "freestanding"
GrowthForm[which(
	(traitdata$Climbing + traitdata$Acaulescence + traitdata$Errect)!=1
	)] <- NA
GrowthForm[which(is.na(traitdata$Errect))] <- NA
GrowthForm[which(traitdata$SpecName=="Licuala bruneiana")] <- "acaulescent"
# Verify:
summary(as.factor(GrowthForm))
# confirmed, all values are as expected.


# Finally, Make a smaller dataframe with only species name and the selected traits.
palmtraits <- data.frame(
	SpecName=traitdata$SpecName,
	GrowthForm=as.factor(GrowthForm),
	MaxStemHeight_m=traitdata$MaxStemHeight_m,
	MaxBladeLength_m=traitdata$Max_Blade_Length_m,
	AverageFruitLength_cm=traitdata$AverageFruitLength_cm,
	AverageFruitWidth_cm=traitdata$AverageFruitWidth_cm
	)
summary(palmtraits)
# It's all there :)
# Write to .csv:
write.csv(palmtraits,file="palmtraits_selected.csv",eol="\r\n",row.names=FALSE)



#########################
# 3. Environmental data #
#########################

# Environmental data was supplied to me in a Excel file.
# I saved this as .csv under the same name.

envdata <- read.csv(file="TDWG_Environment_AllData_2014Dec.csv",header=TRUE)
# Comparison with original Excel file shows a match.
# Do check each variable of interest for consistency, especially for NA values.

# The dataset contains a lot of variables, not all of which we need.
# I will seperate variables of interest into their own data frames.


#----------------------------
# First is TDWG3 information.
#----------------------------
countries <- data.frame(
	AreaCodeL3=envdata$LEVEL_3_CO,	# 3-letter code for each botanical country
	LevelName=envdata$LEVEL_NAME,	# name of botanical country
	Realm=envdata$THREEREALM,	# factor with levels "NewWorld", "OWEast", "OWWest"
	FloraRegion=envdata$REALM_LONG,	# Floristic Region (factor with 8 levels)
	Lat=envdata$LAT,	# Latitude
	Lon=envdata$LONG,	# Longitude
	IsIsland=as.factor(envdata$ISISLAND), # binary: 0 = mainland, 1 = island
	HasPalms=as.factor(envdata$PalmPresAbs),	# Binary: 0 = no palms, 1 = palms present
	PalmRichness=envdata$PALMSR		# Palm species richness
	)
# Integrity check:
summary(countries)
summary(duplicated(countries[,1]))	# only FALSE, so no duplicate codes.
# LevelName has two empty entries. This is not a problem.
# Realm has one empty value, belonging to Antarctica. We have to class this as NA.
countries[which(countries$Realm==""),"Realm"] <- NA
# and drop the empty level
countries <- droplevels(countries)
summary(countries)
# all fine.

# Next we should compare Area codes, HasPalms and PalmRichness
# With our palm presence/absence dataset, in dataframe 'richness'

# HasPalms:
length(richness$tdwg3)						# 198
length(which(countries$HasPalms==1))		# 196. Whoops! Is this data older?
length(which(countries$PalmRichness>0))		# 196. At least that matches.

# PalmRichness:
# First, give 'richness' identical column names
names(richness) <- c("AreaCodeL3","PalmRichness")
# Then compare:
summary(richness$PalmRichness)
# range 1 - 307, median 7, mean 27.04
summary(countries[which(countries$PalmRichness>0),]$PalmRichness)
# range 1 - 275, median 6, mean 25.64
# Clearly the dataset in 'countries' has greater richness.
# Is that data more up-to-date?
# It was supplied to me as the palm distribution data to use, so I will use that.

# Palm richness data from 'richness' is more complete, so we must
# override the IsIsland and HasPalms variables with data from 'richness'.
# First, coerce 'richness' to the same length as 'countries'
richness <- rbind(richness, data.frame(
		AreaCodeL3=rep(NA,length(countries$AreaCodeL3)-length(richness$AreaCodeL3)),
		PalmRichness=rep(NA,length(countries$AreaCodeL3)-length(richness$AreaCodeL3))
		)
	)
# Second, merge by "AreaCodeL3"
a <- merge(countries,richness,by="AreaCodeL3",all.x=TRUE)
# Note this new dataframe has two columns for PalmRichness.
# That is helpful for integrity checks.
summary(a)
# Note that PalmRichness.y has 171 NAs. That is one too many!
length(which(a$PalmRichness.y>0))
# [1] 197. That should be 198. One of the original values is lost.
# Find out which one:
b <- data.frame(
	richness[order(richness$PalmRichness,na.last=TRUE),],
	a[order(a$PalmRichness.y,na.last=TRUE),c("AreaCodeL3","PalmRichness.y")]
	)
b <- data.frame(b, 	diff= b$PalmRichness - b$PalmRichness.y)
# and use your eyes to scan from top to bottom:
b
# In 'richness', twdg3 code goes UGA > VNA > WHM.
# After the merge, the code goes UGA > WHM. So VNA disappeared.
which(richness$AreaCodeL3=="VNA")
# [1] 189
which(countries$AreaCodeL3=="VNA")
# integer(0)	
# WELL THERE IS YOUR PROBLEM! One of the TDWG3 units must be misnamed.
# Investigate! I downloaded the latest tdwg3 database from https://github.com/tdwg/wgsrpd/tree/master/level3
# This database does have a tdwg3 unit called VNA: Venezualan Antilles.
# The database supplied to me by Jun/Daniel however does NOT contain VNA.
# I am forced to conclude that they gave me outdated files.
# Which means I do not have environmental data for VNA either.
# Bottom line: we keep the merged dataframe as-is. We could restore VNA, but we don't have environmental data for it, so there is no point.
# Btw, the palm richness of VNA is 6, so at least we are not missing a major biodiversity hotspot in our study.

# Back to assembling data:
# In the merged dataframe, the PalmRichness data to use is PalmRichness.y.
# HasPalms should be redefined based on it.
# So we overwrite a accordingly:
c <- data.frame(
	a[,c("AreaCodeL3","LevelName","Realm","FloraRegion","Lat","Lon","IsIsland")],
	HasPalms=rep(0,length(a[,1])),
	PalmRichness=a$PalmRichness.y
	)
# Set NA values in PalmRichness to 0, because we know they are 0:
c[which(is.na(c$PalmRichness)),"PalmRichness"] <- 0
# Add 1 values to HasPalms:
c[which(c$PalmRichness>0),"HasPalms"] <- 1

# Done. Write to .csv file:
write.csv(c,file="tdwg3_overview.csv",eol="\r\n",row.names=FALSE)

# And because it might be useful, we want a .csv with only entries
# where PalmRichness > 3, i.e. the subset of data with high enough richness
# for calculating FD indices.
d <- c[which(c$PalmRichness>3),]
write.csv(d,file="tdwg3_overview_PalmsOnly.csv",eol="\r\n",row.names=FALSE)


#---------------------------------------
# Second is non-climatic environmentals.
#---------------------------------------

environ <- data.frame(
	AreaCodeL3=envdata$LEVEL_3_CO,	# 3-letter code for each botanical country
	Area_km2=envdata$AREA_KM2,		# surface area of tdwg3 unit
	AltMean=envdata$srtm_alt_mean,	# Mean surface altitude
	AltRange=envdata$alt_range,		# Altitudinal range
	SoilCount=envdata$soilcount,	# Number of soil types present in tdwg3 unit
	CanHeightMean=envdata$CH_Mean,	# mean canopy height
	CanHeightMin=envdata$CH_Min,	# Minimum canopy height
	CanHeightMax=envdata$CH_Max,	# Maximum canopy height
	CanHeightRange=envdata$CH_Range # Canopy height range (Max - Min)
	)
summary(environ)
# Not bad. A missing value here and there, particularly in canopy height
# Regardless, we first want to subset this to include only tdwg3 units
# where palms occur.
# This we do via merging with the tdwg3 description dataframe
countries <- read.csv(file="tdwg3_overview_PalmsOnly.csv",header=TRUE)
a <- merge(countries,environ,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
a <- a[,c(1,10:17)]
length(a[,"AreaCodeL3"])
# [1] 132	# That's good.
summary(a)
# Still some missing values in canopy height, but not in the other variables.
# So for 14 out of 132 botanical countries, we do not have canopy data.
# That is an issue for later though.

# Write to .csv:
write.csv(a,file="tdwg3_environment.csv",eol="\r\n",row.names=FALSE)


#-------------------------------
# Third is contemporary climate.
#-------------------------------

climate <- data.frame(
	AreaCodeL3=envdata$LEVEL_3_CO,	# 3-letter code for each botanical country
	PrecAnnual=envdata$PREC_Sum,	# (mean?) Annual precipitation
	PrecSeas=envdata$PREC_CV,		# Precipitation Seasonality
	PrecDryQ=envdata$P_drie_quart,	# Precipitation of driest quarter
	MAT=envdata$Tmean_mean,			# Mean Annual Temperature
	TColdQ=envdata$T_cold_quart,	# (mean?) Temperature of coldest quarter
	TempSeas=envdata$Temp_SD,		# Temperature seasonality
	TColdM=envdata$Tmin_cold_month	# Temperature of coldest month
	)
# How do you choose which variables to include? Use your wits.
# As a rule of thumb, you want to include all unique climatic variables.
# Turns out we just use everything :)
# In the end it's only 7 variables. Not many, but it includes the most important
# ones.

# As above, we subset to tdwg3 units where PalmRichness > 3
countries <- read.csv(file="tdwg3_overview_PalmsOnly.csv",header=TRUE)
a <- merge(countries,climate,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
a <- a[,c(1,10:16)]
length(a[,"AreaCodeL3"])
# [1] 132	# That's good.
summary(a)
# No missing values. We can rejoice.

# Write to .csv:
write.csv(a,file="tdwg3_climate.csv",eol="\r\n",row.names=FALSE)


#---------------------------------
# Fourth is LGM climate anomalies.
#---------------------------------

# Keep in mind we are not looking at historical climate,
# but historical climate CHANGE, in the form of anomalies between LGM and present.

anomalies <- data.frame(
	AreaCodeL3=envdata$LEVEL_3_CO,	# 3-letter code for each botanical country
	PrecAnom=envdata$ensLGM_Pano,	# Precipitation anomaly LGM <-> present
									# ensemble mean from CCSM3 + MIROC3.2
	TAnom=envdata$ensLGM_Tano,		# Temperature anomaly LGM <-> present
									# ensemble mean from CCSM3 + MIROC3.2
	PrecAnomCCSM=envdata$ccAnoPmean,# mean prec anom from CCSM3
	TAnomCCSM=envdata$ccAnoTmean,	# mean T anom from CCSM3
	PrecAnomMIROC=envdata$miAnoPmean,	# mean prec anom from MIROC3.2
	TAnomMIROC=envdata$miAnoTmean		#mean T anom from MIROC3.2
	)

# As above, we subset to tdwg3 units where PalmRichness > 3
countries <- read.csv(file="tdwg3_overview_PalmsOnly.csv",header=TRUE)
a <- merge(countries,anomalies,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
a <- a[,c(1,10:15)]
length(a[,"AreaCodeL3"])
# [1] 132	# That's good.
summary(a)
# No missing values. We can rejoice.

# Write to .csv:
write.csv(a,file="tdwg3_climate_anomalies.csv",eol="\r\n",row.names=FALSE)

