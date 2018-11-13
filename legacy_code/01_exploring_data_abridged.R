# Palm Project 01: exploring data
# In which the provided datasets are explored, filtered, and written to .csv files.
#
# The following data files are produced:
# palmtraits_selected.csv - relevant trait data for all palm species.
# palm_tdwg3_presabs.csv - palm presence/absence matrix.
# palm_tdwg3_presabs_subset.csv - same as above but only for PalmRichness > 3.
# tdwg3_overview.csv - descriptors and Palm Richness for all tdwg3 units.
# tdwg3_overview_subset.csv - same as above but only for PalmRichness > 3.
# tdwg3_environment.csv, with relevant environmental variables for tdwg3 units with PalmRichness >3.
# tdwg3_climate.csv, with contemporary climate variables
# tdwg3_climate_anomalies.csv, with LGM <-> present climate anomalies for tdwg3 units with PalmRichness >3.


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
#
# Calculation of FD indices requires more species than traits (S > T).
# Additionally, comparing functional richness values for communities with less than 2^T species is complicated (but not impossible).
#
# We will use 3 traits, therefore only tdwg3 units with at least 4 species can be used, and we'd like to know how many tdwg3 units have richness of 4-8.


# list richness for each tdwg3 unit:
richness <- data.frame(
	tdwg3=levels(palmdist$Area_code_L3),
	richness=as.numeric(by(palmdist, palmdist$Area_code_L3, nrow))
	)
# 'richness' is combined with other data frames below.
# 65 (out of 198) botanical countries have richness < 4.
# 101 (out of 198) botanical countries have richness < 8.
# That means 133 countries are usable, of which 36 have 'low' richness (27%).


# Presence-absence matrix
#########################

# We want to transform the palmdist list into a presence/absence matrix,
# Where rows are tdwg3 units and columns are species.
# Because that is what we need for FD index calculation.

# This is actually really simple:
presabs <- as.data.frame.matrix(table(palmdist[1:2]))
# And add twdg3 unit labels:
presabs <- data.frame(
	AreaCodeL3=levels(palmdist$Area_code_L3),
	presabs
	)

# Write to .csv:
write.csv(presabs,file="palm_tdwg3_presabs.csv",eol="\r\n",row.names=FALSE)


######################
# 2. Palm Traits Data
######################

traitdata <- read.csv(file="PalmTraits_10.csv",header=TRUE)
# Data frame with 31 columns:
#	1. SpecName - binomial species name [sorted alphabetically]
#	2. accGenus - only the genus name
#	3. accSpecies - only the species name
#	4. PalmTribe
#	5. PalmSubfamily
#	6-30 - information on 25 traits
#	31. extra reference notes

# I will combine the three growthform variables into one factor.
# For now, species with missing values or with growthform value 2 will be set to NA.
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

# Finally, Make a smaller dataframe with only species name and the selected traits.
palmtraits <- data.frame(
	SpecName=traitdata$SpecName,
	GrowthForm=as.factor(GrowthForm),
	MaxStemHeight_m=traitdata$MaxStemHeight_m,
	MaxBladeLength_m=traitdata$Max_Blade_Length_m,
	AvgFruitLength_cm=traitdata$AverageFruitLength_cm
	)

# Write palmtraits to .csv:
write.csv(palmtraits,file="palmtraits_selected.csv",eol="\r\n",row.names=FALSE)


#########################
# 3. Environmental data #
#########################

# Environmental data was supplied to me in an Excel file.
# I saved this as .csv under the same name.

envdata <- read.csv(file="TDWG_Environment_AllData_2014Dec.csv",header=TRUE)
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
# Realm has one empty value, belonging to Antarctica. We have to class this as NA.
countries[which(countries$Realm==""),"Realm"] <- NA
# and drop the empty level
countries <- droplevels(countries)

# Palm richness data from 'richness' is more complete, so we must
# override the IsIsland and HasPalms variables with data from 'richness'.
names(richness) <- c("AreaCodeL3","PalmRichness") # for simplicity
# First, coerce 'richness' to the same length as 'countries'
richness <- rbind(richness, data.frame(
		AreaCodeL3=rep(NA,length(countries$AreaCodeL3)-length(richness$AreaCodeL3)),
		PalmRichness=rep(NA,length(countries$AreaCodeL3)-length(richness$AreaCodeL3))
		)
	)
# Second, merge by "AreaCodeL3"
a <- merge(countries,richness,by="AreaCodeL3",all.x=TRUE)
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
write.csv(d,file="tdwg3_overview_subset.csv",eol="\r\n",row.names=FALSE)


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
# We first want to subset this to tdwg3 units with sufficient PalmRichness
# This we do via merging with the tdwg3 description dataframe
countries <- read.csv(file="tdwg3_overview_subset.csv",header=TRUE)
a <- merge(countries,environ,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
a <- a[,c(1,10:ncol(a))]

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
# As above, we subset to tdwg3 units where PalmRichness > 3
countries <- read.csv(file="tdwg3_overview_subset.csv",header=TRUE)
a <- merge(countries,climate,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
a <- a[,c(1,10:ncol(a))]

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
countries <- read.csv(file="tdwg3_overview_subset.csv",header=TRUE)
a <- merge(countries,anomalies,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
a <- a[,c(1,10:ncol(a))]

# Write to .csv:
write.csv(a,file="tdwg3_climate_anomalies.csv",eol="\r\n",row.names=FALSE)


#################################
# short presence/absence matrix #
#################################

# The last thing we need to do is subset the presence/absence matrix
# to include only tdwg3 units with PalmRichness > 3

subpresabs <- merge(countries,presabs,by="AreaCodeL3",all.x=FALSE)
# And subsequent exclusion of columns
subpresabs <- subpresabs[,c(1,10:ncol(subpresabs))]

# Write to .csv:
write.csv(subpresabs,file="palm_tdwg3_presabs_subset.csv",eol="\r\n",row.names=FALSE)


