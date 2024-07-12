setwd("C:/Users/jpt93/Documents/R/Occu-abu model/case_study")
library(AHMbook)
library(jagsUI)
#detection - non detection data
operationmatrix <- read.csv("allNY_2013-2021_7dayocc_americanmarten_detection_nondetection.csv")
yC  <- read.csv("allNY_2013_2021_7dayocc_coyote_counts.csv")
yF <- read.csv("allNY_2013-2021_7dayocc_fisher_counts.csv")
yM <- read.csv("7daycountDetectionsmartenNZ.csv")


# Select NZ years only
NZyears2016 <- which(operationmatrix$year==2016)
NZyears2017 <- which(operationmatrix$year==2017)
NZyears2018 <- which(operationmatrix$year==2018)

NZyears <- c(NZyears2016, NZyears2017, NZyears2018)

operationmatrix <- operationmatrix[NZyears,]
yC <- yC[NZyears,]
yF <- yF[NZyears,]

#remove sites that are all NAs (e.g. not sampled)
row.names(operationmatrix)<- 1:nrow(operationmatrix)
row.names(yC)<- 1:nrow(yC)
row.names(yF)<- 1:nrow(yF)

nacount <- apply(operationmatrix, MARGIN = 1, function(x) sum(is.na(x)))
allnas <- which(nacount == 3)
yM <- yM[-allnas,]
yC <- yC[-allnas,]
yF <- yF[-allnas,]
operationmatrix <- operationmatrix[-allnas,]
y <- operationmatrix[,4:6]


sitenames=unique(yF$sitename)
years=unique(yF$year)
n.year=length(years)
J=rep(NA,n.year)
year.sites=list(3)
for(i in 1:n.year){
  J[i]=length(yF$sitename[yF$year==years[i]])
  year.sites[[i]]=unique(yF$sitename[yF$year==years[i]])
}
maxJ=max(J)
rep(c(0, 1, 2), each = 10)
K=3


#format y into 3D array sites, occasions, years
yF.new=array(0,dim=c(maxJ,K,n.year)) #new data
for(i in 1:nrow(yF)){ #loop through each row
  this.year=which(years==yF$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==yF$sitename[i]) #get site for this row
  yF.new[this.site,1:3,this.year]=as.numeric(yF[i,3:5]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(yF.new[1:J[t],,t], na.rm=TRUE),sum(yF[yF$year==years[t],3:5], na.rm=TRUE)))
}

#format y into 3D array sites, occasions, years
yC.new=array(0,dim=c(maxJ,K,n.year)) #new data
for(i in 1:nrow(yC)){ #loop through each row
  this.year=which(years==yC$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==yC$sitename[i]) #get site for this row
  yC.new[this.site,1:3,this.year]=as.numeric(yC[i,2:4]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(yC.new[1:J[t],,t], na.rm=TRUE),sum(yC[yC$year==years[t],2:4], na.rm=TRUE)))
}


#format y into 3D array sites, occasions, years
yM.new=array(0,dim=c(maxJ,K,n.year)) #new data
for(i in 1:nrow(yM)){ #loop through each row
  this.year=which(years==yM$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==yM$sitename[i]) #get site for this row
  yM.new[this.site,1:3,this.year]=as.numeric(yM[i,2:4]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(yM.new[1:J[t],,t], na.rm=TRUE),sum(yM[yM$year==years[t],2:4], na.rm=TRUE)))
}

K2D <- 1*(!is.na(y))

K2D <- as.data.frame(K2D)
K2D$year <- operationmatrix$year
K2D$sitename <- operationmatrix$sitename

#format y into 3D array sites, occasions, years
K2D.new=array(0,dim=c(maxJ,K,n.year)) #new data
for(i in 1:nrow(K2D)){ #loop through each row
  this.year=which(years==K2D$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==K2D$sitename[i]) #get site for this row
  K2D.new[this.site,1:3,this.year]=as.numeric(K2D[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(K2D.new[1:J[t],,t], na.rm=TRUE),sum(K2D[K2D$year==years[t],1:3], na.rm=TRUE)))
}
idx <- which(K2D.new == 0)

yM.new[idx] <- 0
yF.new[idx] <- 0
yC.new[idx] <- 0

#create occassion matrix
occ<- y
for(i in 1:nrow(y)){
  tmp<- y[i,]
  n.samp<- sum(!is.na(tmp))
  tmp[!is.na(tmp)]<- 1:n.samp
  occ[i,]<- tmp
}



occ <- as.matrix(occ)
occ.c = standardize(occ)

occyear <- cbind(occ.c, yM)

#format occ covs
occ.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(occ)){ #loop through each row
  this.year=which(years==occyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==occyear$sitename[i]) #get site for this row
  occ.new[this.site,1:3,this.year]=as.numeric(occyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(occ.new[1:J[t],,t], na.rm=TRUE),sum(occyear[occyear$year==years[t],1:3],na.rm=TRUE)))
}

occ.new[is.na(occ.new)] <- 0


#dates
#dates
#read in ordinal day data, remove NAs and format
date <- read.csv('juliandays_allNY_2013-2021.csv')
date <- date[NZyears,]
date <- date[-allnas,]

date <- date[3:5]
date <- as.matrix(date)
date.c = standardize(date)

dateyear <- cbind(date.c, yM)
# date[is.na(y)] <- NA # fix the NAs

#for week with year last
date.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(dateyear)){ #loop through each row
  this.year=which(years==dateyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==dateyear$sitename[i]) #get site for this row
  date.new[this.site,1:3,this.year]=as.numeric(dateyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(date.new[1:J[t],,t]),sum(dateyear[dateyear$year==years[t],1:3])))
}


# temp

temp <- read.csv('NZAvgTempByOcc.csv')
temp <- temp[-allnas,]
temp <- temp[4:6]
temp <- as.matrix(temp)
temp.c = standardize(temp)
tempyear <- cbind(temp.c, yM)
# date[is.na(y)] <- NA # fix the NAs

#for week with year last
temp.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(tempyear)){ #loop through each row
  this.year=which(years==tempyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==tempyear$sitename[i]) #get site for this row
  temp.new[this.site,1:3,this.year]=as.numeric(tempyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(temp.new[1:J[t],,t]),sum(tempyear[tempyear$year==years[t],1:3])))
}

# bait

bait <- read.csv('NZ7dayoccbait.csv')
bait <- bait[-allnas,]
bait <- bait[4:6]

baits01 <- model.matrix( ~ bait[,1] - 1, data=bait )
baits02 <- model.matrix( ~ bait[,2] - 1, data=bait )
baits03 <- model.matrix( ~ bait[,3] - 1, data=bait )

deer <- cbind(baits01[,2], baits02[,2], baits03[,2])
beaver <- cbind(baits01[,1], baits02[,1], baits03[,1])
misc <- cbind(baits01[,3], baits02[,3], baits03[,3])
moose <- cbind(baits01[,4], baits02[,4], baits03[,4])

deer <- as.matrix(deer)
deeryear <- cbind(deer, yM)

#for week with year last
deer.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(deeryear)){ #loop through each row
  this.year=which(years==deeryear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==deeryear$sitename[i]) #get site for this row
  deer.new[this.site,1:3,this.year]=as.numeric(deeryear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(deer.new[1:J[t],,t]),sum(deeryear[deeryear$year==years[t],1:3])))
}

beaver <- as.matrix(beaver)
beaveryear <- cbind(beaver, yM)

#for week with year last
beaver.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(beaveryear)){ #loop through each row
  this.year=which(years==beaveryear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==beaveryear$sitename[i]) #get site for this row
  beaver.new[this.site,1:3,this.year]=as.numeric(beaveryear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(beaver.new[1:J[t],,t]),sum(beaveryear[beaveryear$year==years[t],1:3])))
}

moose <- as.matrix(moose)
mooseyear <- cbind(moose, yM)

#for week with year last
moose.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(mooseyear)){ #loop through each row
  this.year=which(years==mooseyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==mooseyear$sitename[i]) #get site for this row
  moose.new[this.site,1:3,this.year]=as.numeric(mooseyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(moose.new[1:J[t],,t]),sum(mooseyear[mooseyear$year==years[t],1:3])))
}

misc <- as.matrix(misc)
miscyear <- cbind(misc, yM)

#for week with year last
misc.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(miscyear)){ #loop through each row
  this.year=which(years==miscyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==miscyear$sitename[i]) #get site for this row
  misc.new[this.site,1:3,this.year]=as.numeric(miscyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(misc.new[1:J[t],,t]),sum(miscyear[miscyear$year==years[t],1:3])))
}
#read in sitecovs

siteCovs15km <- read.csv("allNY__2013-2021_15kmbuffer_allsitecovs.csv")
siteCovs15km <- siteCovs15km[NZyears,]
siteCovs15km <- siteCovs15km[-allnas,]

siteCovs6km <- read.csv("alNY_2013-2021_6kmbuffer_allsitecovs.csv")
siteCovs6km <- siteCovs6km[NZyears,]
siteCovs6km <- siteCovs6km[-allnas,]
head(siteCovs6km)
#create aggregated site covariates
siteCovs15km[,"Open"] <- siteCovs15km[,"Pasture"] + siteCovs15km[,"Cultivated.Crops"] + siteCovs15km[,"Grassland"] + siteCovs15km[,"Barren"] + siteCovs15km[,"Scrub"]
siteCovs15km[,"All.Developed"] <- siteCovs15km[,"Developed.Low"] + siteCovs15km[,"Developed.Mid"] + siteCovs15km[,"Developed.High"]
siteCovs15km[,"Wetlands"]<- siteCovs15km[,"Woody.Wetlands"] + siteCovs15km[,"Herbaceus.Wetlands"]
siteCovs15km[,"Agriculture"] <- siteCovs15km[,"Pasture"] + siteCovs15km[,"Cultivated.Crops"] 
siteCovs15km[,"ConifMixed"] <- siteCovs15km[,"Coniferous"] + siteCovs15km[,"Mixed"] 
siteCovs15km[,"All.Forest"] <- siteCovs15km[,"Coniferous"] + siteCovs15km[,"Mixed"] + siteCovs15km[,"Deciduous"]

siteCovs6km[,"Open"] <- siteCovs6km[,"Pasture"] + siteCovs6km[,"Cultivated.Crops"] + siteCovs6km[,"Grassland"] + siteCovs6km[,"Barren"] + siteCovs6km[,"Scrub"]
siteCovs6km[,"All.Developed"] <- siteCovs6km[,"Developed.Low"] + siteCovs6km[,"Developed.Mid"] + siteCovs6km[,"Developed.High"]
siteCovs6km[,"Wetlands"]<- siteCovs6km[,"Woody.Wetlands"] + siteCovs6km[,"Herbaceus.Wetlands"]
siteCovs6km[,"Agriculture"] <- siteCovs6km[,"Pasture"] + siteCovs6km[,"Cultivated.Crops"] 
siteCovs6km[,"ConifMixed"] <- siteCovs6km[,"Coniferous"] + siteCovs6km[,"Mixed"] 
siteCovs6km[,"All.Forest"] <- siteCovs6km[,"Coniferous"] + siteCovs6km[,"Mixed"] + siteCovs6km[,"Deciduous"]


nsites = nrow(y)
nsurveys = dim(y)[2]
conifmixed15km = array(siteCovs15km$ConifMixed, dim=c(nsites))
deciduous15km = array(siteCovs15km$Deciduous, dim=c(nsites))
deer15km = array(siteCovs15km$deer, dim= c(nsites))
wetlands15km = array(siteCovs15km$Wetlands, dim=c(nsites))
crops15km = array(siteCovs15km$Cultivated.Crops, dim= c(nsites))
forest15km = array(siteCovs15km$All.Forest, dim=c(nsites))
road15km = array(siteCovs15km$road_density, dim=c(nsites))
edge15km = array(siteCovs15km$forest_edgedensity, dim=c(nsites))
snowdepth15km = array(siteCovs15km$snow_depth, dim=c(nsites))
conifmixed6km = array(siteCovs6km$ConifMixed, dim=c(nsites))
deciduous6km = array(siteCovs6km$Deciduous, dim=c(nsites))
wetlands6km = array(siteCovs6km$Wetlands, dim=c(nsites))
crops6km = array(siteCovs6km$Cultivated.Crops, dim= c(nsites))
forest6km = array(siteCovs6km$All.Forest, dim=c(nsites))
road6km = array(siteCovs6km$road_density, dim=c(nsites))
edge6km = array(siteCovs6km$forest_edgedensity, dim=c(nsites))
snowdepth6km = array(siteCovs6km$snow_depth, dim=c(nsites))
days = array(date, dim=c(nsites,nsurveys))
occ =  array(occ, dim=c(nsites,nsurveys))
corerange6km = array(siteCovs6km$core_range, dim=c(nsites))

#center all covariates
forest15km.c = standardize(forest15km)
conifmixed15km.c = standardize(conifmixed15km)
deer15km.c = standardize(deer15km)
deciduous15km.c = standardize(deciduous15km)
road15km.c =  standardize(road15km)
edge15km.c = standardize(edge15km)
days.c =  standardize(days)
crops15km.c = standardize(crops15km)
snowdepth15km.c = standardize(snowdepth15km)
wetlands15km.c = standardize(wetlands15km)

forest6km.c = standardize(forest6km)
conifmixed6km.c = standardize(conifmixed6km)
deciduous6km.c = standardize(deciduous6km)
road6km.c =  standardize(road6km)
edge6km.c = standardize(edge6km)
days.c =  standardize(days)
crops6km.c = standardize(crops6km)
snowdepth6km.c = standardize(snowdepth6km)
wetlands6km.c = standardize(wetlands6km)


#3D format

conifmixed15km.year <- cbind(conifmixed15km.c, yM)
conifmixed15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(conifmixed15km)){ #loop through each row
  this.year=which(years==conifmixed15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==conifmixed15km.year$sitename[i]) #get site for this row
  conifmixed15km.new[this.site,this.year]=as.numeric(conifmixed15km.year[i,1]) #force numeric
}
str(conifmixed15km.new)

deer15km.year <- cbind(deer15km.c, yM)
deer15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(deer15km)){ #loop through each row
  this.year=which(years==deer15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==deer15km.year$sitename[i]) #get site for this row
  deer15km.new[this.site,this.year]=as.numeric(deer15km.year[i,1]) #force numeric
}
str(deer15km.new)

deciduous15km.year <- cbind(deciduous15km.c, yM)
deciduous15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(deciduous15km)){ #loop through each row
  this.year=which(years==deciduous15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==deciduous15km.year$sitename[i]) #get site for this row
  deciduous15km.new[this.site,this.year]=as.numeric(deciduous15km.year[i,1]) #force numeric
}

forest15km.year <- cbind(forest15km.c, yM)
forest15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(forest15km)){ #loop through each row
  this.year=which(years==forest15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==forest15km.year$sitename[i]) #get site for this row
  forest15km.new[this.site,this.year]=as.numeric(forest15km.year[i,1]) #force numeric
}

snowdepth15km.year <- cbind(snowdepth15km.c, yM)
snowdepth15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(snowdepth15km)){ #loop through each row
  this.year=which(years==snowdepth15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==snowdepth15km.year$sitename[i]) #get site for this row
  snowdepth15km.new[this.site,this.year]=as.numeric(snowdepth15km.year[i,1]) #force numeric
}



road15km.year <- cbind(road15km.c, yM)
road15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(road15km)){ #loop through each row
  this.year=which(years==road15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==road15km.year$sitename[i]) #get site for this row
  road15km.new[this.site,this.year]=as.numeric(road15km.year[i,1]) #force numeric
}


edge15km.year <- cbind(edge15km.c, yM)
edge15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(edge15km)){ #loop through each row
  this.year=which(years==edge15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==edge15km.year$sitename[i]) #get site for this row
  edge15km.new[this.site,this.year]=as.numeric(edge15km.year[i,1]) #force numeric
}

wetlands15km.year <- cbind(wetlands15km.c, yM)
wetlands15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(wetlands15km)){ #loop through each row
  this.year=which(years==wetlands15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==wetlands15km.year$sitename[i]) #get site for this row
  wetlands15km.new[this.site,this.year]=as.numeric(wetlands15km.year[i,1]) #force numeric
}

crops15km.year <- cbind(crops15km.c, yM)
crops15km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(crops15km)){ #loop through each row
  this.year=which(years==crops15km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==crops15km.year$sitename[i]) #get site for this row
  crops15km.new[this.site,this.year]=as.numeric(crops15km.year[i,1]) #force numeric
}


conifmixed6km.year <- cbind(conifmixed6km.c, yM)
conifmixed6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(conifmixed6km)){ #loop through each row
  this.year=which(years==conifmixed6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==conifmixed6km.year$sitename[i]) #get site for this row
  conifmixed6km.new[this.site,this.year]=as.numeric(conifmixed6km.year[i,1]) #force numeric
}


deciduous6km.year <- cbind(deciduous6km.c, yM)
deciduous6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(deciduous6km)){ #loop through each row
  this.year=which(years==deciduous6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==deciduous6km.year$sitename[i]) #get site for this row
  deciduous6km.new[this.site,this.year]=as.numeric(deciduous6km.year[i,1]) #force numeric
}

forest6km.year <- cbind(forest6km.c, yM)
forest6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(forest6km)){ #loop through each row
  this.year=which(years==forest6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==forest6km.year$sitename[i]) #get site for this row
  forest6km.new[this.site,this.year]=as.numeric(forest6km.year[i,1]) #force numeric
}

snowdepth6km.year <- cbind(snowdepth6km.c, yM)
snowdepth6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(snowdepth6km)){ #loop through each row
  this.year=which(years==snowdepth6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==snowdepth6km.year$sitename[i]) #get site for this row
  snowdepth6km.new[this.site,this.year]=as.numeric(snowdepth6km.year[i,1]) #force numeric
}



road6km.year <- cbind(road6km.c, yM)
road6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(road6km)){ #loop through each row
  this.year=which(years==road6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==road6km.year$sitename[i]) #get site for this row
  road6km.new[this.site,this.year]=as.numeric(road6km.year[i,1]) #force numeric
}


edge6km.year <- cbind(edge6km.c, yM)
edge6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(edge6km)){ #loop through each row
  this.year=which(years==edge6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==edge6km.year$sitename[i]) #get site for this row
  edge6km.new[this.site,this.year]=as.numeric(edge6km.year[i,1]) #force numeric
}

wetlands6km.year <- cbind(wetlands6km.c, yM)
wetlands6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(wetlands6km)){ #loop through each row
  this.year=which(years==wetlands6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==wetlands6km.year$sitename[i]) #get site for this row
  wetlands6km.new[this.site,this.year]=as.numeric(wetlands6km.year[i,1]) #force numeric
}

crops6km.year <- cbind(crops6km.c, yM)
crops6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(crops6km)){ #loop through each row
  this.year=which(years==crops6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==crops6km.year$sitename[i]) #get site for this row
  crops6km.new[this.site,this.year]=as.numeric(crops6km.year[i,1]) #force numeric
}

corerange6km.year <- cbind(corerange6km, yM)
corerange6km.new=array(0,dim=c(maxJ,n.year)) #new data 
for(i in 1:nrow(corerange6km)){ #loop through each row
  this.year=which(years==corerange6km.year$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==corerange6km.year$sitename[i]) #get site for this row
  corerange6km.new[this.site,this.year]=as.numeric(corerange6km.year[i,1]) #force numeric
}

#reset K to number of occassions
K = 3

#set D to number of sampling days in an occassion
D = 7
# Ds=1


  
  n.chains = 3
  chains = vector("list", n.chains)
  chains2 = vector("list", n.chains)
  for(chain in 1:n.chains){
    source('C:/Users/jpt93/Documents/R/Occu-abu model/github/occupancy_abundance_model/case_study/coyote_fisher_marten_abu_abu_abu_binomialobsmodel_nimble_model.R')
  # Parameters monitored
  parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0I", "alphaI", "beta0I", "betaI", "alpha0S", "alphaS", "beta0S", "betaS", "gamma0DS","gamma0DI","gamma0IS") 
                   # "alpha0D.it", "alpha0I.it", "alpha0S.it", "beta0D.it", "beta0I.it", "beta0S.it")
                   # 'totalnD', 'totalZs', 'totalnI', 'T1obsD', 'T1simD', 'T1obsS', 'T1simS', 'T1obsI', 'T1simI')

#data and constants
Nimdata <- list(yD.new = yC.new, yI.new = yF.new, yS.new = yM.new)

# Nimdata <- list(yD.new = yF.new)

constants=list(J = J, nsurveys = dim(yF.new)[2],nyears = dim(yF.new)[3], corerange = corerange6km.new, deer = deer.new, beaver = beaver.new, moose = moose.new, misc = misc.new, temp = temp.new,  
               DATE = date.new, occ = occ.new, forestD = forest15km.new, deerD = deer15km.new, deciduousD = deciduous15km.new, edgeD = edge15km.new, agricultureD = crops15km.new, snowdepthD = snowdepth15km.new,
               conifmixedI = conifmixed15km.new, deciduousI = deciduous15km.new, agricultureI = crops15km.new, snowdepthI = snowdepth15km.new,
               conifmixedS = conifmixed6km.new, deciduousS = deciduous6km.new, snowdepthS = snowdepth6km.new, K2D = K2D.new, D = D)
#inits
NstNAD <- array(NA, dim=c(maxJ, K, n.year))
for(i in 1:nrow(yC)){ #loop through each row
  this.year=which(years==yC$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==yC$sitename[i]) #get site for this row
  NstNAD[this.site,1:3, this.year]=as.numeric(yC[i,2:4]) #force numeric
}
NstD <-  apply(NstNAD, c(1,3), max, na.rm = TRUE) # Inits for latent N

NstD[NstD == '-Inf'] <- NA

NstNAI <- array(NA, dim=c(maxJ, K, n.year))
for(i in 1:nrow(yF)){ #loop through each row
  this.year=which(years==yF$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==yF$sitename[i]) #get site for this row
  NstNAI[this.site,1:3, this.year]=as.numeric(yF[i,3:5]) #force numeric
}
NstI <-  apply(NstNAI, c(1,3), max, na.rm = TRUE) # Inits for latent N

NstI[NstI == '-Inf'] <- NA


NstNAS <- array(NA, dim=c(maxJ, K, n.year))
for(i in 1:nrow(yM)){ #loop through each row
  this.year=which(years==yM$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==yM$sitename[i]) #get site for this row
  NstNAS[this.site,1:3, this.year]=as.numeric(yM[i,2:4]) #force numeric
}
NstS <-  apply(NstNAS, c(1,3), max, na.rm = TRUE) # Inits for latent N

NstS[NstS == '-Inf'] <- NA



# # if using random effects uncomment out below (but we do not have sufficent data to support yearly random effects) 
# # #random year level effect on psi inits
# beta0S.it  <- rnorm(J)
# beta0D.it  <- rnorm(J)
# beta0I.it  <- rnorm(J)
# 
# #random year level effect on p inits
# alpha0S.it <- rnorm(J)
# alpha0D.it <- rnorm(J)
# alpha0I.it <- rnorm(J)



# #random site level effect on p inits
# eps <- rnorm(totalJ)
# eps.year <- cbind(eps, yF)
# eps.new = array(0,dim=c(maxJ,n.year)) #new data 
# for(i in 1:nrow(yF)){ #loop through each row
#   this.year=which(years==eps.year$year[i]) #get year for this row
#   this.site=which(year.sites[[this.year]]==eps.year$sitename[i]) #get site for this row
#   eps.new[this.site,this.year]=as.numeric(eps.year[i,1]) #force numeric
# }


#inits list

Niminits=list(beta0D=c(0,0,0),beta0S=c(0,0,0),beta0I=c(0,0,0),nD=NstD, nI = NstI, nS=NstS,alpha0S=c(0,0,0), alpha0I=c(0,0,0), alpha0D=c(0,0,0), alphaD=c(0,0,0), alphaI=c(0,0,0,0),alphaS=c(0,0,0,0), betaD=c(0,0),betaI=c(0,0,0), betaS=c(0,0,0), gamma0DS=0, gamma0DI=0, gamma0IS=0)

#thinning rate
nt = 8
nt2 = 200
# 
# nt = 2
# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, monitors2 = parameters, thin=nt, thin2 = nt2,useConjugacy = TRUE, enableWAIC=TRUE)

conf$addSampler(target = c("alpha0D[1]","beta0D[1]", "beta0S[1]"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("alpha0D[2]","beta0D[2]"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("alpha0D[3]","beta0D[3]"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# # #quick test iter
#   n.iter <- 30000

#full run iter
n.iter <- 200000
# Run the model.
start.time2<-Sys.time()
Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com

#get the chains

mvSamples = as.matrix(Cmcmc$mvSamples)
mvSamples2 = as.matrix(Cmcmc$mvSamples2)

chains[[chain]]=mvSamples
chains2[[chain]]=mvSamples2

}

#combine the chains and burn
n.iter=25000
n.burn = 20000
  # plot(mcmc(mvSamples[n.burn:nrow(mvSamples),]))
library(coda)  

# colnames(mvSamples)
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

plot(a)

c=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

plot(c)
gelman <- gelman.diag(a)

# save(chains, file='Coyote-fisher-marten-abu-abu-abu-model-200k-chains.RData')
# save(chains2, file='Coyote-fisher-marten-abu-abu-abu-model-200k-chains2.RData')

load('C:/Users/jpt93/Documents/R/Occu-abu model/case_study/Coyote-fisher-marten-abu-abu-abu-model-200k-chains.RData')

write.csv(gelman$psrf, 'coyote-fisher-marten-abu-abu-abu-200k-gelmanrubin.csv')

a=runjags::combine.mcmc(a)


cor <- cor(a[, c("alpha0D[1]", "alpha0D[2]", "alpha0D[3]", "alpha0I[1]", "alpha0I[2]", "alpha0I[3]", "alpha0S[1]", "alpha0S[2]", "alpha0S[3]", "alphaD[1]", "alphaD[2]",  
                 "alphaD[3]",  "alphaI[1]",  "alphaI[2]","alphaI[3]",  "alphaI[4]",  "alphaS[1]",  "alphaS[2]",  "alphaS[3]", "alphaS[4]", "beta0D[1]",  "beta0D[2]", 
                 "beta0D[3]",  "beta0I[1]", "beta0I[2]",  "beta0I[3]", "beta0S[1]",  "beta0S[2]",  "beta0S[3]", "betaD[1]",   "betaD[2]",   "betaI[1]", "betaI[2]",  
                 "betaI[3]", "betaS[1]",   "betaS[2]",   "betaS[3]",   "gamma0DI",   "gamma0DS","gamma0IS")])

summary(a)

# load('Coyote-fisher-marten-abu-abu-abu-model-100k-chains.RData')
#extract point estimates and credible intervals
esti <- array(NA, dim = c(3, (dim(a)[2])))
colnames(a)

esti[1,1] <- plogis(mean(a[,"alpha0D[1]"]))
esti[c(2,3),1] <- plogis(quantile(a[,"alpha0D[1]"], probs = c(2.5, 97.5)/100))

esti[1,2] <- plogis(mean(a[,"alpha0D[2]"]))
esti[c(2,3),2] <- plogis(quantile(a[,"alpha0D[2]"], probs = c(2.5, 97.5)/100))

esti[1,3] <- plogis(mean(a[,"alpha0D[3]"]))
esti[c(2,3),3] <- plogis(quantile(a[,"alpha0D[3]"], probs = c(2.5, 97.5)/100))

esti[1,4] <- plogis(mean(a[,"alpha0I[1]"]))
esti[c(2,3),4]<- plogis(quantile(a[,"alpha0I[1]"], probs = c(2.5, 97.5)/100))

esti[1,5] <- plogis(mean(a[,"alpha0I[2]"]))
esti[c(2,3),5]<- plogis(quantile(a[,"alpha0I[2]"], probs = c(2.5, 97.5)/100))

esti[1,6] <- plogis(mean(a[,"alpha0I[3]"]))
esti[c(2,3),6]<- plogis(quantile(a[,"alpha0I[3]"], probs = c(2.5, 97.5)/100))

esti[1,7] <- plogis(mean(a[,"alpha0S[1]"]))
esti[c(2,3),7]<- plogis(quantile(a[,"alpha0S[1]"], probs = c(2.5, 97.5)/100))

esti[1,8] <- plogis(mean(a[,"alpha0S[2]"]))
esti[c(2,3),8]<- plogis(quantile(a[,"alpha0S[2]"], probs = c(2.5, 97.5)/100))

esti[1,9] <- plogis(mean(a[,"alpha0S[3]"]))
esti[c(2,3),9]<- plogis(quantile(a[,"alpha0S[3]"], probs = c(2.5, 97.5)/100))

esti[1,10] <- mean(a[,"alphaD[1]"])
esti[c(2,3),10]<- quantile(a[,"alphaD[1]"], probs = c(2.5, 97.5)/100)

esti[1,11] <- mean(a[,"alphaD[2]"])
esti[c(2,3),11]<- quantile(a[,"alphaD[2]"], probs = c(2.5, 97.5)/100)

esti[1,12] <- mean(a[,"alphaD[3]"])
esti[c(2,3),12]<- quantile(a[,"alphaD[3]"], probs = c(2.5, 97.5)/100)

esti[1,13] <- mean(a[,"alphaI[1]"])
esti[c(2,3),13]<- quantile(a[,"alphaI[1]"], probs = c(2.5, 97.5)/100)

esti[1,14] <- mean(a[,'alphaI[2]'])
esti[c(2,3),14]<- quantile(a[,"alphaI[2]"], probs = c(2.5, 97.5)/100)

esti[1,15] <- mean(a[,'alphaI[3]'])
esti[c(2,3),15]<- quantile(a[,"alphaI[3]"], probs = c(2.5, 97.5)/100)

esti[1,16] <- mean(a[,'alphaI[4]'])
esti[c(2,3),16]<- quantile(a[,"alphaI[4]"], probs = c(2.5, 97.5)/100)

esti[1,17] <- mean(a[,"alphaS[1]"])
esti[c(2,3),17]<- quantile(a[,"alphaS[1]"], probs = c(2.5, 97.5)/100)

esti[1,18] <- mean(a[,"alphaS[2]"])
esti[c(2,3),18]<- quantile(a[,"alphaS[2]"], probs = c(2.5, 97.5)/100)

esti[1,19] <- mean(a[,"alphaS[3]"])
esti[c(2,3),19]<- quantile(a[,"alphaS[3]"], probs = c(2.5, 97.5)/100)

esti[1,20] <- mean(a[,"alphaS[4]"])
esti[c(2,3),20]<- quantile(a[,"alphaS[4]"], probs = c(2.5, 97.5)/100)

esti[1,21] <- mean(a[,"beta0D[1]"])
esti[c(2,3),21]<- quantile(a[,"beta0D[1]"], probs = c(2.5, 97.5)/100)

esti[1,22] <- mean(a[,"beta0D[2]"])
esti[c(2,3),22]<- quantile(a[,"beta0D[2]"], probs = c(2.5, 97.5)/100)

esti[1,23] <- mean(a[,"beta0D[3]"])
esti[c(2,3),23]<- quantile(a[,"beta0D[3]"], probs = c(2.5, 97.5)/100)

esti[1,24] <- mean(a[,"beta0I[1]"])
esti[c(2,3),24]<- quantile(a[,"beta0I[1]"], probs = c(2.5, 97.5)/100)

esti[1,25] <- mean(a[,"beta0I[2]"])
esti[c(2,3),25]<- quantile(a[,"beta0I[2]"], probs = c(2.5, 97.5)/100)

esti[1,26] <- mean(a[,"beta0I[3]"])
esti[c(2,3),26]<- quantile(a[,"beta0I[3]"], probs = c(2.5, 97.5)/100)

esti[1,27] <- mean(a[,"beta0S[1]"])
esti[c(2,3),27]<- quantile(a[,"beta0S[1]"], probs = c(2.5, 97.5)/100)

esti[1,28] <- mean(a[,"beta0S[2]"])
esti[c(2,3),28]<- quantile(a[,"beta0S[2]"], probs = c(2.5, 97.5)/100)

esti[1,29] <- mean(a[,"beta0S[3]"])
esti[c(2,3),29]<- quantile(a[,"beta0S[3]"], probs = c(2.5, 97.5)/100)

esti[1,30] <- mean(a[,"betaD[1]"])
esti[c(2,3),30]<- quantile(a[,"betaD[1]"], probs = c(2.5, 97.5)/100)

esti[1,31] <- mean(a[,"betaD[2]"])
esti[c(2,3),31]<- quantile(a[,"betaD[2]"], probs = c(2.5, 97.5)/100)

esti[1,32] <- mean(a[,"betaI[1]"])
esti[c(2,3),32]<- quantile(a[,"betaI[1]"], probs = c(2.5, 97.5)/100)

esti[1,33] <- mean(a[,"betaI[2]"])
esti[c(2,3),33]<- quantile(a[,"betaI[2]"], probs = c(2.5, 97.5)/100)

esti[1,34] <- mean(a[,"betaI[3]"])
esti[c(2,3),34]<- quantile(a[,"betaI[3]"], probs = c(2.5, 97.5)/100)

esti[1,35] <- mean(a[,"betaS[1]"])
esti[c(2,3),35]<- quantile(a[,"betaS[1]"], probs = c(2.5, 97.5)/100)

esti[1,36] <- mean(a[,"betaS[2]"])
esti[c(2,3),36]<- quantile(a[,"betaS[2]"], probs = c(2.5, 97.5)/100)

esti[1,37] <- mean(a[,"betaS[3]"])
esti[c(2,3),37]<- quantile(a[,"betaS[3]"], probs = c(2.5, 97.5)/100)

esti[1,38] <- mean(a[,"gamma0DI"])
esti[c(2,3),38]<- quantile(a[,"gamma0DI"], probs = c(2.5, 97.5)/100)

esti[1,39] <- mean(a[,"gamma0DS"])
esti[c(2,3),39]<- quantile(a[,"gamma0DS"], probs = c(2.5, 97.5)/100)

esti[1,40] <- mean(a[,"gamma0IS"])
esti[c(2,3),40]<- quantile(a[,"gamma0IS"], probs = c(2.5, 97.5)/100)

colnames(a)
# [1] "alpha0D[1]" "alpha0D[2]" "alpha0D[3]" "alpha0I[1]" "alpha0I[2]" "alpha0I[3]" "alpha0S[1]" "alpha0S[2]" "alpha0S[3]" "alphaD[1]"  "alphaD[2]" 
# [12] "alphaD[3]"  "alphaI[1]"  "alphaI[2]"  "alphaI[3]"  "alphaI[4]"  "alphaS[1]"  "alphaS[2]"  "alphaS[3]"  "alphaS[4]"  "beta0D[1]"  "beta0D[2]" 
# [23] "beta0D[3]"  "beta0I[1]"  "beta0I[2]"  "beta0I[3]"  "beta0S[1]"  "beta0S[2]"  "beta0S[3]"  "betaD[1]"   "betaD[2]"   "betaI[1]"   "betaI[2]"  
# [34] "betaI[3]"   "betaS[1]"   "betaS[2]"   "betaS[3]"   "gamma0DI"   "gamma0DS"   "gamma0IS" 

colnames(esti) <- c("alpha0D[1]", "alpha0D[2]", "alpha0D[3]", "alpha0I[1]", "alpha0I[2]", "alpha0I[3]", "alpha0S[1]", "alpha0S[2]", "alpha0S[3]", "alphaD[1]", "alphaD[2]",  
                    "alphaD[3]",  "alphaI[1]",  "alphaI[2]","alphaI[3]",  "alphaI[4]",  "alphaS[1]",  "alphaS[2]",  "alphaS[3]", "alphaS[4]", "beta0D[1]",  "beta0D[2]", 
                    "beta0D[3]",  "beta0I[1]", "beta0I[2]",  "beta0I[3]", "beta0S[1]",  "beta0S[2]",  "beta0S[3]", "betaD[1]",   "betaD[2]",   "betaI[1]", "betaI[2]",  
                    "betaI[3]", "betaS[1]",   "betaS[2]",   "betaS[3]",   "gamma0DI",   "gamma0DS","gamma0IS")

write.csv(esti, 'coyote-fisher-marten-abuabu-posterior-coefs-200k.csv')


# plot these badboys
setwd('C:/Users/jpt93/Documents/R/Occu-abu model/')
coefs <- read.csv('C:/Users/jpt93/Documents/R/Occu-abu model/case_study/coyote-fisher-marten-abuabu-200k-posteriorcoefs-formatted.csv')


#coefs

#colour blind friendly pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# coefs$Parameter <- factor(coefs$Parameter, levels = c("snow depth", "forest edge", "deer availability", "deciduous forest", "conifer-mixed forest", "all forest", "fisher-marten interaction", "coyote-fisher interaction", "coyote-marten interaction"))

#coefs

#colour blind friendly pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


marten <- subset(coefs, Species == "marten")
marten$Parameter <- factor(marten$Parameter, levels = c("fisher abundance", "coyote abundance", "snow depth", "deciduous forest", "conifer-mixed forest", "bait", "linear date", "quadratic date", "occasion"))

fisher <- subset(coefs, Species == "fisher")
fisher$Parameter <- factor(fisher$Parameter, levels = c("coyote abundance", "snow depth", "deciduous forest", "conifer-mixed forest", "minimum temperature", "linear date", "quadratic date", "occasion"))

coyote <- subset(coefs, Species == "coyote")
coyote$Parameter <- factor(coyote$Parameter, levels = c("deer availability", "forest edge", "linear date", "quadratic date", "occasion"))
statecoefs <- subset(coefs, Submodel == "state")
observationcoefs <- subset(coefs, Submodel == "observation")

tiff("marten_abuabuabu_coefs.tiff", units="in", width=4, height=3.5, res=400)

ggplot(data = marten, aes(x = posterior_mode, y = Parameter, group = Submodel, colour = Submodel, fill = Submodel))+geom_point(aes(), size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+ guides(colour = "none", fill= "none")+
  ylab("parameter") + xlab("coefficient estimates") + xlim(-1,2)  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

dev.off()

tiff("coyote_abuabuabu_coefs.tiff", units="in", width=4, height=3.5, res=400)

ggplot(data = coyote, aes(x = posterior_mode, y = Parameter, group = Submodel, colour = Submodel, fill = Submodel))+geom_point(aes(), size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+ guides(colour = "none", fill= "none")+
  ylab("parameter") + xlab("coefficient estimates") + xlim(-0.5,1)  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

dev.off()

tiff("fisher_abuabuabu_coefs.tiff", units="in", width=4, height=3.5, res=400)

ggplot(data = fisher, aes(x = posterior_mode, y = Parameter, group = Submodel, colour = Submodel, fill = Submodel))+geom_point(aes(), size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+ guides(colour = "none", fill= "none")+
  ylab("parameter") + xlab("coefficient estimates") + xlim(-0.5,1)  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

dev.off()

coefs$overlap<- as.factor(coefs$overlap)
coefs$Submodel<- as.factor(coefs$Submodel)



marten <- subset(coefs, Species == "marten")
marten$Parameter <- factor(marten$Parameter, levels = c("fisher abundance", "coyote abundance", "snow depth", "deciduous forest", "conifer-mixed forest", "linear date", "quadratic date", "occasion"))

fisher <- subset(coefs, Species == "fisher")
fisher$Parameter <- factor(fisher$Parameter, levels = c("coyote abundance", "snow depth", "deciduous forest", "conifer-mixed forest", "linear date", "quadratic date", "occasion"))

coyote <- subset(coefs, Species == "coyote")
coyote$Parameter <- factor(coyote$Parameter, levels = c("snow depth", "all forest", "deer availability", "forest edge", "linear date", "quadratic date", "occasion"))
statecoefs <- subset(coefs, Submodel == "state")
observationcoefs <- subset(coefs, Submodel == "observation")

ggplot(data = marten, aes(x = posterior_mode, y = Parameter, group = Submodel, colour = Submodel, fill = Submodel))+geom_point(aes(), size = 2.5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+
  ylab("parameter") + xlab("posterior coefficient estimates") + xlim(-1,2)  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid


ggplot(data = coyote, aes(x = posterior_mode, y = Parameter, group = Submodel, colour = Submodel, fill = Submodel))+geom_point(aes(), size = 2.5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+
  ylab("parameter") + xlab("posterior coefficient estimates") + xlim(-0.5,1)  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

ggplot(data = fisher, aes(x = posterior_mode, y = Parameter, group = Submodel, colour = Submodel, fill = Submodel))+geom_point(aes(), size = 2.5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+
  ylab("parameter") + xlab("posterior coefficient estimates (log-odds scale)") + xlim(-0.5,1)  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

