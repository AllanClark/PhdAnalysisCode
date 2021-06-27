# Author: Allan Clark
# Date 1 January 2021
# The following code fits an Dorazio Royle 2006 species richness model to
# data collected by Marine Drouilly
# The area included here is the Koup
# Denote psi[i,j] as occupancy prob for species i at site j
# p[i,j] is similar
# logit(psi[i,j]) = mean.beta + u[i] + sum(k; beta[i,k]*x[i,k])
# mean.beta denotes overall species level intercept
# u[i] = species specific random intercept
# no site specific random effect is included in the occupancy process
# logit(p[i,j]) = mean.alpha + v[i] + sum(k; alpha[i,k]*w[i,k])
# mean.alpha denotes overall species level intercept
# v[i] = species specific random intercept
# no site specific random effect is included in the detection process
# The variables included in the analysis are as follows:
# occupancy - plain (binary), river indicator (binary), elevation (cont),
# MSAVI2 (cont), Disturbance index (cont), SmallPreyindex (cont),
# PreyIndex (cont), LiveStockIndex (cont)
# detection - trail, plain and river indicator variables
#-------------------------------------------------------------------------------

#Load the libraries

library(reshape)
require(stocc)
require(RColorBrewer)
require(plotly)
require(fields)
require(alphahull)
require(scales)
require(MSO)

options(max.print=1000000000)

#-------------------------------------------------------------------------------

#Functions to perform Multi-species occupancy model
#Here we assume that we know the species richness
#and investigate the regression effects related to variouis covariates

MSOBinlogit<-function(detection.model, occupancy.model,
                      occ.covs, det.covs, Y,
                      prior = NULL,
                      control = NULL,
                      nsitevisits,
                      selection = 1){

  #-----------------------------------------------------------------------------
  #Fits a nonspatial model to Binomial data
  #The assumption is that species richness is known
  #selection =1 ==> return only model selection stats
  #-----------------------------------------------------------------------------

  X <- as.matrix(model.matrix(occupancy.model, occ.covs)) #occupancy design matrix
  V <- as.matrix( model.matrix(detection.model, det.covs) ) #detection design matrix

  ns <- NROW(Y) #number of species
  J <- NCOL(Y) #number of sites

  #The priors used
  a2 <- prior$a2 #residual variance of detection regression effects
  b2 <- prior$b2 #residual variance of occupancy regression effects
  A2 <- prior$A2 #for detection regression effects
  B2 <- prior$B2 #for occupancy regression effects

  #the control elements
  ndraws <- control$ndraws
  percent_burn_in <- control$percent_burn_in
  thin <- control$thin

  z <- t(Y>0) +0 #number of sites by number of species
  siteids <- matrix(1:J, ncol=1)

  p <- z #detection probs
  psi <- p #occupancy probs

  #go to Rcpp and fit the models
  MSOBinlogitcpp(X, V, Y, z, p, psi,
                 nsitevisits,
                 a2, b2, A2, B2,
                 ndraws, percent_burn_in, thin,
                 selection)
}

MSOBinspatlogit<-function(detection.model, occupancy.model, spatial.model,
                          occ.covs, det.covs, Y,
                          site.chars,
                          prior,
                          control){

  #-----------------------------------------------------------------------------
  #Fits a spatial model to Binomial data
  #The assumption is that species richness is known
  #The spatial model is accounted for by using a restricted spatial regression
  #model
  #The arguments used are:
  #detection.model
  #occupancy.model
  #spatial.model
  #occ.covs
  #det.covs
  #Y
  #site.chars
  #prior = NULL,
  #control
  #-----------------------------------------------------------------------------

  #design matrices
  X <- as.matrix(model.matrix(occupancy.model, occ.covs)) #occupancy design matrix
  V <- as.matrix( model.matrix(detection.model, det.covs) ) #detection design matrix

  #some constants
  ns <- NROW(Y) #number of species
  J <- NCOL(Y) #number of sites

  #the control elements
  ndraws <- control$ndraws
  percent_burn_in <- control$percent_burn_in
  thin <- control$thin

  #the site characteristics
  nsitevisits <- site.chars$nsitevisits
  xy <- site.chars$xycords

  #The priors used
  a2 <- prior$a2 #residual variance of detection regression effects
  b2 <- prior$b2 #residual variance of occupancy regression effects
  A2 <- prior$A2 #for detection regression effects
  B2 <- prior$B2 #for occupancy regression effects
  a.tau_i <- prior$a.tau_i
  b.tau_i <- prior$b.tau_i

  #z matrix and probs
  z <- t(Y>0) +0 #number of sites by number of species
  siteids <- matrix(1:J, ncol=1)
  p <- z #detection probs
  psi <- p #occupancy probs

  #-----------------------------------------------------------------------------
  #the spatial bits required
  Q <- stocc::icar.Q(xy = xy, threshold = spatial.model$threshold, rho = 1)
  A <- diag(diag(Q)) - Q
  n <- NCOL(Q)
  P <- diag(n) - X %*% solve(crossprod(X), t(X))
  Op <- (nrow(A)/sum(A)) * (P %*% (A %*% P))
  numSpatre <- spatial.model$moran.cut
  e <- rARPACK::eigs(Op, numSpatre)
  K <- e$vectors[, 1:numSpatre]
  Minv <- crossprod(K, Q) %*% K
  #-----------------------------------------------------------------------------

  #starting values for the spatial tau parameters
  tau_i <- matrix(1:ns, nrow=ns)
  #-----------------------------------------------------------------------------

  #go to Rcpp and fit the models
  MSOBinocclogitcpp(X, V, Y, z,
                    K, Minv,
                    p, psi,
                    tau_i, a.tau_i, b.tau_i,
                    nsitevisits,
                    a2, b2, A2, B2,
                    ndraws, percent_burn_in, thin)
}

#-------------------------------------------------------------------------------

make_coords2 <- function(xy){
  #this assumes that the centroids were calculated by using 2.5/60
  x<-xy[1]
  y<-xy[2]
  #shift<- -0.1533488/2 #half a minute

  up <- -0.2434459/2
  hor <- -0.1533488/2
  a1<-c(x-hor,y+up)
  a2<-c(x+hor, y+up)
  a3<-c(x+hor, y-up)
  a4<-c(x-hor, y-up)
  #a5<-a1
  df<-as.data.frame(rbind(a1,a2,a3,a4))
  df<-cbind(0, 1:4, df)
  colnames(df)<-c("PID","POS","X","Y")
  return(df)
}

image.scale2 <- function(z, zlim, col = heat.colors(12),
                         breaks, horiz=TRUE, ylim=NULL, xlim=NULL, side=2, cex.axis=1, ...)
{
  #Taken from "me nugget" Add in reference!
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directio
    ns
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  #yaxt <- ifelse(horiz, "n", "s")
  yaxt <- "n"
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,8)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)
  #axis(side=side, cex.axis=cex.axis, tick=FALSE)
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }

  for (i in 2:5){text(x=1.5, y=(i-1)*.2, labels=(i-1)*.2, cex=.8)}
  return(col)
}

make_map <- function(map_variable, ngroups_colours, xycords){
  #ngroups_colours <- 300
  #breaks <- seq(0, 1,length.out=ngroups_colours)


  cbPalette <- brewer.pal(n=9, name="YlOrRd")
  map_variable <- c(0, map_variable, 1)
  Colours<-colorRampPalette(cbPalette)(ngroups_colours)[as.numeric(cut(map_variable,ngroups_colours))]

  plot(longitude, latitude, type="n", xaxt="n", yaxt="n", bty="n")
  box()

  for (i in 1:nrow(xycords)){
    xx<-make_coords2( c(xycords[i,]) )
    polygon(xx[,3], xx[,4], col= Colours[i], border=NA)
    #polygon(xx[,3], xx[,4], col= "red", border=NULL)
  }

  # #d <- data.frame(longitude, latitude, z)
  # #qplot(longitude, latitude,colour=z,data=d)
  #
  for (i in 1:(length(bdrpts)-1)){
    segments(x[bdrpts[i]], y[bdrpts[i]], x[bdrpts[i+1]], y[bdrpts[i+1]], lty=3)
  }
  segments(x[bdrpts[i+1]], y[bdrpts[i+1]], x[bdrpts[1]], y[bdrpts[1]], lty=3)
}

make_map2 <- function(map_variable, ngroups_colours, xycords){

  map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  }
  cbPalette <- viridis_pal(direction=-1)(ngroups_colours)

  Colours<-map2color(map_variable, cbPalette,limits=c(0,1))
  #Colours<-colorRampPalette(cbPalette)(ngroups_colours)[as.numeric(cut(map_variable,ngroups_colours))]

  plot(longitude, latitude, type="n", xaxt="n", yaxt="n", bty="n", asp=1)
  box()

  for (i in 1:nrow(xycords)){
    xx<-make_coords2( c(xycords[i,]) )
    polygon(xx[,3], xx[,4], col= Colours[i], border=NA)
    #polygon(xx[,3], xx[,4], col= "red", border=NULL)
  }

  # #d <- data.frame(longitude, latitude, z)
  # #qplot(longitude, latitude,colour=z,data=d)
  #
  for (i in 1:(length(bdrpts)-1)){
    segments(x[bdrpts[i]], y[bdrpts[i]], x[bdrpts[i+1]], y[bdrpts[i+1]], lty=3)
  }
  segments(x[bdrpts[i+1]], y[bdrpts[i+1]], x[bdrpts[1]], y[bdrpts[1]], lty=3)
}

#-------------------------------------------------------------------------------

require(devtools)
#install_github("AllanClark/MSO", auth_token = "cf29ca0d79a4153a11524a4e79e780e910c18132", force=TRUE)

#set the working directory and some file names
office <- FALSE
laptop <- FALSE
Harddrive <- FALSE
OneDrive <- TRUE

if (office==TRUE){
  directory1<-"~/Research/PhD/Paper4"
  setwd(directory1)
}

if (laptop ==TRUE){
  directory1<-"/Users/allanclark/Documents/AllanClark/Stats/Research/Clark_AE/PhD/Paper4/AnalysisExamples/Example1/"
  setwd(paste(directory1,"/RequiredData", sep=""))
}

if (Harddrive ==TRUE){
  directory1<-"E:/Backup7Dec2020/Documents/AllanClark/Stats/Research/Clark_AE/PhD/Paper4/AnalysisExamples/Example1"
  setwd(paste(directory1,"/RequiredData", sep=""))
}

if (OneDrive ==TRUE){
  #Using my OneDrive account at UCT
  directory1<-"C:/Users/User/OneDrive - University of Cape Town/Research/PhDChapt4/Paper4_OfficeComputer Working"
  setwd(directory1)
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#1. Begin 'data wrangling'
#-------------------------------------------------------------------------------

#load the data
#define the data files
occfile<-"Pool6.Data.csv"
stationfile<-"station_data_karoo_final.csv"

#Load the observation data
data <- read.table(occfile, header=TRUE,sep=",",na.strings=TRUE)
names(data)
#-------------------------------------------------------------------------------

#remove some of the stations!
#these stations don't have any pictures!
remove.stations<-c("8","10","15" , "44" , "47" , "48" , "56" , "57" , "119", "127")
for (i in 1:length(remove.stations))
{
    data <- subset(data, Station != remove.stations[i] )
}

remove.species<-c("rodent" )
for (i in 1:length(remove.species))
{
    data <- subset(data, Species != remove.species[i] )
}

length(unique(data$Station)) #322
length(unique(data$Species)) #42
#-------------------------------------------------------------------------------

#List the unique species
u.species = as.character(sort(unique(data$Species)))
cat("\n The names of the species that is kept.")
u.species

#Number of observed species
n=length(u.species)
n

#List the  camera stations
ustations = as.character(sort(unique(data$Station)))
ustations
length(ustations)#should be 322

#Number of camera stations
J=length(ustations)
J
#-------------------------------------------------------------------------------

#Convert the data to a matrix format with the first dimension the stations
#and the second dimension the species.
#ie station number by species type
#records the number days out of the total days that the station was in use
#that a species was recorded
melt.tmp=melt(data,id.var=c("Species", "Station"), measure.var="Count")
Y=cast(melt.tmp, Station ~ Species,sum)
Yorig<-Y
Yorig[,1]
cat("\n Y contains the number of times a species was seen at a camera station.")
cat("\n camera station number down the rows; species name along the column.")
Y<-Y[,-1]
dim(Y) #322 by 42

#what are the rownames of Y?
colnames(Y) #put into alphabetical order
rownames(Y)<-sort(unique(data$Station))
rownames(Y)

#so checking that the order is the same!
cbind(colnames(Y), u.species )

#Load the station covariates

station.cov <- read.csv(stationfile, header=TRUE,sep=";",na.strings=TRUE,check.names=T,
colClasses=c(rep("numeric",2),rep("factor",3),rep("numeric",2),rep("factor",3),
"numeric","factor", rep("numeric",8)) )
station.cov2<-station.cov
names(station.cov)
head(station.cov)
length(unique(station.cov$Station)) #322
length(unique(data$Station)) #322

#Number of days each camera was operating
K=round(station.cov$Nights/6,0)
K
#-------------------------------------------------------------------------------

#Check that the data is correct
minY<-apply(Y,1,min)
minY0<-which(minY<0)
if(length(minY0)>0) cat(c("The following cases should be checked since min
number of counts is too small", minY))

maxY<-apply(Y,1,max)
maxYK<-which((maxY>K) ==1)
if(length(maxYK)>0) cat(c("The following cases should be checked since max
number of counts is too large:", maxYK))

YK_check=cbind(K,apply(Y, 1, max))
YgreaterK<-which( (YK_check[,1]>=YK_check[,2])==0)

if (length(YgreaterK)>0)
{
    Ynew<-Y

    for (i in 1:length(YgreaterK))
    {
        icheck<-which(Ynew[ YgreaterK[i], ]> K[YgreaterK[i]])
        Ynew[YgreaterK[i], icheck ] <- K[YgreaterK[i]]
    }

    YnewK_check=cbind(K,apply(Ynew, 1, max))
    YnewgreaterK<-which( (YnewK_check[,1]>=YnewK_check[,2])==0)
}

#most cameras have very little shots!
apply(Y, 1, sum)

#number of camera shots per species
apply(Y,2,sum)
#-------------------------------------------------------------------------------

#I don't know how to change the latitude and longitude from station.cov
#to numeric values so I am reading them in again after making a new one!!!!
#Some plots regarding the spatial distribution of the detections

#station.cov2 <- read.csv("station_data_karoo_after_deleted_cameras2.csv", header=TRUE,sep=";",na.strings=TRUE,check.names=T,
#                         colClasses=c(rep("numeric",2),rep("factor",3),rep("numeric",2),rep("factor",3),
#                                      "numeric","factor", rep("numeric",8)) )

#Now for each row if the number of observed animals is equal to
#the number of Nights set it equal to the number of nights

#Survey covariate
survey<- as.factor(as.vector(station.cov$SurveyName))
survey.factor<- as.numeric(survey)
survey.levels<-length(levels(as.factor(survey)))

#The different Survey sites
u.Survey<-sort(unique(data$Survey))
u.Station<-sort(unique(data$Station))

#Habitat covariate
habitat<- as.factor(as.vector(station.cov$Habitat))
habitat.factor<- as.numeric(habitat)
habitat.levels<- length(levels(as.factor(habitat)))

#Biome covariate
biome<- as.factor(as.vector(station.cov$Biome))
biome.factor<- as.numeric(biome)
biome.levels<- length(levels(as.factor(biome)))

#Trail covariate (trail/off-trail_)
trail<-station.cov$Trail
trail<- as.factor(as.character(trail))
trail.factor<- as.numeric(trail)
trail.levels<- length(levels(as.factor(trail)))
Trail.ind<-trail.factor
Trail.ind[Trail.ind==2]<-0

DisturbanceIndex<-c(station.cov$DisturbanceIndex)
LivestockIndex<-c(station.cov$LivestockIndex)
MSAVI2<-c(station.cov$MSAVI2)

#Elevation
elevation<-c(station.cov$Elevation)

#identify the camera stations associated with 'Koup' and "Anysberg"
#these are the actual camera station numbers
Koup.Stations<-sort(unique(data[data$Survey=="Koup",][,2]))
Anysberg.Stations<-sort(unique(data[data$Survey=="Anysberg",][,2]))

#now obtain the row numbers associated with the station numbers!
Koup.rows<-which((rownames(Y) %in% Koup.Stations) == 1)
Anysberg.rows<-which((rownames(Y) %in% Anysberg.Stations) == 1)
#-------------------------------------------------------------------------------

#The two data matrices that will be used to undertake the Dorazio
#analysis in order to obtain the species richness at each of the sites
#remove the species from the data that does not have any counts!!
Koup.Data<-t(Y[Koup.rows,])
dim(Koup.Data) #42 by 166
remove.row.Koup<-which(apply(Koup.Data,1, function(x) sum(x>0))==0)
names(Y)[remove.row.Koup] #identify which species was not seen in Koup
Koup.Data<-Koup.Data[-remove.row.Koup,]
Koup.names<-dimnames(Y)[[2]]
Koup.names<-Koup.names[-remove.row.Koup]
length(Koup.names) #so 34 species is seen at Koup!
Koup.names

Anysberg.Data<-t(Y[Anysberg.rows,])
dim(Anysberg.Data) #42 by 156
remove.row.Anysberg<-which(apply(Anysberg.Data,1, function(x) sum(x>0))==0)
names(Y)[remove.row.Anysberg] #identify which species was not seen in Anysberg
#Anysberg.Data<-Anysberg.Data[-remove.row.Anysberg,]
Anysberg.Data<-Anysberg.Data[-remove.row.Anysberg, ]
Anysberg.names<-dimnames(Y)[[2]]
Anysberg.names<-Anysberg.names[-remove.row.Anysberg]
length(Anysberg.names) #so 35 species is seen at Anysberg!
Anysberg.names

#convert the problem to occupancy instead of counts!
Koup.Data2<-(Koup.Data>0)+0
Anysberg.Data2<-(Anysberg.Data>0)+0

#-------------------------------------------------------------------------------

#Habitat covariate
#-----------------
#here I will group kloof and mountain together
habitat2<- as.factor(as.vector(station.cov2$Habitat))
table(habitat2)
habitat2[habitat2=="kloof"]<-"mountain"
table(habitat2)
habitat2<-droplevels(habitat2)
habitat2.factor<- as.numeric(habitat2)
habitat2.levels<- length(levels(as.factor(habitat2)))
#indicator variables for the the different habitat covariates
mountain.index<-rep(0,J); mountain.index[habitat2=="mountain"]<-1
plain.index<-rep(0,J); plain.index[habitat2=="plain"]<-1
river.index<-rep(0,J); river.index[habitat2=="riverine"]<-1

#-------------------------------------------------------------------------------
index1<-167:322 #Anysberg
SmallPreyIndex<-station.cov2$SmallPreyIndex
PreyIndex<-station.cov2$PreyIndex
MSAVI2ind <- as.numeric(scale(MSAVI2[index1],T,T)>0)

#in Anysberg there are only two plant biomes
#Succelent Karoo and Fynbos
SucculentKarooind <- as.numeric(as.numeric(station.cov$Biome[index1])>1)

det.covs<-as.data.frame(cbind(Trail.ind[index1],
                              plain.index[index1],
                              river.index[index1],
                              SucculentKarooind,
                              c(scale(elevation[index1],T,T)),
                              c(scale(DisturbanceIndex[index1],T,T))))
colnames(det.covs)<-c("trail", "plain", "river", "SKaroo", "elevation", "DisturbanceIndex")


occ.covs<-as.data.frame(cbind(
  Trail.ind[index1],
  plain.index[index1],
  river.index[index1],
  c(scale(elevation[index1],T,T)),
  c(scale(MSAVI2[index1],T,T)),
  c(scale(DisturbanceIndex[index1],T,T)), #added
  c(scale(SmallPreyIndex[index1],T,T)), #added
  c(scale(PreyIndex[index1],T,T)), #added
  c(scale(LivestockIndex[index1],T,T)),
  c(scale(station.cov$PreyIndex[index1],T,T)),
  c(scale(station.cov$SmallPreyIndex[index1],T,T)),
  MSAVI2ind,
  SucculentKarooind) ) #only for Anysberg

colnames(occ.covs)<-c("trail", "plain", "river", "elevation", "MSAVI2",
                      "DisturbanceIndex","SmallPreyIndex","PreyIndex",
                      "LivestockIndex",
                      "PreyIndex",
                      "SmlPrey",
                      "MSAVI2ind",
                      "SKaroo")

longitude = c(scale(station.cov$Longitude[index1], T,T))
latitude = c(scale(station.cov$Latitude[index1],T,T))

#-------------------------------------------------------------------------------

#Example of how to fit the nonspatial model
const <- 2.25^2 #a^2
ndraws <- 1000
MSO.priors <- list(a2 = const, b2 = const, A2 = const, B2 = const)
MSO.control <- list( ndraws = ndraws, percent_burn_in =1/3, thin =1)

#running one chian of 30000 takes 3 minutes
t1<-proc.time()
temp <- MSOBinlogit(detection.model = ~ trail,
                    occupancy.model = ~ elevation  ,
                    occ.covs = occ.covs, det.covs = det.covs,
                    Y = Anysberg.Data,
                    prior = MSO.priors,
                    nsitevisits = matrix(K[Anysberg.rows], ncol=1),
                    control = MSO.control,
                    selection=2)
t2<-proc.time()
t2-t1
#-------------------------------------------------------------------------------
#1. End 'data wrangling'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#2. Begin 'Fit all MSO models'
#-------------------------------------------------------------------------------

#list all of the models to be fitted
library(plyr)
library(ggplot2)

all.subsets <- function(set) {
  #See https://www.r-bloggers.com/2012/04/generating-all-subsets-of-a-set/
  n <- length(set)
  bin <- expand.grid(rlply(n, c(F, T)))
  mlply(bin, function(...) { set[c(...)] })
}

get.reg1 <- function(vars) {
  #See https://www.r-bloggers.com/2012/04/generating-all-subsets-of-a-set/
  if (length(vars) == 0) {
    vars = "1"
  }
  formula(paste(" ~ ", paste(vars, collapse = " + ")))
}

get.reg2 <- function(vars) {
  #See https://www.r-bloggers.com/2012/04/generating-all-subsets-of-a-set/
  formula(paste("~", paste(vars, collapse = " + ")))
  #lm(vars.form, data = data.set, subset = train)
}

varlist_occ <- all.subsets(names(occ.covs)[c(4,1,12,10)])
models_occ1 <- llply(varlist_occ, get.reg1)
models_occ <- models_occ1
lo <- length(models_occ1)
rm(models_occ1)

for (i in 2:lo){
  tmp <- list(formula(paste(models_occ[i]," + plain + river", sep="")))
  models_occ <- append(models_occ, tmp)
  names(models_occ)[length(models_occ)] <- paste(lo+i, sep="")
}

models_det <- list()
models_det$'1'<- formula(~1)
models_det$'2'<- formula(~trail)
models_det$'3'<- formula(~plain + river)
models_det$'4'<- formula(~trail + plain + river)

models_occ_list <- NULL
models_det_list <- NULL

models_occ_list_names <- NULL
models_det_list_names <- NULL

icounter <- 0
for (i in 1:31){
  for (j in 1:4){
    icounter <- 1+ icounter
    models_occ_list[icounter] <- i
    models_det_list[icounter] <- j
    models_occ_list_names[icounter] <- models_occ[i]
    models_det_list_names[icounter] <- models_det[j]
  }
}
models_occ_list_names

#remove 'trail' + 'plain+river' from the occupancy process and keep all of the other models.

findp1<-rep(0,124)
findp2<-rep(0,124)

for(i in 1:124){
  findp1[i] <- ifelse (length(grep("trail",models_occ_list_names[[i]]))>0, i, NA)
  findp2[i] <- ifelse (length(grep("plain",models_occ_list_names[[i]]))>0, i, NA)
}
findp1<-na.omit(findp1)
findp2<-na.omit(findp2)
modindex<-1:124
modindex<-modindex[-unique(c(findp1, findp2))]

for (i in modindex){
 print( models_occ_list_names[[i]] )
}

#the models that we keep!
#these are the 32 models in total
modindex

models_det_list_names[modindex]
models_occ_list_names[modindex]

keptmodels <- NULL
for (i in 1:32){
  keptmodels[i] = paste(models_det_list_names[modindex[i]], "+", models_occ_list_names[modindex[i]])
}
keptmodels

paste(models_det_list_names[modindex], "+", models_occ_list_names[modindex])
# paste(models_det_list_names[modindex], "+", models_occ_list_names[modindex])
#  [1] "~1 + ~1"                                                     "~trail + ~1"
#  [3] "~plain + river + ~1"                                         "~trail + plain + river + ~1"
#  [5] "~1 + ~elevation"                                             "~trail + ~elevation"
#  [7] "~plain + river + ~elevation"                                 "~trail + plain + river + ~elevation"
#  [9] "~1 + ~MSAVI2ind"                                             "~trail + ~MSAVI2ind"
# [11] "~plain + river + ~MSAVI2ind"                                 "~trail + plain + river + ~MSAVI2ind"
# [13] "~1 + ~elevation + MSAVI2ind"                                 "~trail + ~elevation + MSAVI2ind"
# [15] "~plain + river + ~elevation + MSAVI2ind"                     "~trail + plain + river + ~elevation + MSAVI2ind"
# [17] "~1 + ~PreyIndex"                                             "~trail + ~PreyIndex"
# [19] "~plain + river + ~PreyIndex"                                 "~trail + plain + river + ~PreyIndex"
# [21] "~1 + ~elevation + PreyIndex"                                 "~trail + ~elevation + PreyIndex"
# [23] "~plain + river + ~elevation + PreyIndex"                     "~trail + plain + river + ~elevation + PreyIndex"
# [25] "~1 + ~MSAVI2ind + PreyIndex"                                 "~trail + ~MSAVI2ind + PreyIndex"
# [27] "~plain + river + ~MSAVI2ind + PreyIndex"                     "~trail + plain + river + ~MSAVI2ind + PreyIndex"
# [29] "~1 + ~elevation + MSAVI2ind + PreyIndex"                     "~trail + ~elevation + MSAVI2ind + PreyIndex"
# [31] "~plain + river + ~elevation + MSAVI2ind + PreyIndex"         "~trail + plain + river + ~elevation + MSAVI2ind + PreyIndex"
#-------------------------------------------------------------------------------

#obtain the information criteria results
#obtain the information criteria results (Bayesian p value,
#sampled Bayesian p value, WAIC, CPO, elppd, pdWAIC)
#Runs 124 models.
#The code can be adapted so that one only runs the required 32 models

require(MSO)
require(foreach)
require(doParallel)

const <- 2.25^2 #a^2
ndraws <- 100000
MSO.priors <- list(a2 = const, b2 = const, A2 = const, B2 = const)
MSO.control <- list(ndraws = ndraws, percent_burn_in =1/3, thin =10)

cl <- makeCluster(16)
registerDoParallel(cl)

t1<-Sys.time()
models_par <- foreach(i=1:124, .inorder=TRUE, .combine = "cbind", .packages = "MSO")%dopar%{
  #One model takes about 12-13 minutes to run using the above settings

  det.mod <- formula(as.character(models_det[models_det_list[i]]))
  occ.mod <- formula(as.character(models_occ[models_occ_list[i]]))


  tmp <- MSOBinlogit(detection.model = formula(det.mod),
                     occupancy.model = occ.mod ,
                     occ.covs = occ.covs, det.covs = det.covs,
                     Y = Anysberg.Data,
                     prior = MSO.priors,
                     nsitevisits = matrix(K[Anysberg.rows], ncol=1),
                     control = MSO.control,
                     selection=1)

  tmp <- append(list(Model_number=i), tmp)
  write(unlist(tmp), paste("file_",i,".txt",sep=""))
  unlist(tmp)
}

t2<-Sys.time()
t2-t1
stopCluster(cl)

#save the file
save.image("model_selection_Anysberg.RData")
#-------------------------------------------------------------------------------
#2. end 'Fit all MSO models'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Begin 'Model selection diagnostics for the models' (ie 124)
# DELETE THIS SECTION OF CODE!
#-------------------------------------------------------------------------------
dontrun <- TRUE

if (dontrun==TRUE){
  cat("\n Model selection for all the models not run")
}else{
  #Model selection results for all 124 models
  load("model_selection_Anysberg.RData")

  #par(mfrow=c(2,3), mar=c(4,4,1,1), pty="s")
  models_par2 <- models_par[,modindex]

  pdf("Bayespvalue.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    plot(models_par[2,], ylab="Bayesian p-value", ylim=c(0,1), xlab="Model #")
    abline(h=0.5, lty=2, col="red", lwd=2)
    abline(h=0.95, lty=2, col="red", lwd=2)
    abline(h=0.05, lty=2, col="red", lwd=2)

    bpv_pos <- c(which( (models_par[2,]<0.55) & (models_par[2,]>0.45) ))
    points(bpv_pos, models_par[2,bpv_pos], pch=19, col="blue")
    waic_10pc <- quantile(models_par[4,], probs=0.1)
    waic_10pc_pos <- c(which( (models_par[4,]<= waic_10pc) ))
    points(waic_10pc_pos, models_par[2,waic_10pc_pos], pch=3, col="red",cex=1.5)
  dev.off()
  #------------------------------------

  pdf("sBayespvalue.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    plot(models_par[3,], ylab="sampled Bayesian p-value", ylim=c(0,1), xlab="Model #")
    abline(h=0.5, lty=2, col="red", lwd=2)
    abline(h=0.95, lty=2, col="red", lwd=2)
    abline(h=0.05, lty=2, col="red", lwd=2)

    sbpv_pos <- c(which( (models_par[3,]<0.55) & (models_par[3,]>0.45) ))
    points(sbpv_pos, models_par[3,sbpv_pos], pch=19, col="blue")
    points(waic_10pc_pos, models_par[3,waic_10pc_pos], pch=3, col="red", cex=1.5)

    points(25, models_par[3,25], pch=21, col="green",cex=1.5)
    points(17, models_par[3,17], pch=21, col="green",cex=1.5)
    text(28, models_par[3,25], label="25")
    text(20, models_par[3,17], label="17")
  dev.off()
  #------------------------------------

  pdf("pvalues.pdf")
    par(mfrow=c(1,2), mar=c(4,4,1,1), pty="s")
    plot(models_par[2,], ylab="Bayesian p-value", ylim=c(0,1), xlab="Model #")
    abline(h=0.5, lty=2, col="red", lwd=2)
    abline(h=0.95, lty=2, col="red", lwd=2)
    abline(h=0.05, lty=2, col="red", lwd=2)

    bpv_pos <- c(which( (models_par[2,]<0.55) & (models_par[2,]>0.45) ))
    points(bpv_pos, models_par[2,bpv_pos], pch=19, col="blue")
    waic_10pc <- quantile(models_par[4,], probs=0.1)
    waic_10pc_pos <- c(which( (models_par[4,]<= waic_10pc) ))
    points(waic_10pc_pos, models_par[2,waic_10pc_pos], pch=3, col="red",cex=1.5)

    plot(models_par[3,], ylab="sampled Bayesian p-value", ylim=c(0,1), xlab="Model #")
    abline(h=0.5, lty=2, col="red", lwd=2)
    abline(h=0.95, lty=2, col="red", lwd=2)
    abline(h=0.05, lty=2, col="red", lwd=2)

    sbpv_pos <- c(which( (models_par[3,]<0.55) & (models_par[3,]>0.45) ))
    points(sbpv_pos, models_par[3,sbpv_pos], pch=19, col="blue")
    points(waic_10pc_pos, models_par[3,waic_10pc_pos], pch=3, col="red", cex=1.5)

    points(25, models_par[3,25], pch=21, col="green",cex=1.5)
    points(17, models_par[3,17], pch=21, col="green",cex=1.5)
    text(28, models_par[3,25], label="25")
    text(20, models_par[3,17], label="17")
  dev.off()
  #------------------------------------

  pdf("WAIC.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    plot(models_par[4,], ylab="WAIC", xlab="Model #")
    waic_5pc <- quantile(models_par[4,], probs=0.05)
    abline(h=waic_10pc, lty=2, col="red")
    points(waic_10pc_pos, models_par[4,waic_10pc_pos], pch=19, col="red")

    points(25, models_par[4,25], pch=21, col="green",cex=1.5)
    points(17, models_par[4,17], pch=21, col="green",cex=1.5)
    text(28, models_par[4,25], label="25")
    text(20, models_par[4,17], label="17")
  dev.off()
  #------------------------------------

  #CPO
  pdf("CPO.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    plot(models_par[5,], ylab="CPO", xlab="Model #")
    cpo_10pc <- quantile(models_par[5,], probs=0.1)
    abline(h=cpo_10pc, lty=2, col="red")
    cpo_10pc_pos <- c(which( (models_par[5,]<= cpo_10pc) ))
    points(cpo_10pc_pos, models_par[5,cpo_10pc_pos], pch=19, col="red")
    points(25, models_par[5,25], pch=21, col="green",cex=1.5)
    points(17, models_par[5,17], pch=21, col="green",cex=1.5)
    text(28, models_par[5,25], label="25")
    text(20, models_par[5,17], label="17")
  dev.off()
  #------------------------------------

  pdf("CPO_WAIC.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    plot(models_par[4,], models_par[5,], xlab="CPO", ylab="WAIC")
  dev.off()
  #------------------------------------

  pdf("sBayespvalue_WAIC.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    plot(models_par[4,], models_par[3,], ylab="sampled Bayesian p-value", xlab="WAIC", type="n")
    abline(h=0.5, col="red", lty=2, lwd=2)
    abline(h=0.6, col="blue", lty=2, lwd=2)
    abline(h=0.4, col="blue", lty=2, lwd=2)
    abline(v=waic_10pc, lty=2, col="red", lwd=2)

    text(models_par[4,-c(17,25)],
         models_par[3,-c(17,25)],
         labels=c(1:16,18:24,26:124))

    points(models_par[4,25], models_par[3,25], pch=19, col="green",cex=1.5)
    points(models_par[4,17], models_par[3,17], pch=19, col="green",cex=1.5)

    text(models_par[4,c(17,25)],
         models_par[3,c(17,25)],
         labels=c("17","25"))
  dev.off()

  #------------------------------------

  #model weights based on WAIC
  pdf("WAICweights.pdf")
    par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
    sf <- .5#0
    uu <- models_par[4,]- min(models_par[4,])
    tt<-exp(-sf*(uu))
    ss<-tt/sum(tt)
    plot(ss, ylab="Model weight (WAIC)", xlab="Model #", type="n")
    text(1:124,
         ss,
         labels=1:124)
  dev.off()
  order(ss, decreasing=TRUE)
  #------------------------------------

  pdf("WAIC_weights.pdf")
    par(mfrow=c(1,2), mar=c(4,4,1,1), pty="s")
    plot(models_par[4,], ylab="WAIC", xlab="Model #")
    waic_5pc <- quantile(models_par[4,], probs=0.05)
    abline(h=waic_10pc, lty=2, col="red")
    points(waic_10pc_pos, models_par[4,waic_10pc_pos], pch=19, col="red")

    points(25, models_par[4,25], pch=21, col="green",cex=1.5)
    points(17, models_par[4,17], pch=21, col="green",cex=1.5)
    text(28, models_par[4,25], label="25")
    text(20, models_par[4,17], label="17")

    plot(ss, ylab="Model weight (WAIC)", xlab="Model #", type="n")
    text(1:124,
         ss,
         labels=1:124)
  dev.off()
}
#-------------------------------------------------------------------------------
# End 'Model selection diagnostics for the models' (ie 124)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 4. Begin 'Model selection diagnostics for the 'modindex' models' (ie 32)
#-------------------------------------------------------------------------------

#Model selection results for all 124 models
load("model_selection_Anysberg.RData")

#Model selection plots for the models included in 'modindex'
#Use these results!

pdf("sampledBayespvals.pdf")
  par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
  plot(models_par[3, modindex], ylab="sampled Bayesian p-value", xlab="Model #")

  for (i in c(4,8,20,24,32)){
    points(i, models_par[3, modindex[i]], pch=19, col="blue",cex=1.5)
    text(i+.5, models_par[3, modindex[i]]+.035, label=as.character(i))
  }

  abline(h=0.5, lty=2, col="red", lwd=1.5)
  abline(h=0.55, lty=2, col="red", lwd=1) +089
  abline(h=0.45, lty=2, col="red", lwd=1)
  abline(h=0.95, lty=2, col="red", lwd=1)
  abline(h=0.05, lty=2, col="red", lwd=1)
dev.off()

pdf("WAIC2.pdf")
  par(mfrow=c(1,1), mar=c(4,4,1,1), pty="s")
  plot(models_par[4, modindex], ylab="WAIC", xlab="Model #")

  for (i in c(4,8,20,24,32)){
    points(i, models_par[4, modindex[i]], pch=19, col="blue",cex=1.5)
    text(i+.5, models_par[4, modindex[i]]+.035, label=as.character(i))
  }
dev.off()

#-------------------------------------------------------------------------------
# 4. End of Model selection work
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 5. Begin: Identify all of the models fitted and outputted to a file
#-------------------------------------------------------------------------------

#some code to extract the model formulae for a model
#as argument is 'model number' (1:124)

model_string<-function(i){
  paste(
    c(paste(as.character( models_occ_list_names[i][[1]]), collapse = ""),
      paste(as.character(models_det_list_names[i][[1]]), collapse = "")),
    collapse="")
}

sink("FulllistofModels.txt")
  cat("\n All models (occupancy ~ detection)")
  cat("\n ***********************\n")

  for (i in 1:124){
    cat("\n ")
    cat(paste("Model ", i, "= ", model_string(i)) )
  }

  cat("\n")
  cat("\n The 32 models used in the paper (occupancy ~ detection)")
  cat("\n *****************************************\n")

  for (i in modindex){
    cat("\n ")
    cat(paste("Model ", i, "= ", model_string(i)) )
  }

sink()

#-------------------------------------------------------------------------------
# 5. End
#-------------------------------------------------------------------------------

#' #-------------------------------------------------------------------------------
#' # 7. Begin. Obtain the waic and loo estimates for 'modindex' models
#' # this set of results uses the following files:
#' # paste("file_", modindex[modelsi],".RData", sep="")
#' # that has already been saved!
#' #-------------------------------------------------------------------------------
#'
#' #Tabulate the models fitted below:
#'
#' logl_MSO<-function(data_i, draws, K){
#'   #return a vector of loglikelihood terms evaluated for the ith observation
#'   #and mcmc draws
#'   #data_i = a constant
#'   #draws=c(psi, p)
#'   #K = the number of visits to a site
#'   #6666 has been hard coded
#'   #the initial models were run using 100 000 iterations; 2/3 iterations were retained
#'   #and the sequence was thinned by retaining every 10 samples
#'   #thats why 6666 MCMC samples were retained.
#'   lbinom_mass.vec(data_i, psi=draws[c(1:6666)], p=draws[-c(1:6666)], K[1])
#' }
#'
#' binom_mass.vec <- Vectorize(binom_mass, vectorize.args = c("psi","p"))
#' lbinom_mass.vec <- Vectorize(lbinom_mass, vectorize.args = c("psi","p"))
#'
#' require(MSO)
#' require(foreach)
#' require(doParallel)
#' require(loo)
#'
#' t1<-Sys.time()
#'
#' Y. <- t(Anysberg.Data)
#'
#' cl <- makeCluster(16)
#' registerDoParallel(cl)
#'
#' models_loo <- foreach(modelsi=1:32, .inorder=TRUE, .packages = c("MSO","loo"))%dopar%{
#'
#'   load(paste("file_", modindex[modelsi],".RData", sep=""))
#'
#'   Draws <- matrix(0, nrow=156*35, 6666*2)
#'   Ydata <- matrix(0, nrow=156*35, ncol=1)
#'   Kdata <- matrix(0, nrow=156*35, ncol=1)
#'   counter<-0
#'   for (i in 1:156){
#'     for (j in 1:35){
#'       counter<-1+counter
#'       Draws[counter,]<- c(tmp$psi_array[i,j,], tmp$p_array[i,j,])
#'       Ydata[counter]<- Y.[i,j]
#'       Kdata[counter]<-K[Anysberg.rows][i]
#'     }
#'   }
#'
#'   #create a matrix to produce the loglikelihood values
#'   #'number of samples' by 'number of observations'
#'   #this could be made more efficient
#'   logl_m <- matrix(0, nrow=6666, ncol=156*35)
#'
#'   for (i in 1:(156*35)){
#'
#'     logl_m[,i]<- logl_MSO(data_i=Ydata[i,1],
#'                           draws=Draws[i,],
#'                           K=Kdata[i,1])
#'   }
#'
#'   waic_mod<-waic(logl_m)
#'   loo_mod<-loo(logl_m)
#'
#'   write(modelsi, paste("loo_",modelsi,".txt", sep="") )
#'
#'   #waic estimates;
#'   #loo estimates;
#'   #as well as how many y[i,j] observations have
#'   #elements have pareto k estimates >0.7;
#'   #the proportion of observations with a pareto k estimate
#'   #greater than 0.7 is also reported;
#'   #if there are any k>0.7 the observation number is recorded
#'   #k>=1 has a '1' in the second column of the output while
#'   #0.7<k<1 has a '0.7' in the second column
#'
#'   index1 <- which(loo_mod$diagnostics$pareto_k<1 & loo_mod$diagnostics$pareto_k>0.7)
#'   if (length(index1)>0){
#'     index1 <- cbind(index1, 0.7)
#'   }else{
#'     index1 <- c(NA, 0.7)
#'   }
#'
#'   index2 <- which(loo_mod$diagnostics$pareto_k>=1)
#'   if (length(index2)>0){
#'     index2 <- cbind(index2, 1)
#'   }else{
#'     index2 <- c(NA, 1)
#'   }
#'
#'   rbind(waic_mod$estimates, loo_mod$estimates,
#'         c(sum(loo_mod$diagnostics$pareto_k>.7),
#'           sum(loo_mod$diagnostics$pareto_k>.7)/(156*35)),
#'         index1,
#'         index2)
#' }
#' stopCluster(cl)
#'
#' t2<-Sys.time()
#' t2-t1
#' #takes about 15 mins
#'
#' #identify which site, species combinations have large pareto k values
#' kpareto_large<-matrix(0, nrow=156, ncol=35)
#'
#' for (imod in 1:32){
#'
#'   if (models_loo[[imod]][7,1]>0){
#'     bigks <- models_loo[[imod]][-c(1:7),1]
#'
#'     for (ik in 1:length(bigks)){
#'
#'       #identify row and column of points with large pareto k's
#'       Nrow <- ceiling(bigks[ik]/35)
#'
#'       if (Nrow == 1){
#'         Ncol <- bigks[ik]
#'       }else{
#'         Ncol <- bigks[ik] - 35*(Nrow-1)
#'       }#endif
#'
#'       kpareto_large[Nrow, Ncol] <- kpareto_large[Nrow, Ncol] + 1
#'     }
#'   }
#' }
#'
#' #
#'
#' par(mfrow=c(1,1), pty="s")
#' image(1:156, 1:35, kpareto_large/sum(kpareto_large), col = hcl.colors(10, "terrain"))
#'
#' require(rgl)
#' plot3d(1:156, 1:35, kpareto_large)
#'
#' loo_model_scores<-NULL
#' for (i in 1:32){ loo_model_scores[i] <- models_loo[[i]][6,1] }
#'
#' #> order(loo_model_scores)
#' #[1]  9 13  1  5 29 21 25 14 10 17  6  2 30 22 15 26  7 31 23 18 16 27  8 11 32 24 19 28 12  3 20  4
#'
#' kpareto<-NULL
#' for (imod  in 1:32){
#'   kpareto[imod] <- models_loo[[imod]][7,1]
#' }
#' kpareto#/(156*35)
#'
#' yy<-NULL
#' for (i in 1:32){ yy[i] <- models_loo[[i]][3,1] }
#' loo_order <-order(yy)
#'
#' #-------------------------------------------------------------------------------
#' #-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 6. Begin. Investigate the variability of various information criteria
# obtain the information criteria results (Bayesian p value,
# sampled Bayesian p value, WAIC, CPO, elppd, pdWAIC)
# Below I assess the variability of the above model selection
# methods by using 10 replications and runs all 124 models.
# this was run on a desktop computer that had 16 cores
# start: 2019/11/18 12:42
# end: 2019/11/20 2:25am
# so this takes a while to run
# but I ran it using all 124 models!!!
#-------------------------------------------------------------------------------

require(MSO)
require(foreach)
require(doParallel)

MonteCarlo_stats <- list(NULL)

for (run_iteration in 1:10){
  cl <- makeCluster(16)
  registerDoParallel(cl)

  const <- 2.25^2 #a^2
  ndraws <- 100000
  MSO.priors <- list(a2 = const, b2 = const, A2 = const, B2 = const)
  MSO.control <- list( ndraws = ndraws, percent_burn_in =1/3, thin =10)

  t1<-Sys.time()
  models_par_vars <- foreach(i=1:124, .inorder=TRUE, .combine = "cbind", .packages = "MSO")%dopar%{

    det.mod <- formula(as.character(models_det[models_det_list[i]]))
    occ.mod <- formula(as.character(models_occ[models_occ_list[i]]))

    tmp <- MSOBinlogit(detection.model = formula(det.mod),
                       occupancy.model = occ.mod ,
                       occ.covs = occ.covs, det.covs = det.covs,
                       Y = Anysberg.Data,
                       prior = MSO.priors,
                       nsitevisits = matrix(K[Anysberg.rows], ncol=1),
                       control = MSO.control,
                       selection=1)

    tmp <- append(list(Model_number=i), tmp)
    #write(unlist(tmp), paste("Var_file_",i,"_run_",1,".txt",sep=""))
    unlist(tmp)
  }
  t2<-Sys.time()
  t2-t1
  stopCluster(cl)

  MonteCarlo_stats[[run_iteration]] <- models_par_vars
  write(run_iteration,paste(run_iteration,".txt"))
}

#save the file
save.image("MonteCarlo_model_selection_Anysberg.RData")
#-------------------------------------------------------------------------------
# 6. End.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Begin. Investigate the variability of the information criteria using 10
# replicate data sets
#-------------------------------------------------------------------------------

load("MonteCarlo_model_selection_Anysberg.RData")
require(xtable)
allmods_sBpval<-matrix(0, 10, 124)
allmods_Bpval<-matrix(0, 10, 124)
allmods_WAIC<-matrix(0, 10, 124)
allmods_elppd<-matrix(0, 10, 124)
allmods_pdWAIC<-matrix(0, 10, 124)
allmods_CPO<-matrix(0, 10, 124)

for (i in 1:10){
  for (j in 1:124){
    allmods_Bpval[i,j]<- MonteCarlo_stats[[i]][2,j]
    allmods_sBpval[i,j]<- MonteCarlo_stats[[i]][3,j]
    allmods_WAIC[i,j]<- MonteCarlo_stats[[i]][4,j]
    allmods_CPO[i,j]<- MonteCarlo_stats[[i]][5,j]
    allmods_elppd[i,j]<- MonteCarlo_stats[[i]][6,j]
    allmods_pdWAIC[i,j]<- MonteCarlo_stats[[i]][7,j]
  }
}

#rm(MonteCarlo_stats)

# load("MonteCarlo_model_selection_Anysberg2.RData")
#
# allmods_sBpval2<-matrix(0, 10, 124)
# allmods_Bpval2<-matrix(0, 10, 124)
# allmods_WAIC2<-matrix(0, 10, 124)
# allmods_elppd2<-matrix(0, 10, 124)
# allmods_pdWAIC2<-matrix(0, 10, 124)
#
# for (i in 1:10){
#   for (j in 1:124){
#     allmods_Bpval2[i,j]<- MonteCarlo_stats[[i]][2,j]
#     allmods_sBpval2[i,j]<- MonteCarlo_stats[[i]][3,j]
#     allmods_WAIC2[i,j]<- MonteCarlo_stats[[i]][4,j]
#     allmods_elppd2[i,j]<- MonteCarlo_stats[[i]][6,j]
#     allmods_pdWAIC2[i,j]<- MonteCarlo_stats[[i]][7,j]
#   }
# }
#
# allmods_sBpval <- rbind(allmods_sBpval, allmods_sBpval2)
# allmods_Bpval <- rbind(allmods_Bpval, allmods_pval2)
# allmods_WAIC <- rbind(allmods_WAIC, allmods_WAIC2)
# allmods_elppd <- rbind(allmods_WAIC, allmods_elppd2)
# allmods_pdWAIC <- rbind(allmods_pdWAIC, allmods_pdWAIC2)

#here
par(mfrow=c(1,1), pty="m")

#decided to keep the order from 1:32
WAIC_index <- order(apply(allmods_WAIC[,modindex],2,median)) #order(yy)#

colindex3 <- rep(NA,32)
colindex3[c(1,4,9,11,19,21,22,23,25,26,28,31)] <- "grey"
colindex3[c(1,9,11)] <- "red"

#plot of the boxplots of the sample Bayesian p-values
boxplot(allmods_sBpval[,modindex][,WAIC_index],
        ylim=c(0,1),
        col=colindex3,
        names=WAIC_index,
        xlab="Model number (sorted according to LOO (increasing order))",
        ylab="Sampled Bayesian p-value (based on 10 MCMC runs)",
        cex.lab=1.5, las=1, cex.axis=.75)

#plot of the boxplots of the Bayesian p-values
boxplot(allmods_Bpval[,modindex][,WAIC_index],
        ylim=c(0,1),
        col=colindex3,
        names=WAIC_index,
        xlab="Model number (sorted according to LOO (increasing order))",
        ylab="Sampled Bayesian p-value (based on 10 MCMC runs)",
        cex.lab=1.5, las=1, cex.axis=.75)

#plot of the boxplots of the Bayesian p-values
boxplot(allmods_CPO[,modindex],
        col=colindex3,
        names=1:32,
        xlab="Model number (sorted according to LOO (increasing order))",
        ylab="Sampled Bayesian p-value (based on 10 MCMC runs)",
        cex.lab=1.5, las=1, cex.axis=.75)


med_BPval <- apply(allmods_Bpval[,modindex], 2, median)

#-------------------------------------------------------------------------------

# perform Anderson-Darling tests using the 10 samples of the sampled Bayesian p
# values in order to assess whether they are approximately uniformly distributed
# Cheng and Spring 1987
# A Test to Identify the Uniform Distribution, with Applications to
# Probability Plotting and Other Distributions
# is also to test as an alternative

wu <- function(ns){
  #Used to obtain the sampling distribution of Shapiro Wilks test statistic
  #used to test whether x is uniformly distributed
  #See Cheng and Spring 1987
  #A Test to Identify the Uniform Distribution, with Applications to
  #Probability Plotting and Other Distributions

  x <- runif(ns)
  n <- length(x)
  t1 <- (max(x) - min(x))*(n+1)/(n-1)
  t1 <- t1^2
  t2 <- sum((x-mean(x))^2)
  t1/t2
}

wut <- function(x, distr){
  #Used to calculate the Shapiro Wilks test statistic
  #used to test whether x is uniformly distributed
  #See Cheng and Spring 1987
  #A Test to Identify the Uniform Distribution, with Applications to
  #Probability Plotting and Other Distributions

  #x <- runif(10)
  n <- length(x)
  t1 <- (max(x) - min(x))*(n+1)/(n-1)
  t1 <- t1^2
  t2 <- sum((x-mean(x))^2)
  stat=t1/t2

  c(stat, sum(distr>stat)/length(distr))
}

#The sampling distribution of a uniform distribution based on 10 observations
wudist <- replicate(1000000, wu(10))
hist(wudist, breaks=100)

#Tests
require(goftest)
ResultsAD<-matrix(0, 32, 2)
ResultsWu<-matrix(0, 32, 2)

temp <- allmods_sBpval[,modindex]

#the ad.test uses the Braun's method of randomly calculating the pvalue
#of the Anderson Darling test. There are choose(10,5)=252 different random
#splits of 5. Here I use 500 random splits understanding that some of the
#splits will be the same
#This takes a few minutes to run
#Note that these results can be random since the test statistic is constructed
#by using random splits of the data!
for (repl in 1:500){
  for (i in 1:32){
    tempv <- temp[,i]
    test <- wut(tempv, wudist)
    ResultsWu[i, ] <- ResultsWu[i, ] + test

    #test <- cvm.test(tempv, null = "punif", estimated=TRUE)
    #ResultsWu[i, ] <- c(test$statistic, test$p.value)

    #Anderson Darling test for uniform distribution
    test <-ad.test(tempv, null = "punif", estimated=TRUE)
    ResultsAD[i, ] <- ResultsAD[i, ] +c(test$statistic, test$p.value)
  }

  #cat("\n repl=", repl)
}

ResultsWu <- ResultsWu/500
ResultsWu <- cbind(ResultsWu, 1:32)
ResultsWu <- ResultsWu[,c(3,1,2)]

rownames(ResultsWu) <- paste(models_det_list_names[modindex], "+", models_occ_list_names[modindex])
#ResultsWu

ResultsAD <- ResultsAD/500
ResultsAD <- cbind(ResultsAD, 1:32)
ResultsAD <- ResultsAD[,c(3,1,2)]

rownames(ResultsAD) <- paste(models_occ_list_names[modindex],  models_det_list_names[modindex])
#ResultsAD

ResultsWu[order(ResultsWu[,3]), ]
ResultsAD[order(ResultsAD[,3]), ]

pdf("ModelSelection_sBpvals.pdf")
  par(mfrow=c(1,1), pty="m")
  colindexAD <- rep("NA",32)
  colindexAD[which(round(ResultsAD[,3],2)>=.3) ] <- "grey"

  boxplot(allmods_sBpval[,modindex],
          ylim=c(0,1),
          col=colindexAD,
          names=1:32,
          xlab="Model number",
          ylab="Sampled Bayesian p-value",
           xaxt="n")
  #axis(side=1, at=seq(from=2, to=32, by=2), labels = FALSE)
  text(x=seq(from=1, to=32, by=1),  par("usr")[3],
       labels = seq(from=1, to=32, by=1), pos = 1, xpd = TRUE,
       cex=.75)

  #for (i in 1:7){abline(v=4.5+4*(i-1), lty=2)}
dev.off()

xtable(round(ResultsAD[which(ResultsAD[,3]>=0.3),c(1,3)], digits=3))
ADfail2rej <- round(ResultsAD[which(ResultsAD[,3]>=0.3),], digits=3)
round(ResultsAD, digits=2)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Tabulate the Bayesian p values, sampled bayesian p values and waic values
#-------------------------------------------------------------------------------

require(xtable)
model.selection.table <- cbind(
  cbind( apply( t(allmods_Bpval)[modindex,], 1, mean),
         apply( t(allmods_Bpval)[modindex,], 1, sd )/sqrt(10),
         apply( t(allmods_Bpval)[modindex,], 1, median) ),

  cbind( apply( t(allmods_sBpval)[modindex,], 1, mean),
         apply( t(allmods_sBpval)[modindex,], 1, sd )/sqrt(10),
         apply( t(allmods_sBpval)[modindex,], 1, median) ),

  cbind( apply( t(allmods_WAIC)[modindex,], 1, mean),
         apply( t(allmods_WAIC)[modindex,], 1, sd )/sqrt(10),
         apply( t(allmods_WAIC)[modindex,], 1, median) )
)


model.selection.table <- round(model.selection.table, digits=2)
model.selection.table[,2] <- round(apply( t(allmods_Bpval)[modindex,], 1, sd )/sqrt(10), 3)

colnames(model.selection.table)<-rep(c("Mean", "SE", "Median"), 3)

model.selection.table.rows <- matrix("1", nrow=32, ncol=2)
colnames(model.selection.table.rows) <- c("Occ", "Det")

models_occ_list_names <- NULL
models_det_list_names <- NULL

icounter <- 0
for (i in 1:31){
  for (j in 1:4){
    icounter <- 1+ icounter
    models_occ_list[icounter] <- i
    models_det_list[icounter] <- j
    models_occ_list_names[icounter] <- models_occ[i]
    models_det_list_names[icounter] <- models_det[j]
  }
}

#list of occupancy and detection models considered for the analysis
models_det_list_names2 <- models_det_list_names[modindex]
models_occ_list_names2 <- models_occ_list_names[modindex]

for (i in 1:32){
  model.selection.table.rows[i,1]<- as.character(models_occ_list_names2[[i]][2])
  model.selection.table.rows[i,2]<- as.character(models_det_list_names2[[i]][2])
}

boot.modelselection<-noquote( cbind(model.selection.table.rows, model.selection.table) )

#the models fitted - their names
#for all the models
print.xtable( xtable(cbind(round(1:32,0), model.selection.table),
                     digits=c(0,0,2,3,rep(2,7))), include.rownames=FALSE)

model.selection.table[ADfail2rej[,1],9]
save.image("Section8.RData")
# 8. End.
#-------------------------------------------------------------------------------


# I am here
#-------------------------------------------------------------------------------
# 9. Begin. RJMCMC code
#-------------------------------------------------------------------------------

#Here I code up a RJMCMC algorithm
load("Section8.RData")

models_occ_list_names <- NULL
models_det_list_names <- NULL

icounter <- 0
for (i in 1:31){
  for (j in 1:4){
    icounter <- 1+ icounter
    models_occ_list[icounter] <- i
    models_det_list[icounter] <- j
    models_occ_list_names[icounter] <- models_occ[i]
    models_det_list_names[icounter] <- models_det[j]
  }
}

#list of occupancy and detection models considered for the analysis
models_det_list_names2 <- models_det_list_names[modindex]
models_occ_list_names2 <- models_occ_list_names[modindex]
#-------------------------------------------------------------------------------

require(mvnfast)

#Functions required to perform the RJMCMC
log_var_dens <- function(x, mu_x, a2, var_max, ncores=1){
  #x = simulations for a variable (alpha_i or beta_i)
  #mu_x = posterior mean of 'x'
  #a2 = variance used for 'x'
  #var_max = the maximum number of covariates per process

  dim_mu <- length(x)
  #mu <- matrix(mu_x, ncol=1)
  sigma <- diag(dim_mu)*a2
  temp1 <- dmvn(x, mu_x, sigma, log = TRUE, ncores = ncores)
  #cat("\n temp1 done \n")

  temp2 <- 0
  dim_diff<-var_max-dim_mu

  if (dim_diff>0){
    temp2 <- dnorm(rnorm(dim_diff), log = TRUE)
  }

  return(temp1 + temp2)
}

# #Y=Anysberg.Data,
# idx = k.rand.seq[1]
# alpha=tmp$alpha[,,idx]
# beta=tmp$beta[,,idx]
# tau_alpha=tmp$tau_alpha[idx]
# tau_beta=tmp$tau_beta[idx]
# Xdet=det.covs2
# Xocc=occ.covs2
# det.mat=det.mat
# occ.mat=occ.mat
# alpha_max=4
# beta_max=4
# M=M
# #mu_alpha...
# mu_alpha=tmp$mu_alpha[,idx]
# mu_beta=tmp$mu_beta[,idx]
# const=const
# ncores=1

M_l_contr <- function(Y,
                      alpha, beta,
                      tau_alpha, tau_beta,
                      Xdet, Xocc,
                      det.mat, occ.mat,
                      alpha_max, beta_max,
                      M,
                      #for mu_alpha contribution
                      mu_alpha,
                      mu_beta,
                      const,
                      ncores=1){
  #calculates the 'log joint pmf' contribution to Pr(Mi|.) using simulated samples
  #Y = data (matrix of number of detections)
  #alpha = simulated (MCMC) detection regression effects for all species
  #beta = simulated (MCMC) occupancy regression effects for all species
  #Xocc = the design matrix associated with the occupancy model
  #Xdet = the design matrix associated with the detection model
  #det.mat = a matrix of 0's and 1's that indicates what all the detection models are
  #occ.mat= a matrix of 0's and 1's that indicates what all the occupancy models are
  #alpha_max = the maximum number of alpha regression coefs
  #beta_max = the maximum number of beta regression coefs
  #M = the model being iterated over
  #mu_alpha =
  #mu_beta =
  #const =
  #ncores = 1
  #a number of objects below not defined as arguments should be in the workspace!
  #-----------------------------------------------------------------------------

  #iterated model
  #mod_coefs = index indicating which model to use
  #det.covs = detection cov matrix
  #identify the names of the detection covariates in the model used
  det.mod.vars <- colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))

  #identify which of of the coefficients in model Ml is in model M=M
  if (length(det.mod.vars)==1){
    #if only an intercept is fitted
    det_index <- (det.mat[M,]==1) #M = the model number
    theta_a_index <- 1

  }else{
    rownames(alpha) <- det.mod.vars
    det_index <- (det.mat[M,]==1) #M = the model number
    theta_a_index <- which( (det_index + (det.mat[which(modindex==mod_coefs),]))==2)
    temp<-sapply(colnames(det.mat)[theta_a_index], function(x)which(colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))==x))
    theta_a_index <- temp
  }
  #-----------------------------------------------------------------------------

  #identify the names of the occupancy covariates in the model used
  occ.mod.vars<-colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))

  if (length(occ.mod.vars)==1){
    occ_index <- (det.mat[M,]==1)
    theta_b_index <- 1
  }else{
    occ_index <- (occ.mat[M,]==1)
    rownames(beta) <- occ.mod.vars
    theta_b_index <- which( (occ_index + (occ.mat[which(modindex==mod_coefs),]))==2)
    temp<-sapply(colnames(occ.mat)[theta_b_index], function(x)which(colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))==x))
    theta_b_index <- temp
  }
  #-----------------------------------------------------------------------------

  # cat("\n det_index=", det_index)
  # cat("\n theta_a_index=", theta_a_index)
  #
  # cat("\n occ_index=", occ_index)
  # cat("\n theta_b_index=", theta_b_index)
  # cat("\n")
  #-----------------------------------------------------------------------------

  #calculation of the log likelihood contribution
  #cat("\n start logl ...")
  #'35' is the number of species in my application
  #place N(0,1) random variables in theta_
  theta_a1<- matrix(rnorm(alpha_max*35),nrow=alpha_max, ncol=35)
  theta_b1<- matrix(rnorm(beta_max*35),nrow=beta_max, ncol=35)

  #replace appropriate values of theta_ with the alpha values
  if (length(det.mod.vars)==1){
    #if the intercept model is selected
    #here 'alpha' is a vector
    theta_a1[theta_a_index,] <- alpha[theta_a_index]
  }else{
    #here 'alpha' is a matrix
    theta_a1[theta_a_index,] <- alpha[theta_a_index,]
  }

  if (length(occ.mod.vars)==1){
    #if the intercept model is selected
    #'beta' is a vector
    theta_b1[theta_b_index,] <- beta[theta_b_index]
  }else{
    theta_b1[theta_b_index,] <- beta[theta_b_index,]
  }
  #-----------------------------------------------------------------------------

  Y <- t(Y) #nsites by species

  #occupancy probs
  if (NCOL(theta_b1[occ_index,])==1){
    #only one variable in the iterated model
    psi_ij <- 1/(1+exp(-matrix(Xocc[,occ_index], ncol=1)%*%matrix(theta_b1[occ_index,], nrow=1))) #nsites by species
  }else{
    psi_ij <- 1/(1+exp(-Xocc[,occ_index]%*%theta_b1[occ_index,])) #nsites by species
  }

  #detection probs
  if (NCOL(theta_a1[det_index,])==1){
    #only one variable in the iterated model
    p_ij <- 1/(1+exp(-matrix(Xdet[,det_index], ncol=1)%*%matrix(theta_a1[det_index,], nrow=1)))
  }else{
    p_ij <- 1/(1+exp(-Xdet[,det_index]%*%theta_a1[det_index, ])) #nsites by species
  }

  logl <-  sum( log( dbinom(Y,  K2, p_ij)*psi_ij + (1-psi_ij)*(Y==0) ))
  #cat("\n end logl ...")
  #cat("\n")
  #-----------------------------------------------------------------------------

  #mu_alpha, u_alpha, mu_beta and u_beta contributions
  #inefficient code but it works!
  #cat("\n start mu_ contributions ...")
  theta_a2<- matrix(rnorm(alpha_max),nrow=alpha_max, ncol=1)
  theta_a2[theta_a_index] <- mu_alpha[theta_a_index]

  dim_mu <- length(theta_a2[theta_a_index])
  mu <- matrix(0, dim_mu, ncol=1)
  sigma <- diag(dim_mu)*const

  #cat("\n mu=", mu)
  #cat("\n sigma=", sigma)
  #cat("\n theta_a2[theta_a_index]=", theta_a2[theta_a_index])
  temp1 <- dmvn(theta_a2[theta_a_index], mu, sigma, log = TRUE, ncores = ncores)

  temp2 <- 0
  dim_diff<-alpha_max-dim_mu

  if (dim_diff>0){
    temp2 <- dnorm(rnorm(dim_diff), log = TRUE)
  }

  mu_alpha_c <- sum(temp1) + sum(temp2)
  #-----------------------------------

  theta_b2<- matrix(rnorm(beta_max),nrow=beta_max, ncol=1)
  theta_b2[theta_b_index,] <- mu_beta[theta_b_index]

  dim_mu <- length(theta_b2[theta_b_index])
  mu <- matrix(0, dim_mu, ncol=1)
  sigma <- diag(dim_mu)*const
  temp3 <- dmvn(theta_b2[theta_b_index], mu, sigma, log = TRUE, ncores = ncores)

  temp4 <- 0
  dim_diff<-beta_max-dim_mu

  if (dim_diff>0){
    temp4 <- dnorm(rnorm(dim_diff), log = TRUE)
  }

  mu_beta_c <- sum(temp3) + sum(temp4)
  #cat("\n end mu_ contributions ...")
  #cat("\n")
  #-----------------------------------


  #alpha_i, beta_i, u_alpha_i and u_beta_i contributions
  #cat("\n start _i contributions ...")
  if (length(theta_a_index)==1){
    #cat("\n alpha1")
    #cat("\n theta_a1[theta_a_index,]=",theta_a1[theta_a_index,])
    #cat("\n theta_a2[theta_a_index,]=",theta_a2[theta_a_index,])
    #cat("\n tau_alpha[idx]=",tau_alpha[idx])

    mu_alpha_i_c<-sum(apply(matrix(theta_a1[theta_a_index,], nrow=1),
                            2,
                            FUN=log_var_dens,
                            mu_x=matrix(theta_a2[theta_a_index], nrow=1),
                            a2=1/tau_alpha,
                            var_max=alpha_max,
                            ncores))
  }else{
    #if the iterated model has more than 1 variable
    #cat("\n alpha2")
    mu_alpha_i_c<-sum(apply(theta_a1[theta_a_index,],
                            2,
                            FUN=log_var_dens,
                            mu_x=matrix(theta_a2[theta_a_index], ncol=1),
                            a2=1/tau_alpha,
                            var_max=alpha_max,
                            ncores))
  }

  if (length(theta_b_index)==1){
    #cat("\n beta1")
    mu_beta_i_c<-sum(apply(matrix(theta_b1[theta_b_index,], nrow=1),
                           2,
                           FUN=log_var_dens,
                           matrix(theta_b2[theta_b_index], nrow=1),
                           1/tau_beta,
                           beta_max,
                           ncores))
  }else{
    #if the iterated model has more than 1 variable
    #cat("\n beta2")
    mu_beta_i_c<-sum(apply(theta_b1[theta_b_index,],
                           2,
                           FUN=log_var_dens,
                           matrix(theta_b2[theta_b_index], ncol=1),
                           1/tau_beta,
                           beta_max,
                           ncores))
  }
  #cat("\n end _i contributions ...")
  #cat("\n")

  return(logl + mu_alpha_c + mu_beta_c + mu_alpha_i_c + mu_beta_i_c)
}

#-------------------------------------------------------------------------------

#list of occupancy and detection models considered for the analysis
models_det_list_names2 <- models_det_list_names[modindex]
models_occ_list_names2 <- models_occ_list_names[modindex]

#create an matrix of indices for the different models
det.mat<-matrix(0,32,4)
det.mat[1,]<-c(1,0,0,0)
det.mat[2,]<-c(1,1,0,0)
det.mat[3,]<-c(1,0,1,1)
det.mat[4,]<-c(1,1,1,1)
det.mat[5,]<-c(1,0,0,0)
det.mat[6,]<-c(1,1,0,0)
det.mat[7,]<-c(1,0,1,1)
det.mat[8,]<-c(1,1,1,1)
det.mat[9,]<-c(1,0,0,0)
det.mat[10,]<-c(1,1,0,0)
det.mat[11,]<-c(1,0,1,1)
det.mat[12,]<-c(1,1,1,1)
det.mat[13,]<-c(1,0,0,0)
det.mat[14,]<-c(1,1,0,0)
det.mat[15,]<-c(1,0,1,1)
det.mat[16,]<-c(1,1,1,1)
det.mat[17,]<-c(1,0,0,0)
det.mat[18,]<-c(1,1,0,0)
det.mat[19,]<-c(1,0,1,1)
det.mat[20,]<-c(1,1,1,1)
det.mat[21,]<-c(1,0,0,0)
det.mat[22,]<-c(1,1,0,0)
det.mat[23,]<-c(1,0,1,1)
det.mat[24,]<-c(1,1,1,1)
det.mat[25,]<-c(1,0,0,0)
det.mat[26,]<-c(1,1,0,0)
det.mat[27,]<-c(1,0,1,1)
det.mat[28,]<-c(1,1,1,1)
det.mat[29,]<-c(1,0,0,0)
det.mat[30,]<-c(1,1,0,0)
det.mat[31,]<-c(1,0,1,1)
det.mat[32,]<-c(1,1,1,1)
colnames(det.mat) <- c("(Intercept)", "trail", "plain", "river")

occ.mat<-matrix(0,32,4)
occ.mat[1,]<-c(1,0,0,0)
occ.mat[2,]<-c(1,0,0,0)
occ.mat[3,]<-c(1,0,0,0)
occ.mat[4,]<-c(1,0,0,0)
occ.mat[5,]<-c(1,1,0,0)
occ.mat[6,]<-c(1,1,0,0)
occ.mat[7,]<-c(1,1,0,0)
occ.mat[8,]<-c(1,1,0,0)
occ.mat[9,]<-c(1,0,1,0)
occ.mat[10,]<-c(1,0,1,0)
occ.mat[11,]<-c(1,0,1,0)
occ.mat[12,]<-c(1,0,1,0)
occ.mat[13,]<-c(1,1,1,0)
occ.mat[14,]<-c(1,1,1,0)
occ.mat[15,]<-c(1,1,1,0)
occ.mat[16,]<-c(1,1,1,0)
occ.mat[17,]<-c(1,0,0,1)
occ.mat[18,]<-c(1,0,0,1)
occ.mat[19,]<-c(1,0,0,1)
occ.mat[20,]<-c(1,0,0,1)
occ.mat[21,]<-c(1,1,0,1)
occ.mat[22,]<-c(1,1,0,1)
occ.mat[23,]<-c(1,1,0,1)
occ.mat[24,]<-c(1,1,0,1)
occ.mat[25,]<-c(1,0,1,1)
occ.mat[26,]<-c(1,0,1,1)
occ.mat[27,]<-c(1,0,1,1)
occ.mat[28,]<-c(1,0,1,1)
occ.mat[29,]<-c(1,1,1,1)
occ.mat[30,]<-c(1,1,1,1)
occ.mat[31,]<-c(1,1,1,1)
occ.mat[32,]<-c(1,1,1,1)
colnames(occ.mat)<- c( "(Intercept)", "elevation" , "MSAVI2ind", "PreyIndex"  )

#the occupancy and detection covariates used
occ.covs2 <- as.matrix(cbind(1, occ.covs[,c("elevation", "MSAVI2", "PreyIndex")]))
det.covs2 <- as.matrix(cbind(1, det.covs[,c("trail", "plain","river")]))
#------------------------------------

#This is a single pass: with a given model start, obtain the posterior prob
#of each model using the RJMCMC algorithm
# #the maximum number of coefficients in detection and occupancy process
# alpha_max <- 4
# beta_max <- 4
# ncores <- 1
#
# #the data set of detections
# K2 <- matrix(K[Anysberg.rows], nrow=156, ncol=35)
#
# #the betas comes from the 'tmp' object
# #this comes from modindex
# if (exists("tmp")==TRUE){
#   rm("tmp")
# }else{
#   paste(" 'tmp' does not exist. Nothing is deleted.")
# }
# #------------------------------------
#
# #serial implementation of the RJMCMC algorithm
# #6666 simulations were kept!
# #use the full model as a start of the RJMCMC algorithm
# mod_coefs <- 53
# load(paste("file_",mod_coefs,".RData", sep=""))
# colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
# colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))
#
# if (mod_coefs==1){
#   set.seed(2) #used for mod_coefs=1
#   model_sel <-rep(which(modindex==1), 6666)
#   mod <- 1
# }
#
# if (mod_coefs==56){
#   set.seed(20) #used for mod_coefs=56
#   model_sel <-rep(which(modindex==56), 6666)
#   mod <- 56
# }
#
# if (mod_coefs==17){
#   set.seed(200) #used for mod_coefs=17
#   model_sel <-rep(which(modindex==17), 6666)
#   mod <- 17
# }
#
# if (mod_coefs==53){
#   set.seed(2000) #used for mod_coefs=29
#   model_sel <-rep(which(modindex==53), 6666)
#   mod <- 53
# }
#
#
# k.rand.seq <- sample(2:6666, 6665)
# tmp.sum <- rep(0, 32)
#
# icounter<-0
# mod_coefs_b <- mod_coefs
#
# t1 <- Sys.time()
#
# for (idx in k.rand.seq){
#
#   icounter<- icounter+1
#   #iterate over all of the models
#   for (M in 1:32){
#
#     tmp.sum[M]<-M_l_contr(Y=Anysberg.Data,
#                           alpha=tmp$alpha[,,idx], beta=tmp$beta[,,idx],
#                           tau_alpha=tmp$tau_alpha[idx], tau_beta=tmp$tau_beta[idx],
#                           Xdet=det.covs2, Xocc=occ.covs2,
#                           det.mat=det.mat, occ.mat=occ.mat,
#                           alpha_max=4, beta_max=4,
#                           M=M,
#                           #mu_alpha...
#                           tmp$mu_alpha[,idx],
#                           tmp$mu_beta[,idx],
#                           const,
#                           ncores=1)
#   }
#
#   P.M <- exp(tmp.sum-max(tmp.sum) -log(sum(exp(tmp.sum-max(tmp.sum)))) )
#   model_sel[idx]<- sample(1:32, 1, prob=P.M)
#
#   mod_coefs<-modindex[model_sel[idx]]
#
#   if (mod_coefs_b != mod_coefs){
#     load(paste("file_",mod_coefs,".RData", sep=""))
#     mod_coefs_b <- mod_coefs
#   }
#
#   if (icounter==50*round(icounter/50,0)) { cat("\n Iteration = ", icounter) }
# }
# t2<-Sys.time()
# t2-t1
#
# if (mod==1){
#   #starting value is the null model
#   pmods_1<-table(model_sel[floor(6666/2):6666])/length(floor(6666/2):6666)
#   plot(pmods_1)
#   #save the file
#   save.image("RJMCMC2_Anysberg_1.RData")
# }
#
# if (mod==56){
#   #starting values is the full model
#   pmods_56<-table(model_sel[floor(6666/2):6666])/length(floor(6666/2):6666)
#   plot(pmods_56)
#   #save the file
#   save.image("RJMCMC2_Anysberg_56.RData")
# }
#
# if (mod==17){
#   #starting values is the full model
#   pmods_17<-table(model_sel[floor(6666/2):6666])/length(floor(6666/2):6666)
#   plot(pmods_17)
#   #save the file
#   save.image("RJMCMC2_Anysberg_17.RData")
# }
#
# if (mod==53){
#   #this is model 29
#   #starting values is the full model
#   pmods_53<-table(model_sel[floor(6666/2):6666])/length(floor(6666/2):6666)
#   plot(pmods_53)
#   #save the file
#   save.image("RJMCMC2_Anysberg_53.RData")
# }
#
# #------------------------------------

# #alternative method
# #here we attempt to estimate the transition prob matrix
#
# #the maximum number of coefficients in detection and occupancy process
# alpha_max <- 4
# beta_max <- 4
# ncores <- 1
#
# #the data set of detections
# K2 <- matrix(K[Anysberg.rows], nrow=156, ncol=35)
#
# #the betas comes from the 'tmp' object
# #this comes from modindex
# if (exists("tmp")==TRUE){
#   rm("tmp")
# }else{
#   paste(" 'tmp' does not exist. Nothing is deleted.")
# }
#
# #serial implementation of the RJMCMC algorithm
# #6666 simulations were kept!
# #use the full model as a start of the RJMCMC algorithm
# mod_coefs <- 1
# load(paste("file_",mod_coefs,".RData", sep=""))
# colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
# colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))
#
# #model_sel <-rep(which(modindex==mod_coefs), 6666)
#
# k.rand.seq <- sample(2:6666, 6665)
# tmp.sum <- rep(0, 32)
#
# icounter<-0
# #mod_coefs_b <- mod_coefs
#
# t1 <- Sys.time()
#
# P.M <- rep(0, 32)
#
# #for (idx in k.rand.seq){
# for (idx in 1:100){
#   icounter<- icounter+1
#
#   #iterate over all of the models
#   for (M in 1:32){
#
#     tmp.sum[M]<-M_l_contr(Y=Anysberg.Data,
#                           alpha=tmp$alpha[,,idx], beta=tmp$beta[,,idx],
#                           tau_alpha=tmp$tau_alpha[idx], tau_beta=tmp$tau_beta[idx],
#                           Xdet=det.covs2, Xocc=occ.covs2,
#                           det.mat=det.mat, occ.mat=occ.mat,
#                           alpha_max=4, beta_max=4,
#                           M=M,
#                           #mu_alpha...
#                           tmp$mu_alpha[,idx],
#                           tmp$mu_beta[,idx],
#                           const,
#                           ncores=1)
#   }
#
#   P.M <- exp(tmp.sum-max(tmp.sum) -log(sum(exp(tmp.sum-max(tmp.sum)))) ) + P.M
#
#   if (icounter==50*round(icounter/50,0)) { cat("\n Iteration = ", icounter) }
# }
# t2<-Sys.time()
# t2-t1
#-------------------------------------------------------------------------------

#parallel implementation of the RJMCMC algorithm
#here we attempt to estimate the transition prob matrix

t1.RJMCMC <- Sys.time()

#the maximum number of coefficients in detection and occupancy process
alpha_max <- 4
beta_max <- 4
ncores <- 1
TransMatrix <- matrix(0, nrow=32, ncol=32)

#the data set of detections
K2 <- matrix(K[Anysberg.rows], nrow=156, ncol=35)

require(foreach)
require(doParallel)

cl <- parallel::makeCluster(15)
doParallel::registerDoParallel(cl)

t1<-Sys.time()
par(mfrow=c(6,6), mar=c(1,0,2,0))

#here I calculate the transition probability matrix for moving from one model to the next
for (iMod in 1:32){
  t1mod<-Sys.time()

  #the betas comes from the 'tmp' object
  #this comes from modindex
  if (exists("tmp")==TRUE){
    rm("tmp")
  }else{
    paste(" 'tmp' does not exist. Nothing is deleted.")
  }

  #load the conditioned model
  mod_coefs <- modindex[iMod]
  load(paste("file_",mod_coefs,".RData", sep=""))
  colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
  colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))

  #each row of the TransMatrix uses a different random seed
  set.seed(iMod)
  k.rand.seq <- sample(2:6666, 6665)

  tmp.sum <- rep(0, 32)

  #do this in parallel
  P.M <- foreach(idx=k.rand.seq, .combine=rbind, .inorder = TRUE,.packages=c("mvnfast"))%dopar% {

    #iterate over all of the models
    for (M in 1:32){

      tmp.sum[M] <- M_l_contr(Y=Anysberg.Data,
                            alpha=tmp$alpha[,,idx], beta=tmp$beta[,,idx],
                            tau_alpha=tmp$tau_alpha[idx], tau_beta=tmp$tau_beta[idx],
                            Xdet=det.covs2, Xocc=occ.covs2,
                            det.mat=det.mat, occ.mat=occ.mat,
                            alpha_max=4, beta_max=4,
                            M=M,
                            #mu_alpha...
                            tmp$mu_alpha[,idx],
                            tmp$mu_beta[,idx],
                            const,
                            ncores=1)
    }#end M
    exp(tmp.sum-max(tmp.sum) -log(sum(exp(tmp.sum-max(tmp.sum)))) )
  }#end idx

  TransMatrix[iMod, ] <- apply(P.M, 2, mean)
  t2mod<-Sys.time()

  cat("\n iMod = ", iMod, "; Time = ", t2mod - t1mod )
  plot(TransMatrix[iMod,], type="b", main=iMod, ylim=c(0,1), xaxt="n", yaxt="n")
}
t2<-Sys.time()
t2-t1

stopCluster(cl)

# #See https://youtu.be/6tiryQpHYSk
#Calculate the stationary distribution!
Q <- t(TransMatrix) - diag(1, nrow=32)
Q[1,] <- 1
b <- matrix(0, nrow=32)
b[1] <- 1
probs = solve(Q, b)
cat("\n Posterior model probs = ")
print( cbind(1:32,round(probs,4)) )

model.order <- order(round(probs,4), decreasing =TRUE)
round(probs,4)[model.order]
model.order

t2.RJMCMC <- Sys.time()
RJMCMC.timer <- t2.RJMCMC - t1.RJMCMC
RJMCMC.timer
#-------------------------------------------------------------------------------

#Lets have a look at the top models

#load the conditioned model
#post prob for model 6 = 0.85
#post prob for model 8 = 0.15
mod_coefs <- modindex[6]
load(paste("file_",mod_coefs,".RData", sep=""))
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))

posttables <- function(tmp){

  t1 <- cbind( apply(tmp$mu_alpha, 1, mean),
              apply(tmp$mu_alpha, 1, sd),
              t(apply(tmp$mu_alpha, 1, function(x) quantile(x, probs=c(0.025, 0.975)))) )

  t2 <- cbind( apply(tmp$mu_beta, 1, mean),
               apply(tmp$mu_beta, 1, sd),
               t(apply(tmp$mu_beta, 1, function(x) quantile(x, probs=c(0.025, 0.975)))) )

  rbind(t1,t2)
}

require(xtable)
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))

xtable( posttables(tmp) )

#-------------------------------------------------------------------------------

#Plots of the transition prob matrix
# require(lattice)
# breaks <- seq(0, 1,length.out=200)
# cbPalette <- brewer.pal(n=9, name="YlOrRd")
# Colours<-colorRampPalette(cbPalette)
# levelplot(TransMatrix, xlab="From", ylab="To", col.regions=Colours)
#
# require(rgl)
# plot3d(x=1:32, y=1:32, z=TransMatrix, type="h")
#
#
#
#
# par(mfrow=c(6,6), mar=c(1,0,2,0))
# for (iMod in 1:32){
#   plot(TransMatrix[iMod,], type="b", main=iMod, ylim=c(0,1), xaxt="n", yaxt="n")
# }
#
# #t2<-Sys.time()
# #t2-t1
#-------------------------------------------------------------------------------
# 9. End. RJMCMC code
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 10. Begin Post RJMCMC analysis
#-------------------------------------------------------------------------------

#Lets plot the map of the occupancy probability.
#The RJMCMC indicated that model 8 is best!

#MAPS USING STANDARDISED VARIABLES
longitude = c(scale(station.cov$Longitude[index1], T,T))
latitude = c(scale(station.cov$Latitude[index1],T,T))

x<-longitude
y<-latitude
survlocs <- cbind(x,y)

# #Obtain the basic x grid and y space distance
xspace<- 0.02/1 #xy[53,2]-xy[52,2]= 0.0179 #(consecutive x points)
yspace<- 0.02/1 #xy[5,3] - xy[10,3] = 0.0182 #(consecutive y points)
# #create a grid of points

#makes every possible combination of coordinates from this grid.list:
candidates<- make.surface.grid( list( seq( min(x)-3*xspace,max(x)+3*xspace,0.05 ),
                                      seq(min(y)-3*yspace,max(y)+3*yspace,0.05)))

plot(x, y, pch=".", xlab="Longitude (standardized)", ylab="Latitude (standardized)", type="n")
#points(candidates[,1], candidates[,2], pch=19, cex=.1)
points(x,y,pch=3,col="red", cex=.75)

#approximate border
bdr <- ahull(survlocs[,1], survlocs[,2], alpha = .25)
bdrpts <- bdr$arcs[,7]

for (i in 1:(length(bdrpts)-1)){
  segments(x[bdrpts[i]], y[bdrpts[i]], x[bdrpts[i+1]], y[bdrpts[i+1]], lty=3)
}
segments(x[bdrpts[i+1]], y[bdrpts[i+1]], x[bdrpts[1]], y[bdrpts[1]], lty=3)

#--------------------------------------------------------------

#plot of the maps
cbPalette <- viridis_pal(direction=-1)(100)
breaks <- seq(0, 1,length.out=100)

#load the results for model 6
load("file_6.RData")
tmp6 <- tmp
rm(tmp)

mod_coefs <- modindex[6]
#load(paste("file_",mod_coefs,".RData", sep=""))
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))

#---------------------------------

load("file_8.RData")
tmp8 <- tmp
rm(tmp)

mod_coefs <- modindex[8]
#load(paste("file_",mod_coefs,".RData", sep=""))
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))
#---------------------------------

load("file_22.RData")
tmp22 <- tmp
rm(tmp)

mod_coefs <- modindex[22]
#load(paste("file_",mod_coefs,".RData", sep=""))
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))
#---------------------------------

load("file_24.RData")
tmp24 <- tmp
rm(tmp)

mod_coefs <- modindex[24]
#load(paste("file_",mod_coefs,".RData", sep=""))
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))
#---------------------------------

load("file_5.RData")
tmp5 <- tmp
rm(tmp)

mod_coefs <- modindex[5]
#load(paste("file_",mod_coefs,".RData", sep=""))
colnames(model.matrix(models_occ_list_names[mod_coefs][[1]], occ.covs))
colnames(model.matrix(models_det_list_names[mod_coefs][[1]], det.covs))
#---------------------------------

#Posterior means using BMA
#.8543
#.1442
#0.0011
#0.0002
#0.0001

#detection process
#trail
tm_BMAd <- apply(tmp6$mu_alpha, 1, mean)[2]*.8543 + apply(tmp8$mu_alpha, 1, mean)[2]*.1442 +
apply(tmp22$mu_alpha, 1, mean)[2]*.0011 + apply(tmp24$mu_alpha, 1, mean)[2]*.0002

#plain
pm_BMAd <- apply(tmp8$mu_alpha, 1, mean)[3]*.1442 + apply(tmp24$mu_alpha, 1, mean)[3]*.0002

#river
rm_BMAd <- apply(tmp8$mu_alpha, 1, mean)[4]*.1442 + apply(tmp24$mu_alpha, 1, mean)[4]*.0002

#---------------------------------

#occupancy process
#elevation
em_BMAo <- apply(tmp6$mu_beta, 1, mean)[2]*.8543 + apply(tmp8$mu_beta, 1, mean)[2]*.1442 +
  apply(tmp22$mu_beta, 1, mean)[2]*.0011 + apply(tmp24$mu_beta, 1, mean)[2]*.0002 +
  apply(tmp5$mu_beta, 1, mean)[2]*.0002

#preyindex
pm_BMAo <- apply(tmp22$mu_beta, 1, mean)[3]*.0011 + apply(tmp24$mu_beta, 1, mean)[3]*.0002

#---------------------------------

#calculate the variances using BMA
#detection process
#trail
tv_BMAd <- (apply(tmp6$mu_alpha, 1, var)[2] + apply(tmp6$mu_alpha, 1, mean)[2]^2)*.8543 + (apply(tmp8$mu_alpha, 1, var)[2] + apply(tmp8$mu_alpha, 1, mean)[2]^2)*.1442 +
  (apply(tmp22$mu_alpha, 1, var)[2] + apply(tmp22$mu_alpha, 1, mean)[2]^2)*.0011 + (apply(tmp24$mu_alpha, 1, var)[2] + apply(tmp24$mu_alpha, 1, mean)[2])*.0002 - tm_BMAd^2

#plain
pv_BMAd <- (apply(tmp8$mu_alpha, 1, var)[3] + apply(tmp8$mu_alpha, 1, mean)[3]^2)*.1442 + (apply(tmp24$mu_alpha, 1, var)[3] + apply(tmp24$mu_alpha, 1, mean)[3]^2)*.0002 - pm_BMAd^2

#river
rv_BMAd <- (apply(tmp8$mu_alpha, 1, var)[4] + apply(tmp8$mu_alpha, 1, mean)[4]^2)*.1442 + (apply(tmp24$mu_alpha, 1, var)[4] + apply(tmp24$mu_alpha, 1, mean)[4]^2)*.0002 - rm_BMAd^2

#---------------------------------

#detection process
#elevation
ev_BMAo <- (apply(tmp6$mu_beta, 1, var)[2] + apply(tmp6$mu_beta, 1, mean)[2]^2)*.8543 + (apply(tmp8$mu_beta, 1, var)[2] + apply(tmp8$mu_beta, 1, mean)[2]^2)*.1442 +
  (apply(tmp22$mu_beta, 1, var)[2] + apply(tmp22$mu_beta, 1, mean)[2]^2)*.0011 + (apply(tmp24$mu_beta, 1, var)[2] + apply(tmp24$mu_beta, 1, mean)[2]^2)*.0002 +
  (apply(tmp5$mu_beta, 1, var)[2] + apply(tmp5$mu_beta, 1, mean)[2]^2)*.0002 - em_BMAo^2

#preyindex
pv_BMAo <- (apply(tmp22$mu_beta, 1, var)[3] + apply(tmp22$mu_beta, 1, mean)[3]^2)*.0011 + (apply(tmp24$mu_beta, 1, var)[3] + apply(tmp24$mu_beta, 1, mean)[3]^2)*.0002 - pm_BMAo^2

#---------------------------------

xtable(
rbind(
  cbind(c(tm_BMAd, pm_BMAd, rm_BMAd), c(tv_BMAd, pv_BMAd, rv_BMAd)^.5, c(1-.0001, 1-.8543-.0011, 1-.8543-.0011)),
cbind(c(em_BMAo, pm_BMAo), c(ev_BMAo, pv_BMAo)^.5, c(1, 1-.8543-.1442-.0001))
))
#---------------------------------

par(mfrow=c(6,6), pty="s", mai=c(0,0,0,0))
#conditional occupancy maps
for (ii in 1:35){
  tmpz <- tmp6$z[,ii,]*.8543 + tmp8$z[,ii,]*.1442 + tmp22$z[,ii,]*.0011 +tmp24$z[,ii,]*.0002 +tmp5$z[,ii,]*.0001
  mapv<-apply(tmpz, 1, mean)
  make_map2(mapv, 100, cbind(longitude, latitude))
  legend("bottomright", legend=paste(round(mean(mapv),2)), bty="n")
  legend("topleft", legend=paste(ii), bty="n")
  pp<-which(Anysberg.Data[ii,]>0)
  points(longitude[pp], latitude[pp], pch=20, col="red")
}

is<-image.scale2(mapv, col=cbPalette[-1],
                 breaks=breaks, horiz=F, xlab="", ylab="",
                 side=2, cex.axis=1)
legend("topright", legend="Detection", pch=20, col="red",
       bty="n")
#---------------------------------

par(mfrow=c(6,6), pty="s", mai=c(0,0,0,0))
#conditional detection prob maps

#md<-NULL
for (ii in 1:35){
  tmpd <- tmp6$p_array[,ii,]*.8543 + tmp8$p_array[,ii,]*.1442 + tmp22$p_array[,ii,]*.0011 +tmp24$p_array[,ii,]*.0002 +tmp5$p_array[,ii,]*.0001
  mapv<-apply(tmpd, 1, mean)
  #md[ii]<-mean(mapv)
  make_map2(mapv, 100, cbind(longitude, latitude))
  legend("bottomright", legend=paste(round(mean(mapv),2)), bty="n")
  legend("topleft", legend=paste(ii), bty="n")
  pp<-which(Anysberg.Data[ii,]>0)
  points(longitude[pp], latitude[pp], pch=20, cex=.35, col="red")
}
is<-image.scale2(mapv, col=cbPalette[-1],
                 breaks=breaks, horiz=F, xlab="", ylab="",
                 side=2, cex.axis=1)

legend("topright", legend="Detection", pch=15, col="red",
       bty="n")
#---------------------------------

#maps for the different trophic guilds
#below I plot the occupancy probs associated with the different guilds
#setwd("~/Documents/AllanClark/Stats/Research/Drouilly_M/Analysis/Dorazio_occupancy_Models/RequiredData")
species.cov <- read.csv("animal_covariates.csv", header=TRUE,sep=";")

#setwd("/Users/allanclark/Documents/AllanClark/Stats/Research/Clark_AE/PhD/Paper4/AnalysisExamples/Example1/Code")

species.cov.index<-1:dim(species.cov)[1]

species.cov.index.carnivore<-species.cov.index[as.character(species.cov$Trophic.category)=="carnivore"]
num.carnivore<-length(species.cov.index.carnivore)
species.cov.index.carnivore
#8  9 12 16 18 19 24 26 37 38 39 are the carnivores in the covariates data set
#they are not all observed at Anysberg
#8  9 12 18 19 24 26 37 38 39 are observed at Anysberg
#the particular species that are carnivores as found in the species.cov object
#as.character(species.cov$Species[species.cov.index.carnivore])

species.cov.index.carnivore<-species.cov.index.carnivore[as.character(species.cov$Species[species.cov.index.carnivore]) %in%Anysberg.names]
num.obs.carnivore<-length(species.cov.index.carnivore)
#the particular carnovores oberved at Anysberg!
carnivore.names<-as.character(species.cov$Species[species.cov.index.carnivore])

#this indicates the associated index
carnivore.index<-NULL
for (i in 1:length(carnivore.names)){
  carnivore.index[i]<-which(Anysberg.names==carnivore.names[i])
}
carnivore.index

# par(mfrow=c(2,2), pty="s")
# #mapv<-apply(occ.jags.M4$sims.list$Z[ , , carnivore.index], 2, mean)
# mapv<-apply(temp.spat2$z[ , carnivore.index, ], 1, mean)
# make_map(mapv, 200, cbind(longitude, latitude))
# mean.carnivore<-mean(mapv)
# #apply(occ.jags.M4$sims.list$Z[ , , c(11,18)], 2, mean) -apply(rbind(occ.jags.M4$sims.list$Z[ , , 11], occ.jags.M4$sims.list$Z[ , , 18]), 2, mean)

#--------------------------------------

#identify which species are herbivores, herbivores, omnivores and insectivores
species.cov.index.herbivore<-species.cov.index[as.character(species.cov$Trophic.category)=="herbivore"]
num.herbivore<-length(species.cov.index.herbivore)

#which of the herbivores are observed at Anysberg
species.cov.index.herbivore<-species.cov.index.herbivore[as.character(species.cov$Species[species.cov.index.herbivore]) %in%Anysberg.names]
num.obs.herbivore<-length(species.cov.index.herbivore)
herbivore.names<-as.character(species.cov$Species[species.cov.index.herbivore])

#this indicates the associated index
herbivore.index<-NULL
for (i in 1:length(herbivore.names)){
  herbivore.index[i]<-which(Anysberg.names==herbivore.names[i])
}

# mapv<-apply(temp.spat2$z[ , herbivore.index, ], 1, mean)
# #mapv<-apply(occ.jags.M4$sims.list$Z[ , , herbivore.index], 2, mean)
# make_map(mapv, 200, cbind(longitude, latitude))
# mean.herbivore<-mean(mapv)
#--------------------------------------

#identify which species are omnivores, omnivores, omnivores and insectivores
species.cov.index.omnivore<-species.cov.index[as.character(species.cov$Trophic.category)=="omnivore"]
num.omnivore<-length(species.cov.index.omnivore)

#which of the omnivores are oberved at Anysberg
species.cov.index.omnivore<-species.cov.index.omnivore[as.character(species.cov$Species[species.cov.index.omnivore]) %in%Anysberg.names]
num.obs.omnivore<-length(species.cov.index.omnivore)
omnivore.names<-as.character(species.cov$Species[species.cov.index.omnivore])

#this indicates the associated index
omnivore.index<-NULL
for (i in 1:length(omnivore.names)){
  omnivore.index[i]<-which(Anysberg.names==omnivore.names[i])
}

# mapv<-apply(temp.spat2$z[ , omnivore.index, ], 1, mean)
# #mapv<-apply(occ.jags.M4$sims.list$Z[ , , omnivore.index], 2, mean)
# make_map(mapv, 200, cbind(longitude, latitude))
# mean.omnivore<-mean(mapv)
#--------------------------------------

#identify which species are insectivores, insectivores, insectivores and insectivores
species.cov.index.insectivore<-species.cov.index[as.character(species.cov$Trophic.category)=="insectivore"]
num.insectivore<-length(species.cov.index.insectivore)

#which of the insectivores are oberved at Anysberg
species.cov.index.insectivore<-species.cov.index.insectivore[as.character(species.cov$Species[species.cov.index.insectivore]) %in%Anysberg.names]
num.obs.insectivore<-length(species.cov.index.insectivore)
insectivore.names<-as.character(species.cov$Species[species.cov.index.insectivore])

#this indicates the associated index
insectivore.index<-NULL
for (i in 1:length(insectivore.names)){
  insectivore.index[i]<-which(Anysberg.names==insectivore.names[i])
}

# mapv<-apply(temp.spat2$z[ , insectivore.index, ], 1, mean)
# #mapv<-apply(occ.jags.M4$sims.list$Z[ , , insectivore.index], 2, mean)
# make_map(mapv, 200, cbind(longitude, latitude))
# mean.insectivore<-mean(mapv)
#--------------------------------------

guildindex <- c(carnivore.index, herbivore.index, insectivore.index, omnivore.index, c(9, 10, 30, 32 ))
guildcols<- rep("brown", length(guildindex))
guildcols[carnivore.index]<-"red"
guildcols[herbivore.index]<-"forestgreen"
guildcols[insectivore.index]<-"black"
guildcols[omnivore.index]<-"blue"
#guildcols[-c(1:31)]<-"brown"

#note the order!!!!
Anysberg.names[guildindex]
Anysberg.names.common<-rep(1,35)
Anysberg.names.common[1]="Wildcat"
Anysberg.names.common[2]="Badger"
Anysberg.names.common[3]="Leopard"
Anysberg.names.common[4]="Otter"
Anysberg.names.common[5]="Marsh Mongoose"
Anysberg.names.common[6]="Grey Mongoose"
Anysberg.names.common[7]="Polecat"
Anysberg.names.common[8]="Jackal"
Anysberg.names.common[9]="Caracal"
Anysberg.names.common[10]="Hyaena"
Anysberg.names.common[11]="Oryx"
Anysberg.names.common[12]="Eland"
Anysberg.names.common[13]="Kudu"
Anysberg.names.common[14]="Red Hartebeest"
Anysberg.names.common[15]="Zebra"
Anysberg.names.common[16]="Hare"
Anysberg.names.common[17]="Klipspringer"
Anysberg.names.common[18]="Grey Rhebok"
Anysberg.names.common[19]="Hyrax"
Anysberg.names.common[20]="Smith's rock hare"
Anysberg.names.common[21]="Steenbok"
Anysberg.names.common[22]="Duiker"
Anysberg.names.common[23]="Sprinbok"
Anysberg.names.common[24]="Aardvark"
Anysberg.names.common[25]="Aardwolf"
Anysberg.names.common[26]="Bat-eared fox"
Anysberg.names.common[27]="Vervet Monkey"
Anysberg.names.common[28]="Genet"
Anysberg.names.common[29]="Porcupine"
Anysberg.names.common[30]="Chacma baboon"
Anysberg.names.common[31]="Bushpig"
Anysberg.names.common[32]="Black Korhaan"
Anysberg.names.common[33]="Karoo Korhaan"
Anysberg.names.common[34]="Spurfowl"
Anysberg.names.common[35]="Francolin"

print(xtable(cbind(Anysberg.names[guildindex], Anysberg.names.common)), include.rownames=FALSE)
#--------------------------------------

#This is the Bayesian model averaged map!
pdf("AnysbergOccupancyMapBMA.pdf")
par(mfrow=c(6,6), pty="s", mai=c(0,0,0,0), mar=c(0.25,0,2,0))
counter<-0

for (ii in guildindex){
  counter<- 1 + counter
  tmpz <- tmp6$z[,ii,]*.8543 + tmp8$z[,ii,]*.1442 + tmp22$z[,ii,]*.0011 +tmp24$z[,ii,]*.0002 +tmp5$z[,ii,]*.0001
  mapv<- apply(tmpz, 1, mean)
  make_map2(mapv, 100, cbind(longitude, latitude))
  cat("\n ii=", ii, " ", fivenum(mapv))

  legend("bottomright", legend=paste("(",round(mean(mapv),2),")", sep=""), bty="n", cex=1.1)
  #legend("topleft", legend=paste(ii), bty="n", text.col=guildcols[ii])
  pp<-which(Anysberg.Data[ii,]>0)
  points(longitude[pp], latitude[pp], pch=20, cex=.2, col="red")
  title(main=Anysberg.names.common[counter], col.main=guildcols[ii], cex.main=1.15)
  #print(guildcols[ii])
}

is<-image.scale2(mapv, col=cbPalette[-1],
                 breaks=breaks, horiz=F, xlab="", ylab="",
                 side=2, cex.axis=1)
legend("topright", legend="Detection", pch=20, col="red",
       bty="n", cex=1)
dev.off()
#--------------------------------------

#plot the regression effects

pred.ints<-function(x,y, numspecies, c.lab=.9)
{
  x.q1<-c(quantile(x, c(0.025, 0.975)))
  x.q2<-t(apply(y, 2, function(x) c(quantile(x, c(0.025, 0.975)))))

  min1<-min(x.q1, x.q2)
  max1<-max(x.q1, x.q2)

  plot(1:numspecies, type="n", xlim=c(min1, max1), yaxt='n', xlab="Estimates",
       cex.lab=c.lab, ylab="", xaxt="n")
  axis(side=1, cex.axis=c.lab)

  polygon(c(x.q1[1],x.q1[1], x.q1[2], x.q1[2]), c(0, numspecies+1, numspecies+1, 0),
          col = "lightgrey", border = NA)
  abline(v=0, lty=2)
  abline(v=mean(x), col="red", lwd=1.5)
  abline(v=x.q1[1], lty=2, col="red", lwd=1.5)
  abline(v=x.q1[2], lty=2, col="red", lwd=1.5)

  for(i in 1:numspecies)
  {
    if (x.q2[i,2]< 0){
      segments(x.q2[i,1], i, x.q2[i,2], i, col="blue", lwd=2)
      points( mean(y[,i]), i, pch=19)
    }else
      if (x.q2[i,1]> 0){
        segments(x.q2[i,1], i, x.q2[i,2], i, col="blue", lwd=2)
        points(mean(y[,i]),i, pch=19)
      }else
      {
        segments(x.q2[i,1], i, x.q2[i,2], i, col="black")
        points(mean(y[,i]), i, pch=19, col="black")
      }
  }

  box(col="black")
}

trailBMA <- t(tmp6$alpha[2,,]*.8543 + tmp8$alpha[2,,]*.1442 + tmp22$alpha[2,,]*.0011 +tmp24$alpha[2,,]*.0002)
trailBMA <- trailBMA[, guildindex]

elevBMA <- t(tmp6$beta[2,,]*.8543 + tmp8$beta[2,,]*.1442 + tmp22$beta[2,,]*.0011 +tmp24$beta[2,,]*.0002 + tmp5$beta[2,,]*.0001)
elevBMA <- elevBMA[, guildindex]

numspecies<-dim(Anysberg.Data)[1] #35

pdf("Anysberg_covariateeffects.pdf")
  par(mar=c(3,.1,.5,0.1), mai=c(.8,2,.35,0.1))
  ll <- layout(matrix(1:2, nrow=1, byrow=TRUE),
               widths=c(9,5))

  xx1 <- tmp6$mu_alpha[2,]*.8543 + tmp8$mu_alpha[2,]*.1442 + tmp22$mu_alpha[2,]*.0011 +tmp24$mu_alpha[2,]*.0002
  pred.ints(x= xx1 , y=trailBMA, numspecies, .8)

  for (i in 1:35){
    axis(2, at=i, las=1, cex.axis=0.8,
         col.axis = guildcols[guildindex[i]],
         labels=Anysberg.names.common[i])
  }
  title("(a)", cex.main=1)

  par(mar=c(3,.1,.5,0.1), mai=c(.8,0,.35,0.1))
  xx2 <- tmp6$mu_beta[2,]*.8543 + tmp8$mu_beta[2,]*.1442 + tmp22$mu_beta[2,]*.0011 +tmp24$mu_beta[2,]*.0002 + tmp5$mu_beta[2,]*.0001
  pred.ints(x= xx2 , y=elevBMA, numspecies, .8)

  title("(b)", cex.main=1)
dev.off()

#-------------------------------------------------------------------------------
# 10. End Post RJMCMC analysis
#-------------------------------------------------------------------------------

save.image("Section10.RData")


