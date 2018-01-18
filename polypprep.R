#Suzan Shahrestani
#Date: 08132017
#Polyp Distribution Analysis
#Presence/Absence
#52 sites
#Explore 3 variables: residence time july/winter and variability
#GLM model selection
#Error Analysis
#Model Predictions : Hotspots for jellyfish polyps

rm(list = ls())#clear ws
options(warn=-1)#surpress warnings
#Get Data
setwd('/polyp/')
load(file = 'data.RData')

#install.packages("dplyr", "jpeg", "raster", "mgcv", "boot", "resample")
#Required Pagackages
library(dplyr)
library(jpeg)
library(raster)
library(mgcv)
library(boot)
library(resample)

################################################################################
#################    Rasterising Residence - Variable Images            ########
################################################################################
#Extracting vales from jpeg images based on color
#JULY RESIDENCE TIME
rasjuly<-raster()
extent(rasjuly)<-extent(res_july)
dim(rasjuly)<-dim(res_july)
rasjuly[]<-getValues(res_july)
rasjuly[(rasjuly==5)]<-300
rasjuly[(rasjuly==1)]<-40
rasjuly[(rasjuly==0)]<-NA
rasjuly[(rasjuly==6)]<-280
rasjuly[(rasjuly==7)]<-40
rasjuly[(rasjuly==2)]<-80
rasjuly[(rasjuly==4)]<-200
rasjuly[(rasjuly==8)]<-160
rasjuly[(rasjuly==3)]<-120
rasjuly[(rasjuly==9)]<-20

#VARIABILITY RESIDENCE TIME
rasvar<-raster()
extent(rasvar)<-extent(res_var)
dim(rasvar)<-dim(res_var)
rasvar[]<-getValues(res_var)
rasvar[(rasvar==0)]<-NA
rasvar[(rasvar==1)]<-0
rasvar[(rasvar==2)]<-10
rasvar[(rasvar==3)]<-20
rasvar[(rasvar==4)]<-30
rasvar[(rasvar==5)]<-50
rasvar[(rasvar==6)]<-40
rasvar[(rasvar==7)]<-NA
rasvar[(rasvar==8)]<-NA

#JANUARY RESIDENCE TIME
rasjan<-rasjuly-rasvar

#SALINITY
ras_sal<-crop(ras_sal, extent(res_july))
rassal<-raster(ras_sal)
extent(rassal)<-extent(ras_sal)
dim(rassal)<-dim(ras_sal)
rassal[]<-getValues(ras_sal)
rassal[(rassal==12)]<-NA
rassal[(rassal==9)]<-18.1
rassal[(rassal==1)]<-0
rassal[(rassal==2)]<-0.6
rassal[(rassal==7)]<-12.6
rassal[(rassal==6)]<-10.1
rassal[(rassal==4)]<-5.1
rassal[(rassal==3)]<-2.6
rassal[(rassal==11)]<-24.1
rassal[(rassal==5)]<-7.6
rassal[(rassal==8)]<-15.1
rassal[(rassal==10)]<-21.1

################################################################################
###################         Georeferencing the images           ################
################################################################################
#Adding a coordinate system to the images
#JAN RESIDENCE
#Image lat and lon extent (bounding box)
lonextent = 77.366809-75.501253
latextent = 39.594411-36.789872
#Setting the resolution(pixel units)
resx = lonextent/extent(res_var)[2]
resy = latextent/extent(res_var)[4]
#carrying over values
rasjanv<-rasjan
#modifying  extent and resolution
extent(rasjan) <- c(-77.366809,-75.501253, 39.594411,36.789872)
res(rasjan)<-c(resx,resy)
#***Check that the dimensions are the same as the png dimensions
#changing the projection
projection(rasjan) <- CRS("+proj=longlat +datum=WGS84")
#adding values back in
rasjan[]<-values(rasjanv)

#JULY RESIDENCE
rasjulyv<-rasjuly
rasjuly<-raster()
extent(rasjuly)<-extent(rasjan)
dim(rasjuly)<-dim(rasjan)
res(rasjuly)<-res(rasjan)
rasjuly[]<-values(rasjulyv)

#VARIABILITY RESIDENCE
rasjvarv<-rasvar
rasvar<-raster()
extent(rasvar)<-extent(rasjuly)
dim(rasvar)<-dim(rasjuly)
res(rasvar)<-res(rasjuly)
rasvar[]<-values(rasjvarv)

#SALINITY
rassalv<-rassal
rassal<-raster()
extent(rassal)<-extent(rasjuly)
dim(rassal)<-dim(rasjuly)
res(rassal)<-res(rasjuly)
rassal[]<-values(rassalv)

################################################################################
#######	Extracting values
################################################################################
#Extracting information from  images, with regard to 52 site locations
#use_ df all
df_all$res_jul<-extract(rasjuly, cbind(df_all$lon, df_all$lat))
df_all$res_jan<-extract(rasjan, cbind(df_all$lon, df_all$lat))
df_all$res_var<-extract(rasvar, cbind(df_all$lon, df_all$lat))
df_all$res_sal<-extract(rassal, cbind(df_all$lon, df_all$lat))

################################################################################
#########		GLM 		                                
################################################################################
#GLM for presence/absence  data
##Model Selection - with Cargo counts and Settlement towers
#variables: waterscape, shore, bayarea, res_jul, res_jan, resvar
fit1 <- glm(polyps ~  res_jan,
		family = binomial(link = "probit"),  data = df_all)
fit2 <- glm(polyps ~  res_jul,
		family = binomial(link = "probit"),  data = df_all)
fit3 <- glm(polyps ~  res_jul + res_jan,
		family = binomial(link = "probit"),  data = df_all)
fit4 <- glm(polyps ~  res_jul*res_jan,
		family = binomial(link = "probit"),  data = df_all)

summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)

c(
AIC(fit1),
AIC(fit2),
AIC(fit3),
AIC(fit4)
)

######################CROSS - VALIDATION
set.seed(112086)
# leave-one-out and 10-fold cross-validation prediction error for
# the nodal data set.  Since the response is a binary variable an
# appropriate cost function is
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
(cv.err <- cv.glm(df_all, fit1, cost, K =10)$delta)
(cv.err <- cv.glm(df_all, fit2, cost, K =10)$delta)
(cv.err <- cv.glm(df_all, fit3, cost, K =10)$delta)
(cv.err <- cv.glm(df_all, fit4, cost, K =10)$delta)

#*************Computationally Demanding
# cv1<-NULL
# for(i in 1:1000){
# cv1[i]<-cv.glm(df_all, fit1, cost, K =10)$delta
# }
#
# cv2<-NULL
# for(i in 1:1000){
# cv2[i]<-cv.glm(df_all, fit2, cost, K =10)$delta
# }
#
# cv3<-NULL
# for(i in 1:1000){
# cv3[i]<-cv.glm(df_all, fit3, cost, K =10)$delta
# }
#
# cv4<-NULL
# for(i in 1:1000){
# cv3[i]<-cv.glm(df_all, fit4, cost, K =10)$delta
# }
#
# mean(cv1)
# mean(cv2)
# mean(cv3)
# mean(cv4)

######################GOODNESS OF FIT - Preidction Error
#Jack-Knife estimations
glmfunc1<-function(df){
	mean(glm(polyps ~  res_jan ,
		family = binomial(link = "probit"),  data = df)$fitted.values)}
glmfunc2<-function(df){
	mean(glm(polyps ~  res_jul ,
		family = binomial(link = "probit"),  data = df)$fitted.values)}
glmfunc3<-function(df){
	mean(glm(polyps ~  res_jul*res_jan ,
		family = binomial(link = "probit"),  data = df)$fitted.values)}

jackdf<-with(df_all, data.frame(polyps, res_jul, res_jan))
####
args<-as.list(c('polyps', 'res_july', 'res_jan'))
jackknife(jackdf, glmfunc1, trace=FALSE)
jackknife(jackdf, glmfunc2, trace=FALSE)
jackknife(jackdf, glmfunc3, trace=FALSE)

################################################################################
###############################    SPATIAL PREDICTION   ########################
################################################################################
#Predict polyp presence in all of Chesapeake Bay
stras<-stack(rasjan, rasjuly)
names(stras)<-c('res_jan', 'res_jul')
pred <- predict(stras, se.fit=TRUE, fit3, type='response')
raspolyp<-raster()
extent(raspolyp)<-extent(stras)
dim(raspolyp)<-dim(stras)
raspolyp[]<-values(pred)
################################################################################
#######################       Applying Salinity Threshold       ################
################################################################################
raspolyp[(rassal<= 2 )] <- 0
raspolypv<-raspolyp
################################################################################
#######################      PLOTTING       ####################################
################################################################################
breakpoints <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
colors<-c("#3288bd","#99d594","#e6f598",'#ffffbf',"#fee08b","#fc8d59","#d53e4f")
plot(raspolyp,breaks=breakpoints,col=colors, colNA='darkgray')

dfc<-filter(df_all, scientist=='cargo')
dfa<-filter(dfc, polyps==0)
dfp<-filter(dfc, polyps==1)
points(dfp$lon, dfp$lat, cex=1.5, pch=19, col='black')
points(dfa$lon, dfa$lat, cex=1.5, pch=1, col='black')
points(dfa$lon, dfa$lat, cex=1.5, pch=20, col='white')

#Suzan
dfs<-filter(df_all, scientist=='suzan')
dfa<-filter(dfs, polyps==0)
dfp<-filter(dfs, polyps==1)
points(dfp$lon, dfp$lat, cex=1.5, pch=18, col='black')
points(dfa$lon, dfa$lat, cex=1.5, pch=5, col='black')
points(dfa$lon, dfa$lat, cex=1.5, pch=18, col='white')
