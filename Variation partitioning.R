#-----------------------------------------------------------------------------------------------
# Script: gam for multivariate data: example
# Project: METACOM 1
# Date: 25/03/2021
#-----------------------------------------------------------------------------------------------

library(vegan)
library(mgcv)
library(gllvm)
library(doParallel)
library(foreach)
library(spdep)
library(adespatial)

source("Gam_functions.R")

# load data
xy <- read.csv2("XY.csv", header = TRUE, row.names = 1, sep = ";") #XY matrix
env <- read.csv2("env.csv", header = TRUE, row.names = 1, sep = ";") #environmental matrix
sp <- read.csv2("species.csv", header = TRUE, row.names = 1, sep = ";") #Species matrix 
#------------------------------------------------------------------------------------------------
# Arrange data
sp <- sp[,colSums(sp) != 0]
sp<-decostand(sp, method="pa")
gll.hat <- gll.hat.fun(sp = sp, family = "binomial", max.lv = 3) # latent variables

# PCA environment
pca.env <- as.data.frame(scores(rda(env, scale = TRUE), choices = c(1:3), display = "si")) #3 first PCs

#-------------------------------------------------------------------------------------------------
# Forward selection

## environment
r2.all <- gam.multi.par(sp = gll.hat, pred = pca.env) #R2 with all variables
r2a.thresh <- 1 - (1 - r2.all) * (nrow(sp) - 1)/(nrow(sp) - ncol(pca.env) - 1) #Adjusted R2 with all variables (Threshold)
fw <- fw.sel.gam.par(sp = gll.hat, pred = pca.env, r2a.thresh = r2a.thresh) #Forward selection with double-stopping criterion
fw #selected variables
env.sel <- pca.env[, fw$variable] #Selected variables matrix
env.sel <- as.data.frame(env.sel) #As dataframe
R2X1<-gam.multi.par(sp = gll.hat, pred = env.sel)#R2 with the selected variables
R2X1
R2aX1<-adjR2X1(sp=gll.hat,env=env.sel,xy=xy, nrepet=10)#Adjusted R2 using MSR method
R2aX1

# space
r2.all <- gam.multi.par(sp = gll.hat, pred = xy) #R2 with all variables
r2a.thresh <- 1 - (1 - r2.all) * (nrow(sp) - 1)/(nrow(sp) - ncol(xy) - 1) #Adjusted R2 with all variables (Threshold)
fw <- fw.sel.gam.par(sp = gll.hat, pred = xy, r2a.thresh = r2a.thresh) #Forward selection with double-stopping criterion
fw #Forward selection with double-stopping criterion
spa.sel <- xy[, fw$variable] #Selected variables matrix
spa.sel <- as.data.frame(spa.sel) #As data.frame
R2X2<-gam.multi.par(sp = gll.hat, pred = spa.sel)#R2 with the selected variables
R2X2
R2aX2 <- 1 - (nrow(gll.hat) - 1) / (nrow(gll.hat) - ncol(xy) - 1) * (1 - R2X2)#Adjusted R2 using Ezequiel formula

#environment+space
all<-cbind(env.sel,spa.sel) #Combination of environmental and spatial variables
R2X1X2 <- gam.multi.par(sp = gll.hat, pred = all) #R2 with the selected variables
R2X1X2
R2aX1X2<-adjR2X1X2(sp=gll.hat,env.sel,spa=spa.sel,xy=xy, nrepet=10) #Adjusted R2 using MSR method

#Partitioning
a <- 1 - (1-(R2X1X2 - R2X2)) / (1-mean(R2aX1X2 - R2X2)) #PURE ENVIRONMENT
b <- R2aX1 -  a #Environment-Space OVERLAP
c <- R2aX2 - b #PURE SPACE
d <- 1 - (a + b + c) #RESIDUALS

vp <- c(a, b, c, d) #Vector of the fractions
names(vp) <- letters[1:length(vp)]
vp #Results of variation partitioning
