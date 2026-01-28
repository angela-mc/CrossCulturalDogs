
rm(list=ls())
library("ggmcmc")
library(brms)
library(rstan)
library(car)
library(pROC)
library(performance)

##############
##### data ###
##############

# IMPORT positive care, negative care, and personhood from Supplementary Data 1
  # save it in df with society name as column

df[is.na(df$personhood),]$personhood<- 0 # make personhood NA = 0

read.csv(file="Rscripts/toshare/data/SupplementaryData2.csv",stringsAsFactors = F)->restdata
restdata<-restdata[match(df$v0.2.socName, restdata$v0.2.socName),]
restdata$poscare<-df$poscare
restdata$negcare<-df$negcare
restdata$personhood<-df$personhood
rm(df)
df<-restdata
rm(restdata)


#############################
##### space and phylogeny ###
#############################

load("Rscripts/toshare/data/SupplementaryData3.rds")
nametree<-nametree_jeg
nametree<-ladderize(nametree)

df$Language_ID <- df$v0.2.socName
df$Language_ID2 <- df$v0.2.socName

#checking for duplicated locations and jittering them
duplicate_coords = df[duplicated(df[,c("Longitude", "Latitude")]) | duplicated(df[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = df$Language_ID %in% duplicate_coords
df$Latitude[duplicate_rowid] = jitter(df$Latitude[duplicate_rowid], factor = 1)
df$Longitude[duplicate_rowid] = jitter(df$Longitude[duplicate_rowid], factor = 1)

phylo_covar_mat <- ape::vcv(nametree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# phylo_covar_mat->p1
# phylo_covar_mat<- ape::vcv.phylo(nametree, corr=TRUE)
# phylo_covar_mat->p2
kappa = 1 # smoothness parameter as recommended by Dinnage et al. (2020)
phi = c(1, 1) # Sigma parameter. First value is not used.

source("Rscripts/toshare/scripts/varcov.spatial_function.R")
spatial_covar_mat = varcov.spatial(df[,c("Longitude", "Latitude")], cov.pars = phi, kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(df$Language_ID2, df$Language_ID2)
spatial_covar_mat <- spatial_covar_mat / max(spatial_covar_mat)


#################################
##### prepare data for models ###
#################################

source("Rscripts/toshare/scripts/dfprep.R")
colnames(df)
Fcols01<-c("Latitude", "Longitude","Language_ID", "Language_ID2", 
           "FarmingPropensity", "AnimalHusbandry", 
           "MeanAnnualTemperature", "NumberFunctionsIfAllAreNonNA", "NumberParagraphs")
prior <- c(#set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("exponential(1)", class = "sd"))


# POS CARE data
df1<-dfprep(df, "poscare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat,Fcols01=Fcols01)$df1
pmat<-dfprep(df, "poscare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$phylo_covar_mat
spmat<-dfprep(df, "poscare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$spatial_covar_mat
fit_all <- brm(formula=focalcol ~ NumberFunctionsIfAllAreNonNA +MeanAnnualTemperature+ FarmingPropensity +AnimalHusbandry+NumberParagraphs+
                 (1|gr(Language_ID, cov = pmat)) + (1|gr(Language_ID2, cov=spmat)), data = list(df1),
               prior = prior, family = "bernoulli", control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE),
               cores=4, data2 = list(pmat = pmat, spmat = spmat), iter=7000,thin=10) # default is warmup = iter/2
summary(fit_all)
r2_bayes (fit_all)

# NEG CARE data
df1<-dfprep(df, "negcare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat,Fcols01=Fcols01)$df1
pmat<-dfprep(df, "negcare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$phylo_covar_mat
spmat<-dfprep(df, "negcare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$spatial_covar_mat
fit_all <- brm(formula=focalcol ~ NumberFunctionsIfAllAreNonNA +MeanAnnualTemperature+ FarmingPropensity +AnimalHusbandry+NumberParagraphs+
                 (1|gr(Language_ID, cov = pmat)) + (1|gr(Language_ID2, cov=spmat)), data = list(df1),
               prior = prior, family = "bernoulli", control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE),
               cores=4, data2 = list(pmat = pmat, spmat = spmat), iter=7000,thin=10) # default is warmup = iter/2
summary(fit_all)
r2_bayes (fit_all)


# PERSONHOOD data
df1<-dfprep(df, "personhood", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat,Fcols01=Fcols01)$df1
pmat<-dfprep(df, "personhood", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$phylo_covar_mat
spmat<-dfprep(df, "personhood", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$spatial_covar_mat
fit_all <- brm(formula=focalcol ~ NumberFunctionsIfAllAreNonNA +MeanAnnualTemperature+ FarmingPropensity +AnimalHusbandry+NumberParagraphs+
                 (1|gr(Language_ID, cov = pmat)) + (1|gr(Language_ID2, cov=spmat)), data = list(df1),
               prior = prior, family = "bernoulli", control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE),
               cores=4, data2 = list(pmat = pmat, spmat = spmat), iter=7000,thin=10) # default is warmup = iter/2
summary(fit_all)
r2_bayes (fit_all)


#########################
##### function models ###
#########################

Fcols01<-c("Latitude", "Longitude","Language_ID", "Language_ID2", 
           "FarmingPropensity", "AnimalHusbandry", 
           "MeanAnnualTemperature", "NumberParagraphs",
           "v2.1.hunt_01", "v2.2.defense_01", "v2.3.guardHerds_01" , "v2.4.herding_01", "newcarry_01")
colnames(df)


# POS CARE data
df1<-dfprep(df, "poscare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat,Fcols01=Fcols01)$df1
pmat<-dfprep(df, "poscare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$phylo_covar_mat
spmat<-dfprep(df, "poscare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$spatial_covar_mat
fit_all <- brm(formula=focalcol ~ v2.1.hunt_01 + v2.2.defense_01 + v2.3.guardHerds_01 +v2.4.herding_01 + newcarry_01 +MeanAnnualTemperature+ FarmingPropensity +AnimalHusbandry+NumberParagraphs+
                 (1|gr(Language_ID, cov = pmat)) + (1|gr(Language_ID2, cov=spmat)), data = list(df1),
               prior = prior, family = "bernoulli", control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE),
               cores=4, data2 = list(pmat = pmat, spmat = spmat), iter=7000,thin=10) # default is warmup = iter/2
summary(fit_all)
r2_bayes (fit_all)


# NEG CARE data
df1<-dfprep(df, "negcare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat,Fcols01=Fcols01)$df1
pmat<-dfprep(df, "negcare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$phylo_covar_mat
spmat<-dfprep(df, "negcare", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$spatial_covar_mat
fit_all <- brm(formula=focalcol ~ v2.1.hunt_01 + v2.2.defense_01 + v2.3.guardHerds_01 +v2.4.herding_01 + newcarry_01 +MeanAnnualTemperature+ FarmingPropensity +AnimalHusbandry+NumberParagraphs+
                 (1|gr(Language_ID, cov = pmat)) + (1|gr(Language_ID2, cov=spmat)), data = list(df1),
               prior = prior, family = "bernoulli", control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE),
               cores=4, data2 = list(pmat = pmat, spmat = spmat), iter=7000,thin=10) # default is warmup = iter/2
summary(fit_all)
r2_bayes (fit_all)


# personhood data
df1<-dfprep(df, "personhood", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat,Fcols01=Fcols01)$df1
pmat<-dfprep(df, "personhood", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$phylo_covar_mat
spmat<-dfprep(df, "personhood", phylo_covar_mat = phylo_covar_mat, spatial_covar_mat = spatial_covar_mat, Fcols01 = Fcols01)$spatial_covar_mat
fit_all <- brm(formula=focalcol ~ v2.1.hunt_01 + v2.2.defense_01 + v2.3.guardHerds_01 +v2.4.herding_01 + newcarry_01 +MeanAnnualTemperature+ FarmingPropensity +AnimalHusbandry+NumberParagraphs+
                 (1|gr(Language_ID, cov = pmat)) + (1|gr(Language_ID2, cov=spmat)), data = list(df1),
               prior = prior, family = "bernoulli", control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE),
               cores=4, data2 = list(pmat = pmat, spmat = spmat), iter=7000,thin=10) # default is warmup = iter/2
summary(fit_all)
r2_bayes (fit_all)
