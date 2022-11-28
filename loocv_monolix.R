# Load my usual libraries that I use in every R script
# If running in Rstudio, Click session -> set working directory -> to source file location

library(ggplot2)
library(dplyr)
library(plyr)
library(readxl)
library(ggplot2)
library(lme4)
library(lmerTest)
library(knitr)
library(kableExtra)
library(tidyverse)
library(scales)
library(gtools)
library(utils)

# Load lixoft connectors (This library has to be manually installed and is included with Monolix)
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix")

# Initialize empty LOOCV data frame
LOOCV <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(LOOCV) <- c('Model', 'Virus', 'sq_error', 'iD')

# Load in SHIV data
origDataSHIV <-
  read.csv(file = "data/VOP PVL_Rebound_cw_lod.csv", header =
             TRUE, sep = ",")
getActiveDocumentContext()
for (i in 1:length(origDataSHIV$animal_id)) {
  # Load path that has all of my monolix files
  demoPath = getwd()
  # Load the specific monolix file that I've used
  project <- paste0(demoPath, "Damped Gomperts - SHIV.mlxtran")
  # Load the project and run population parameter estimation, so R can actually interface with the project
  loadProject(projectFile = project)
  
  runPopulationParameterEstimation()
  
  # remove the ith data point and sepearate the data frames into the removed data and the rest of the data
  crossvalidData <- origDataSHIV[-c(i),]
  removedData <- origDataSHIV[i,]
  #Write the new data frame with the removed data point to a new csv
  write.csv(crossvalidData, file = "data/VOP PVL_Rebound_cw_lod_crossvalid.csv", row.names = FALSE)
  #get the path of the new data frame with the removed data point
  filtereddatfile <-
    "data/VOP PVL_Rebound_cw_lod_crossvalid.csv"
  #Load the original data set from the monolix file
  BaseData <- getData()
  #Replace the original data set with the dataset with the removed data point,
  # and use the same headers and observation types as the original file
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  #Run population and individual (conditional) parameter estimation with the new 
  # data set
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  #Get the new parameter set from the monolix results
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$animal_id,]
  
  #Run the model with the new parameter set and evaluate the model at the time
  # from the removed data set
  V0 <- newParams$V0
  lambda0 <- newParams$lambda0
  beta <- newParams$beta
  gamma <- newParams$gamma
  tau <- removedData$week_post_treatment
  
  vl <- V0 * (exp((
    lambda0 * (-1 + exp(-(tau) * beta)) * beta + gamma *
      (1 - exp(-(tau) * beta) - (tau) * beta)
  ) / beta ^ 2))
  
  #calculate the squared error, and add that error to the LOOCV data frame
  sq_error <- (log10(removedData$viral_load) - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model),] <-
    c("DG", 'SHIV', as.numeric(sq_error), removedData$animal_id)
}
#Rinse and repeat this step for all models and all data sets
for (i in 1:length(origDataSHIV$animal_id)) {
  demoPath = getwd()
  project <- paste0(demoPath, "Two Stage Decay - SHIV.mlxtran")
  loadProject(projectFile = project)
  runPopulationParameterEstimation()
  
  crossvalidData <- origDataSHIV[-c(i), ]
  removedData <- origDataSHIV[i, ]
  write.csv(crossvalidData, file = "data/VOP PVL_Rebound_cw_lod_crossvalid.csv", row.names = FALSE)
  filtereddatfile <-
    "data/VOP PVL_Rebound_cw_lod_crossvalid.csv"
  BaseData <- getData()
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$animal_id, ]
  
  V0 <- newParams$V0
  lambda1 <- newParams$lambda1
  lambda2 <- newParams$lambda2
  cp <- newParams$cp
  tau <- removedData$week_post_treatment
  
  vl <-
    1 * (tau <= cp) * V0 * exp(-lambda1 * (tau)) + 1 * (tau > cp) * V0 * exp(-lambda1 * (cp)) * exp(-lambda2 * (tau - cp))
  
  sq_error <- (log10(removedData$viral_load) - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model), ] <-
    c("TSD", 'SHIV', as.numeric(sq_error), removedData$animal_id)
}

for (i in 1:length(origDataSHIV$animal_id)) {
  demoPath = getwd()
  project <- paste0(demoPath, "Two Stage Decay New - SHIV.mlxtran")
  loadProject(projectFile = project)
  runPopulationParameterEstimation()
  
  crossvalidData <- origDataSHIV[-c(i), ]
  removedData <- origDataSHIV[i, ]
  write.csv(crossvalidData, file = "data/VOP PVL_Rebound_cw_lod_crossvalid.csv", row.names = FALSE)
  filtereddatfile <-
    "data/VOP PVL_Rebound_cw_lod_crossvalid.csv"
  BaseData <- getData()
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$animal_id, ]
  
  V0 <- newParams$V0
  NKT <- newParams$NKT
  c <- newParams$c
  delta <- newParams$delta
  mu_m <- newParams$mu_m
  tau <- removedData$week_post_treatment
  
  vl <-
    V0 * ((1 - (NKT / (c - delta)) - ((c - NKT) / (c - mu_m))) * exp(-c * tau) + (NKT /
                                                                                    (c - delta)) * exp(-delta * tau) + ((c - NKT) / (c - mu_m)) * exp(-mu_m *
                                                                                                                                                        tau))
  
  sq_error <- (log10(removedData$viral_load) - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model), ] <-
    c("TSD", 'SHIV', as.numeric(sq_error), removedData$animal_id)
}

## SIV Portion

origDataSIV <-
  read.csv(file = "data/InfantR01 - SIV.csv", header =
             TRUE, sep = ",")

for (i in 1:length(origDataSIV$animal_id)) {
  demoPath = getwd()
  project <- paste0(demoPath, "Damped Gomperts - SIV.mlxtran")
  loadProject(projectFile = project)
  runPopulationParameterEstimation()
  
  crossvalidData <- origDataSIV[-c(i), ]
  removedData <- origDataSIV[i, ]
  write.csv(crossvalidData, file = "data/InfantR01 - SIV_crossvalid.csv", row.names = FALSE)
  filtereddatfile <-
    "data/InfantR01 - SIV_crossvalid.csv"
  BaseData <- getData()
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$animal_id, ]
  
  V0 <- newParams$V0
  lambda0 <- newParams$lambda0
  beta <- newParams$beta
  gamma <- newParams$gamma
  tau <- removedData$weeks_post_treatment
  
  vl <- V0 * (exp((
    lambda0 * (-1 + exp(-(tau) * beta)) * beta + gamma *
      (1 - exp(-(tau) * beta) - (tau) * beta)
  ) / beta ^ 2))
  
  sq_error <- (log10(removedData$viral_load) - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model), ] <-
    c("DG", 'SIV', as.numeric(sq_error), removedData$animal_id)
}

for (i in 1:length(origDataSIV$animal_id)) {
  demoPath = getwd()
  project <- paste0(demoPath, "Two Stage Decay New - SIV.mlxtran")
  loadProject(projectFile = project)
  runPopulationParameterEstimation()
  
  crossvalidData <- origDataSIV[-c(i), ]
  removedData <- origDataSIV[i, ]
  write.csv(crossvalidData, file = "data/InfantR01 - SIV_crossvalid.csv", row.names = FALSE)
  filtereddatfile <-
    "data/InfantR01 - SIV_crossvalid.csv"
  BaseData <- getData()
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$animal_id, ]
  
  V0 <- newParams$V0
  NKT <- newParams$NKT
  c <- newParams$c
  delta <- newParams$delta
  mu_m <- newParams$mu_m
  tau <- removedData$weeks_post_treatment
  
  vl <-
    V0 * ((1 - (NKT / (c - delta)) - ((c - NKT) / (c - mu_m))) * exp(-c * tau) + (NKT /
                                                                                    (c - delta)) * exp(-delta * tau) + ((c - NKT) / (c - mu_m)) * exp(-mu_m *
                                                                                                                                                        tau))
  
  sq_error <- (log10(removedData$viral_load) - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model), ] <-
    c("TSD", 'SIV', as.numeric(sq_error), removedData$animal_id)
}

## HIV Portion

origDataHIV <-
  read.csv(file = "data/humanData.csv", header =
             TRUE, sep = ",")

for (i in 1:length(origDataHIV$human_id)) {
  demoPath = getwd()
  project <- paste0(demoPath, "Damped Gomperts - HIV.mlxtran")
  loadProject(projectFile = project)
  runPopulationParameterEstimation()
  
  crossvalidData <- origDataHIV[-c(i), ]
  removedData <- origDataHIV[i, ]
  write.csv(crossvalidData, file = "data/humanData_crossvalid.csv", row.names = FALSE)
  filtereddatfile <-
    "data/humanData_crossvalid.csv"
  BaseData <- getData()
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$human_id, ]
  
  V0 <- newParams$V0
  lambda0 <- newParams$lambda0
  beta <- newParams$beta
  gamma <- newParams$gamma
  tau <- removedData$time_relative
  
  vl <- V0 * (exp((
    lambda0 * (-1 + exp(-(tau) * beta)) * beta + gamma *
      (1 - exp(-(tau) * beta) - (tau) * beta)
  ) / beta ^ 2))
  
  sq_error <- (removedData$log_vLoad - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model), ] <-
    c("DG", 'HIV', as.numeric(sq_error), removedData$human_id)
}

for (i in 1:length(origDataHIV$human_id)) {
  demoPath = getwd()
  project <- paste0(demoPath, "Two Stage Decay New - HIV.mlxtran")
  loadProject(projectFile = project)
  runPopulationParameterEstimation()
  
  crossvalidData <- origDataHIV[-c(i), ]
  removedData <- origDataHIV[i, ]
  write.csv(crossvalidData, file = "data/humanData_crossvalid.csv", row.names = FALSE)
  filtereddatfile <-
    "data/humanData_crossvalid.csv"
  BaseData <- getData()
  setData(filtereddatfile,
          BaseData$headerTypes,
          BaseData$observationTypes)
  
  runPopulationParameterEstimation()
  runConditionalModeEstimation()
  
  crossvalidParams <-
    getEstimatedIndividualParameters(method = "conditionalMode")
  newParams = crossvalidParams$conditionalMode[crossvalidParams$conditionalMode$id == removedData$human_id, ]
  
  V0 <- newParams$V0
  NKT <- newParams$NKT
  c <- newParams$c
  delta <- newParams$delta
  mu_m <- newParams$mu_m
  tau <- removedData$time_relative
  
  vl <-
    V0 * ((1 - (NKT / (c - delta)) - ((c - NKT) / (c - mu_m))) * exp(-c * tau) + (NKT /
                                                                                    (c - delta)) * exp(-delta * tau) + ((c - NKT) / (c - mu_m)) * exp(-mu_m *
                                                                                                                                                        tau))
  
  sq_error <- (removedData$log_vLoad  - log10(vl)) ^ 2
  LOOCV[1 + length(LOOCV$Model), ] <-
    c("TSD", 'HIV', as.numeric(sq_error), removedData$human_id)
}

#Finally, write the LOOCV data frame to a csv, including the standard error
# for each data/model combination for use in plots.
write.csv(LOOCV, file = "crossvalidData.csv", row.names = FALSE)
cdata <- ddply(
  LOOCV,
  c("Model", "Virus"),
  summarise,
  N    = length(sq_error),
  loocv = mean(as.numeric(sq_error)),
  sd   = sd(sq_error),
  se   = sd / sqrt(N)
)