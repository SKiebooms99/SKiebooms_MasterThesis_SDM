### Calibration des modèles ###

setwd("Z:/projects-unil/GEN4MIG/Solene")

##Libraries
library(biomod2)
library(raster)
library(foreign)
library(dismo)

NomSp<-as.character("Dryopteris")  # Insérer l'espèce à analyser

DataSp <- read.csv(paste0("Data_",NomSp,"_Occ_Bp_ClimTopo_EU2.csv"),row.names=1) #Données d'occurence spécifique + background points
Climtopo_sub <- stack("ClimTopo_sub_Current.tif") # Variables environnementales sélectionnees
resp <- read.table(paste0("resp_",NomSp,".txt"))[,1] #doit etre un vecteur numerique à la fin avec des 1 et NA
resp[resp=="presence"] = "1" 
resp[resp=="absence"] = NA 
resp <- as.numeric(resp)

env <- raster::extract(Climtopo_sub,DataSp)


DataSpEU<-BIOMOD_FormatingData(resp.var = resp,
                               expl.var = env,
                               resp.xy = DataSp,
                               resp.name = NomSp,
                               PA.nb.rep = 1,
                               PA.nb.absences = 10000,
                               PA.strategy = 'random')

ModOpt<-BIOMOD_ModelingOptions(GBM = list(n.trees = 5000),GAM = list(algo = "GAM_mgcv"))

SpModelsEU<-BIOMOD_Modeling(data=DataSpEU, models=c("GLM","GBM","GAM"),
                            models.options = ModOpt,NbRunEval = 10,DataSplit = 70,
                            models.eval.meth = c("TSS", "ROC"),
                            VarImport = 10,SaveObj=TRUE,do.full.models = FALSE,
                            Prevalence = 0.5)

# Vérifier importance des variables:
Sp_models_var_import_ClimTopo <- get_variables_importance(SpModelsEU) 
VarImport <- apply(Sp_models_var_import_ClimTopo, c(1, 2), mean) #The higher the score, the more important the variable
write.csv(VarImport,paste0("OutputsCluster/",NomSp,"_VariableImportance.csv"))

# Response curves
Sp_glm <- BIOMOD_LoadModels(SpModelsEU, models = 'GLM')
Sp_gbm <- BIOMOD_LoadModels(SpModelsEU, models = 'GBM')
Sp_gam <- BIOMOD_LoadModels(SpModelsEU, models = 'GAM')

glm_eval_strip <- biomod2::response.plot2(models = Sp_glm,
                                          Data = get_formal_data(SpModelsEU,
                                                                 'expl.var'),
                                          show.variables = get_formal_data(SpModelsEU,
                                                                           'expl.var.names'),
                                          do.bivariate = FALSE,
                                          fixed.var.metric = 'mean',
                                          legend = FALSE,
                                          display_title = FALSE,
                                          data_species = get_formal_data(SpModelsEU,'resp.var'))
gbm_eval_strip <-
  biomod2::response.plot2(models =Sp_gbm,
                          Data = get_formal_data(SpModelsEU,
                                                 'expl.var'),
                          show.variables = get_formal_data(SpModelsEU,
                                                           'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(SpModelsEU,
                                                         'resp.var'))
gam_eval_strip <-
  biomod2::response.plot2(models = Sp_gam,
                          Data = get_formal_data(SpModelsEU,
                                                 'expl.var'),
                          show.variables= get_formal_data(SpModelsEU,
                                                          'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(SpModelsEU,
                                                         'resp.var'))


# Scores des modèles:
Sp_models_scores_EU <- get_evaluations(SpModelsEU)
Sp_models_scores_EU
TSS_Sp <- as.data.frame(Sp_models_scores_EU["TSS",,,,])
AUC_Sp <- as.data.frame(Sp_models_scores_EU["ROC",,,,])
write.csv(TSS_Sp,paste0("OutputsCluster/",NomSp,"_model_EU_TSS.csv"))
write.csv(AUC_Sp,paste0("OutputsCluster/",NomSp,"_model_EU_AUC.csv"))

## Ensemble Modelling ##

Sp_ensemble_models_EU <-BIOMOD_EnsembleModeling(modeling.output = SpModelsEU,
                                                em.by = 'all',eval.metric = 'TSS',
                                                eval.metric.quality.threshold = 0,
                                                models.eval.meth = c("TSS", "ROC"),
                                                prob.mean = FALSE,
                                                prob.cv = F,
                                                committee.averaging = F,
                                                prob.mean.weight = TRUE,
                                                VarImport = 0)
Sp_ensemble_models_scores_EU <- get_evaluations(Sp_ensemble_models_EU) #fiabilité 
WMeanEval <- Sp_ensemble_models_scores_EU[[1]]
write.csv(WMeanEval,paste0("OutputsCluster/",NomSp,"_EvaluationEnsemble_WMean.csv"))
removeTmpFiles(h=0.1)
