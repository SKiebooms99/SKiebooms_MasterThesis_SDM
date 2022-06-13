### ModelAdeno Europe ###

setwd("Z:/projects-unil/GEN4MIG/Solene")

library(biomod2)
library(raster)
library(foreign)
library(dismo)
library(ggplot2)
library(reshape2)
library(terra)
library(dplyr)

## télécharger les données (pour l'Europe) ##

# Espèce:
Adeno<-read.table("1_Occurences_data/Occ_Adenostyle_glabra_alpina/AdenostylesPV.csv",header = TRUE, sep = ",", row.names=1)
Adeno<-Adeno[,c("decimalLon","decimalLat")] #2010
# Supprimer les occurences se trouvant à une altitude aberrante (celles se trouvant en dehors du range altitudinal de l'sp):
altitude<-raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif") # doi:10.1038/sdata.2018.40.
AltPoints<-extract(altitude,Adeno)
#View(AltPoints)
Adeno <- Adeno[AltPoints >= 900,] #On retire les points en dessous de 900m -> reste 1694

# Charger les variables sélectionnées (ClimTopo_sub):
ClimTopo_sub<-stack("ClimTopo_sub_Current.tif")

# Générer des background points (Bp) pour l'Europe (et juste dans les biomes comprenant la majorité des points d'occurence)
biome <- raster("shapeEurope/GB_EBR_map_BiomesUE.tif")
crs(biome)
Occ_AdenoProj <- SpatialPointsDataFrame(coords = Adeno,data = Adeno,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
Occ_AdenoProj <- spTransform(Occ_AdenoProj,CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
Occ_AdenoProj <- as.data.frame(Occ_AdenoProj)[,-c(1:2)]
biomePoints <- extract(biome, Occ_AdenoProj)
BiomeToKeep <- as.numeric(names(table(biomePoints))[table(biomePoints)>= (5*nrow(Occ_AdenoProj)/100)]) # on excluera les biomes contenant moins de 5% des occurences
BiomeToKeep

biomeToSample <- rasterToPolygons(biome,fun=function(x){x==1 | x==7}) #sélectionner le(s) biome(s) à garder
Points <- extract(ClimTopo_sub, Adeno)
Bp <- spsample(biomeToSample, 20000, type="random")
Bp <- spTransform(Bp,CRS("+proj=longlat +datum=WGS84 +no_defs "))
Bp <- as.data.frame(Bp)
BpEnv <- raster::extract(ClimTopo_sub,Bp)
Bp <- na.omit(cbind(Bp,BpEnv)) #retirer les NA dans les lignes
#CheckPointInTheWater <- extract(subset(ClimTopo_sub,1),Bp)
#Bp <- Bp[!is.na(CheckPointInTheWater),]
Bp <- Bp[1:10000,1:2]
write.csv(Bp,"BpCoordinates_Adeno_EU2.csv")


#Fusionner les données de présences et les background points:
Bp<-read.csv("BpCoordinates_Adeno_EU2.csv", row.names = 1)
colnames(Adeno) = c("x","y")
Occ_Adeno <- gridSample(Adeno,ClimTopo_sub) # Ne garder qu'une occurence max par pixel -> 1226
Data_Adeno<-rbind(Occ_Adeno,Bp)
write.csv(Data_Adeno, "Data_Adeno_Occ_Bp_ClimTopo_EU2.csv")

###
#Graphique: densité des points (présence et Bp) pour chaque variable pour l'sp considérée
env <- raster::extract(ClimTopo_sub,Data_Adeno)
resp <- c(rep("presence", nrow(Occ_Adeno)), rep("background",10000))
tab <- cbind.data.frame(resp=resp,env)
tab$resp=as.factor(tab$resp)
tab <- melt(tab,id="resp")
rm(resp)
png("DensityPlotsVariableEnv_Adeno_EU.png",width= 5000,height = 8000,res=300)
p <- ggplot(data=tab,aes(x=value,color=resp,fill=resp)) + geom_density(alpha=0.1) +facet_wrap(~variable,scale= "free")
p
dev.off()
###

write.table(resp,"resp_Adeno.txt") # pour le cluster de calibration

## -> Direction cluster pour la calibration et la génération du modèle (Code "Calibration.R")
##################################################################################################

# Formater les données pour utiliser le package biomod2:
env <- raster::extract(ClimTopo_sub,Data_Adeno)
resp <- c(rep(1, nrow(Occ_Adeno)), rep(NA,10000))

Data_Adeno_EU<-BIOMOD_FormatingData(resp.var = resp,
                                    expl.var = env,
                                    resp.xy = Data_Adeno,
                                    resp.name = "AdenoModelEU",
                                    PA.nb.rep = 1,
                                    PA.nb.absences = 10000,
                                    PA.strategy = 'random')

Adeno_opt<-BIOMOD_ModelingOptions(GBM = list(n.trees = 5000),GAM = list(algo = "GAM_mgcv"))

Adeno_models_EU<-BIOMOD_Modeling(data=Data_Adeno_EU, models=c("GLM","GBM","GAM"),
                                 models.options = Adeno_opt,NbRunEval = 10,DataSplit = 70,
                                 models.eval.meth = c("TSS", "ROC"),
                                 VarImport = 10,SaveObj=TRUE,do.full.models = FALSE,
                                 Prevalence = 0.5)

# Vérifier importance des variables:
Adeno_models_var_import_ClimTopo <- get_variables_importance(Adeno_models_EU) 
apply(Adeno_models_var_import_ClimTopo, c(1, 2), mean) #The higher the score, the more important the variable
# Response curves
Adeno_glm <- BIOMOD_LoadModels(Adeno_models_EU, models = 'GLM')
Adeno_gbm <- BIOMOD_LoadModels(Adeno_models_EU, models = 'GBM')
Adeno_gam <- BIOMOD_LoadModels(Adeno_models_EU, models = 'GAM')

glm_eval_strip <- biomod2::response.plot2(models = Adeno_glm,
                                          Data = get_formal_data(Adeno_models_EU,
                                                                 'expl.var'),
                                          show.variables = get_formal_data(Adeno_models_EU,
                                                                           'expl.var.names'),
                                          do.bivariate = FALSE,
                                          fixed.var.metric = 'mean',
                                          legend = FALSE,
                                          display_title = FALSE,
                                          data_species = get_formal_data(Adeno_models_EU,'resp.var'))
gbm_eval_strip <-
  biomod2::response.plot2(models =Adeno_gbm,
                          Data = get_formal_data(Adeno_models_EU,
                                                 'expl.var'),
                          show.variables = get_formal_data(Adeno_models_EU,
                                                           'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(Adeno_models_EU,
                                                         'resp.var'))
gam_eval_strip <-
  biomod2::response.plot2(models = Adeno_gam,
                          Data = get_formal_data(Adeno_models_EU,
                                                 'expl.var'),
                          show.variables= get_formal_data(Adeno_models_EU,
                                                          'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(Adeno_models_EU,
                                                         'resp.var'))

# Scores des modèles:
Adeno_models_scores_EU <- get_evaluations(Adeno_models_EU)
Adeno_models_scores_EU
TSS_Adeno <- as.data.frame(Adeno_models_scores_EU["TSS",,,,])
AUC_Adeno <- as.data.frame(Adeno_models_scores_EU["ROC",,,,])
write.csv(TSS_Adeno,"Adeno_model_EU_TSS.csv")
write.csv(AUC_Adeno,"Adeno_model_EU_AUC.csv")

## Ensemble Modelling ##

Adeno_ensemble_models_EU <-BIOMOD_EnsembleModeling(modeling.output = Adeno_models_EU,
                                                   em.by = 'all',eval.metric = 'TSS',
                                                   eval.metric.quality.threshold = 0,
                                                   models.eval.meth = c("TSS", "ROC"),
                                                   prob.mean = FALSE,
                                                   prob.cv = TRUE,
                                                   committee.averaging = TRUE,
                                                   prob.mean.weight = TRUE,
                                                   VarImport = 0)
(Adeno_ensemble_models_scores_EU <- get_evaluations(Adeno_ensemble_models_EU)) #fiabilité 
TSS_Adeno_EM <- Adeno_ensemble_models_EU[[2]]
AUC_Adeno_EM <- Adeno_ensemble_models_EU[[3]]
write.csv(TSS_Adeno_EM,"Adeno_Emodel_EU_TSS.csv")
write.csv(AUC_Adeno_EM,"Adeno_Emodel_EU_AUC.csv")
##############################################################################################
### Projections ###

Adeno_models_EU<-get(load("Adeno/Adeno.1649146839.models.out"))
Adeno_ensemble_models_EU<-get(load("Adeno/Adeno.1649146839ensemble.models.out"))

## Projection actuelle ##
bio2<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/1981-2010/bio/CHELSA_bio2_1981-2010_V.2.1.tif")
bio6<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/1981-2010/bio/CHELSA_bio6_1981-2010_V.2.1.tif")
bio12<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif")
bio15<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/1981-2010/bio/CHELSA_bio15_1981-2010_V.2.1.tif")
slope<-raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

bio2<-crop(bio2, Mask_EU)
bio6<-crop(bio6, Mask_EU)
bio12<-crop(bio12, Mask_EU)
bio15<-crop(bio15, Mask_EU)
slope<-crop(slope, Mask_EU)

new_env_current<-stack(bio2,bio6,bio12,bio15,slope)

ForWater <- raster("Shapes_Masks/shapeEurope/ForNAWater.tif")
ForWater <- terra::crop(ForWater,Mask_EU)
new_env_current <- mask(new_env_current,ForWater)

names(new_env_current)<-c("ClimTopo_sub_Current.1","ClimTopo_sub_Current.2","ClimTopo_sub_Current.3","ClimTopo_sub_Current.4","ClimTopo_sub_Current.5")
Adeno_models_proj_current_EU <- BIOMOD_Projection(modeling.output = Adeno_models_EU,
                                                   new.env = stack(new_env_current),
                                                   proj.name = "Adeno_Current_EU_2",
                                                   binary.meth = "TSS",
                                                   selected.models="all")

#BIOMOD_EnsembleProjection
Adeno_ensemble_models_proj_current_EU <- BIOMOD_EnsembleForecasting(EM.output = Adeno_ensemble_models_EU,
                                                                     projection.output = Adeno_models_proj_current_EU,
                                                                     binary.meth = "TSS",
                                                                     selected.models="all")
plot(Adeno_ensemble_models_proj_current_EU)
# Save:
cartes_Adeno <- stack("Adeno/proj_Adeno_Current_EU_2/proj_Adeno_Current_EU_2_Adeno_ensemble.grd")
writeRaster(cartes_Adeno$Adeno_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"Adeno/Adeno_Current_EU_weighted_mean_2.tif")
carte_AdenoEU <- raster("Adeno/Adeno_Current_EU_weighted_mean_2.tif")
plot(carte_AdenoEU)

## Projections futures ##

Mask_EU<-rast("Shapes_Masks/shapeEurope/Mask_EU.tif")
Mask_EU<-terra::project(Mask_EU, "epsg:4326")
#couper le masque à 35° Est (réduire surface = réduire nombre de pixels = réduire temps de calculs):
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)

ForWater <- raster("Shapes_Masks/shapeEurope/ForNAWater.tif")
ForWater <- terra::crop(ForWater,Mask_EU)

# GCM = MPI-ESM

bio2_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio2_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio2_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio2_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio6_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio6_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio12_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio12_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio15_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio15_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
slope<- raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

#Couper les cartes sur la zone d'étude:
bio2_RCP26<-crop(bio2_RCP26, Mask_EU)
bio2_RCP85<-crop(bio2_RCP85, Mask_EU)
bio6_RCP26<-crop(bio6_RCP26, Mask_EU)
bio6_RCP85<-crop(bio6_RCP85, Mask_EU)
bio12_RCP26<-crop(bio12_RCP26, Mask_EU)
bio12_RCP85<-crop(bio12_RCP85, Mask_EU)
bio15_RCP26<-crop(bio15_RCP26, Mask_EU)
bio15_RCP85<-crop(bio15_RCP85, Mask_EU)
slope <-crop(slope, extent(Mask_EU))

# Stack des variables (+ reprendre la pente):
BioClim_RCP26<-stack(bio2_RCP26,bio6_RCP26,bio12_RCP26,bio15_RCP26,slope)
BioClim_RCP85<-stack(bio2_RCP85,bio6_RCP85,bio12_RCP85,bio15_RCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26 <- mask(BioClim_RCP26,ForWater)
BioClim_RCP85 <- mask(BioClim_RCP85,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
#names(BioClim_RCP26)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21","CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")
names(BioClim_RCP26)<-c("ClimTopo_sub_Current.1","ClimTopo_sub_Current.2","ClimTopo_sub_Current.3","ClimTopo_sub_Current.4","ClimTopo_sub_Current.5")

# Projection future optimiste
new_env<-stack(BioClim_RCP26)
Adeno_models_proj_2100_RCP26 <-BIOMOD_Projection(modeling.output = Adeno_models_EU,
                                                       new.env = new_env,
                                                       proj.name = "Adeno_2100_RCP26_T2",
                                                       binary.meth = "TSS",
                                                       selected.models = 'all')

Adeno_ensemble_models_proj_2100_RCP26 <- BIOMOD_EnsembleForecasting(EM.output = Adeno_ensemble_models_EU,
                                                                          projection.output = Adeno_models_proj_2100_RCP26,
                                                                          binary.meth = "TSS",
                                                                          selected.models = 'all')
plot(Adeno_ensemble_models_proj_2100_RCP26)
#Save:
cartes_Adeno_RCP26 <- stack("Adeno/proj_Adeno_2100_RCP26_T2/proj_Adeno_2100_RCP26_T2_Adeno_ensemble.grd")
writeRaster(cartes_Adeno_RCP26$Adeno_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"Adeno/Adeno_2100_RCP26_MPI_weighted_mean_T2.tif")
carte_AdenoEU_RCP26 <- raster("Adeno/Adeno_2100_RCP26_MPI_weighted_mean_T2.tif")
plot(carte_AdenoEU_RCP26)

# Projection future pessimiste
new_env_85<-stack(BioClim_RCP85)
names(new_env_85)<-c("ClimTopo_sub_Current.1","ClimTopo_sub_Current.2","ClimTopo_sub_Current.3","ClimTopo_sub_Current.4","ClimTopo_sub_Current.5")

Adeno_models_proj_2100_RCP85 <-BIOMOD_Projection(modeling.output = Adeno_models_EU,
                                                       new.env = new_env_85,
                                                       proj.name = "Adeno_2100_RCP85_MPI_T2",
                                                       binary.meth = "TSS",
                                                       selected.models = 'all')

Adeno_ensemble_models_proj_2100_RCP85 <- BIOMOD_EnsembleForecasting(EM.output = Adeno_ensemble_models_EU,
                                                                          projection.output = Adeno_models_proj_2100_RCP85,
                                                                          binary.meth = "TSS",
                                                                          selected.models = 'all')
plot(Adeno_ensemble_models_proj_2100_RCP85)
#Save:
cartes_Adeno_RCP85_MPI <- stack("Adeno/proj_Adeno_2100_RCP85_MPI_T2/proj_Adeno_2100_RCP85_MPI_T2_Adeno_ensemble.grd")
writeRaster(cartes_Adeno_RCP85_MPI$Adeno_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"Adeno/Adeno_2100_RCP85_MPI_weighted_mean_T2.tif")
carte_AdenoEU_RCP85_MPI <- raster("Adeno/Adeno_2100_RCP85_MPI_weighted_mean_T2.tif")
plot(carte_AdenoEU_RCP85_MPI)

# GCM = UKESM

bio2_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio2_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio2_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio2_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio6_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio6_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio6_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio6_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio12_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio12_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio12_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio12_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio15_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio15_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
slope<-raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

#Couper les cartes sur la zone d'étude:
bio2_UKRCP26<-crop(bio2_UKRCP26, Mask_EU)
bio2_UKRCP85<-crop(bio2_UKRCP85, Mask_EU)
bio6_UKRCP26<-crop(bio6_UKRCP26, Mask_EU)
bio6_UKRCP85<-crop(bio6_UKRCP85, Mask_EU)
bio12_UKRCP26<-crop(bio12_UKRCP26, Mask_EU)
bio12_UKRCP85<-crop(bio12_UKRCP85, Mask_EU)
bio15_UKRCP26<-crop(bio15_UKRCP26, Mask_EU)
bio15_UKRCP85<-crop(bio15_UKRCP85, Mask_EU)
slope<-crop(slope, Mask_EU)

# Stack des variables (+ reprendre la pente):
BioClim_RCP26_UKESM<-stack(bio2_UKRCP26,bio6_UKRCP26,bio12_UKRCP26,bio15_UKRCP26,slope)
BioClim_RCP85_UKESM<-stack(bio2_UKRCP85,bio6_UKRCP85,bio12_UKRCP85,bio15_UKRCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26_UKESM <- mask(BioClim_RCP26_UKESM,ForWater)
BioClim_RCP85_UKESM <- mask(BioClim_RCP85_UKESM,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_RCP26_UKESM)<-c("ClimTopo_sub_Current.1","ClimTopo_sub_Current.2","ClimTopo_sub_Current.3","ClimTopo_sub_Current.4","ClimTopo_sub_Current.5")

# Projection future optimiste
new_env_26_UK<-stack(BioClim_RCP26_UKESM)
Adeno_models_proj_2100_RCP26_UKESM <-BIOMOD_Projection(modeling.output = Adeno_models_EU,
                                                             new.env = new_env_26_UK,
                                                             proj.name = "Adeno_2100_RCP26_UKESM_T2",
                                                             binary.meth = "TSS",
                                                             selected.models = 'all')

Adeno_ensemble_models_proj_2100_RCP26_UKESM <- BIOMOD_EnsembleForecasting(EM.output = Adeno_ensemble_models_EU,
                                                                                projection.output = Adeno_models_proj_2100_RCP26_UKESM,
                                                                                binary.meth = "TSS",
                                                                                selected.models = 'all')
plot(Adeno_ensemble_models_proj_2100_RCP26_UKESM)
#Save:
cartes_Adeno_RCP26_UKESM <- stack("Adeno/proj_Adeno_2100_RCP26_UKESM_T2/proj_Adeno_2100_RCP26_UKESM_T2_Adeno_ensemble.grd")
writeRaster(cartes_Adeno_RCP26_UKESM$Adeno_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"Adeno/Adeno_2100_RCP26_UKESM_weighted_mean_T2.tif",overwrite=TRUE)
carte_AdenoEU_RCP26_UKESM <- raster("Adeno/Adeno_2100_RCP26_UKESM_weighted_mean_T2.tif")
plot(carte_AdenoEU_RCP26_UKESM)

# Projection future pessimiste
new_env_85_UKESM<-stack(BioClim_RCP85_UKESM)
names(new_env_85_UKESM)<-c("ClimTopo_sub_Current.1","ClimTopo_sub_Current.2","ClimTopo_sub_Current.3","ClimTopo_sub_Current.4","ClimTopo_sub_Current.5")
Adeno_models_proj_2100_RCP85_UKESM <-BIOMOD_Projection(modeling.output = Adeno_models_EU,
                                                             new.env = new_env_85_UKESM,
                                                             proj.name = "Adeno_2100_RCP85_UKESM_T2",
                                                             binary.meth = "TSS",
                                                             selected.models = 'all')

Adeno_ensemble_models_proj_2100_RCP85_UKESM <- BIOMOD_EnsembleForecasting(EM.output = Adeno_ensemble_models_EU,
                                                                                projection.output = Adeno_models_proj_2100_RCP85_UKESM,
                                                                                binary.meth = "TSS",
                                                                                selected.models = 'all')
plot(Adeno_ensemble_models_proj_2100_RCP85_UKESM)
#Save:
cartes_Adeno_RCP85_UKESM <- stack("Adeno/proj_Adeno_2100_RCP85_UKESM_T2/proj_Adeno_2100_RCP85_UKESM_T2_Adeno_ensemble.grd")
writeRaster(cartes_Adeno_RCP85_UKESM$Adeno_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"Adeno/Adeno_2100_RCP85_UKESM_weighted_mean_T2.tif")
carte_AdenoEU_RCP85_UKESM <- raster("Adeno/Adeno_2100_RCP85_UKESM_weighted_mean_T2.tif")
plot(carte_AdenoEU_RCP85_UKESM)


## Projection passée: LGM

# Télécharger variables bioclim pour le LGM (PMIP3, CCSM4):
bio2_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_02.tif")
bio6_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_06.tif")
bio12_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_12.tif")
bio15_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_15.tif")

# + recalculer la pente (fait dans SAGAGis):
slope_LGM<-raster("3_VariablesEnvi_horsBioclim/lgm/slope_lgm.tif")

# Couper les cartes sur la zone d'étude:
bio2_LGM<-crop(bio2_LGM, Mask_EU)
bio6_LGM<-crop(bio6_LGM, Mask_EU)
bio12_LGM<-crop(bio12_LGM, Mask_EU)
bio15_LGM<-crop(bio15_LGM, Mask_EU)
slope_LGM<-crop(slope_LGM, Mask_EU)
slope_LGM <- resample(slope_LGM,bio15_LGM)

# Convertir les unités (pour que ce soit les mêmes que pour les proj current et futures):
bio2_LGM<-(bio2_LGM/10)-273.15
bio6_LGM<-(bio6_LGM/10)-273.15
bio12_LGM<-bio12_LGM/10
bio15_LGM<-bio15_LGM/10

# Stack des variables (+ reprendre la pente):
BioClim_LGM<-stack(bio2_LGM,bio6_LGM,bio12_LGM,bio15_LGM,slope_LGM)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_LGM)<-c("ClimTopo_sub_Current.1","ClimTopo_sub_Current.2","ClimTopo_sub_Current.3","ClimTopo_sub_Current.4","ClimTopo_sub_Current.5")

# Projection passée (LGM: -21 000 BP):
new_env_LGM<-stack(BioClim_LGM)
Adeno_models_proj_LGM <-BIOMOD_Projection(modeling.output = Adeno_models_EU,
                                           new.env = new_env_LGM,
                                           proj.name = "Adeno_LGM_2",
                                           binary.meth = "TSS",
                                           selected.models = 'all')

Adeno_ensemble_models_proj_LGM <- BIOMOD_EnsembleForecasting(EM.output = Adeno_ensemble_models_EU,
                                                              projection.output = Adeno_models_proj_LGM,
                                                              binary.meth = "TSS",
                                                              selected.models = 'all')
plot(Adeno_ensemble_models_proj_LGM)
#Save:
cartes_Adeno_LGM <- stack("Adeno/proj_Adeno_LGM_2/proj_Adeno_LGM_2_Adeno_ensemble.grd")
writeRaster(cartes_Adeno_LGM$Adeno_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"Adeno/Adeno_LGM_weighted_mean_2.tif")
carte_Adeno_LGM <- raster("Adeno/Adeno_LGM_weighted_mean_2.tif")
plot(carte_Adeno_LGM)

##############################################################################################################################################
### Comparer les distributions ###

# Charger les projections (carteWM):
Proj_Current<-raster("Adeno/Adeno_Current_weightedmean_2.tif")
Proj_LGM_Adeno<-raster("Adeno/Adeno_LGM_weighted_mean_2.tif")
Proj_MPI26<-raster("Adeno/Adeno_2100_RCP26_MPI_weighted_mean_T2.tif")
Proj_MPI85<-raster("Adeno/Adeno_2100_RCP85_MPI_weighted_mean_T2.tif")
Proj_UK26<-raster("Adeno/Adeno_2100_RCP26_UKESM_weighted_mean_T2.tif")
Proj_UK85<-raster("Adeno/Adeno_2100_RCP85_UKESM_weighted_mean_T2.tif")

# Pour Proj_LGM: retirer les zones sous glace comme étant favorable:
Glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # charger carte présence de glaciers
# Couper carte glacier sur l'EU:
Mask_EU<-rast("Y_Shapes_Masks/shapeEurope/ForNaWater.tif")
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
Glacier<-crop(Glacier,Mask_EU)
glacier<-resample(Glacier,Proj_LGM_Adeno,method='ngb')
values(glacier)[is.na(values(glacier))] = 0 #sinon ça vire là où il y a des NA et donc quand il n'y a pas de glacier
glacierAbsence <- 1-glacier # pour mettre la favorabilité à 0 quand présence de glacier
Proj_LGM_HG_Adeno<-Proj_LGM_Adeno*glacierAbsence
plot(Proj_LGM_HG_Adeno) 
writeRaster(Proj_LGM_HG_Adeno,"6_Cartes_Projections/Adeno_LGM_WM_2.tif")

## Binariser les cartes de projection (en sélectionnant le treshold qui maximise le TSS, cf. Sp_Model_TSS):
# Treshold ici = 560
Adeno_bin_current<-BinaryTransformation(Proj_Current, 560)
Adeno_bin_2100_MPI26<-BinaryTransformation(Proj_MPI26, 560)
Adeno_bin_2100_MPI85<-BinaryTransformation(Proj_MPI85, 560)
Adeno_bin_2100_UK26<-BinaryTransformation(Proj_UK26, 560)
Adeno_bin_2100_UK85<-BinaryTransformation(Proj_UK85, 560)
Adeno_bin_LGM<-BinaryTransformation(Proj_LGM, 560)
plot(Adeno_bin_current)
plot(Adeno_bin_2100_MPI26)
plot(Adeno_bin_2100_MPI85)
plot(Adeno_bin_2100_UK26)
plot(Adeno_bin_2100_UK85)
plot(Adeno_bin_LGM)
writeRaster(Adeno_bin_current,"Adeno/Adeno_bin_current.tif")
writeRaster(Adeno_bin_2100_MPI26,"Adeno/Adeno_bin_2100_MPI26.tif")
writeRaster(Adeno_bin_2100_MPI85,"Adeno/Adeno_bin_2100_MPI85.tif")
writeRaster(Adeno_bin_2100_UK26,"Adeno/Adeno_bin_2100_UK26.tif")
writeRaster(Adeno_bin_2100_UK85,"Adeno/Adeno_bin_2100_UK85.tif")
writeRaster(Adeno_bin_LGM,"Adeno/Adeno_bin_LGM.tif")

## Rangeshifts:

# Current-LGM (par étapes pour comprendre)
#diff<-2*Pedicularis_bin_current + Pedicularis_bin_LGM 
#plot(diff)  # diff = 1: présence qu'au LGM, diff=2: présence qu'au présent, diff=3: présence aux 2 périodes
#valPixel<-freq(diff)
#valPixel # affiche nb de pixels de présence dans chaque catégorie (1, 2 ou 3)
#valPixel <- valPixel[-5,] # retirer les NA
# % de diff:
#PixelConserve <- valPixel[valPixel[,1]==3,"count"]/ (valPixel[valPixel[,1]==3,"count"] + valPixel[valPixel[,1]==2,"count"]) # % de l'aire de distribution conservé
#PixelPresenceOnly <- valPixel[valPixel[,1]==2,"count"]/ (valPixel[valPixel[,1]==3,"count"] + valPixel[valPixel[,1]==2,"count"]) # ou, en complémentaire, % de l'aire perdue
#DistPastParRapPres <- (valPixel[valPixel[,1]==1,"count"] +valPixel[valPixel[,1]==3,"count"])/ (valPixel[valPixel[,1]==3,"count"] + valPixel[valPixel[,1]==2,"count"]) # % de l'aire actuelle occupée au passé

# Comparer chaque proj au current (utilisation de la fonction BIOMOD_RangeSize, plus rapide):

# Current-LGM
Shift_LGM_current<-BIOMOD_RangeSize(Adeno_bin_current,Adeno_bin_LGM)
Shift_LGM_current<-Shift_LGM_current$Compt.By.Models
Shift_LGM_current<-as.data.frame(Shift_LGM_current)
#plot(Shift_LGM_current_$Diff.By.Pixel) # 0 = pixels non-occupé, 1 = pixels occupés slmt ds la nouvelle période
# -1 = pixels occupés aux 2 périodes, -2 = pixels perdus entre current et nouvelle période
diff<-2*Adeno_bin_current + Adeno_bin_LGM 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Adeno_Diff_LGM.tif")

# Current-MPI26
Shift_current_2100_MPI26 <- BIOMOD_RangeSize(Adeno_bin_current, Adeno_bin_2100_MPI26)
Shift_current_2100_MPI26<-Shift_current_2100_MPI26$Compt.By.Models
Shift_current_2100_MPI26<-as.data.frame(Shift_current_2100_MPI26)
# Pour le plot (plus lisible que plot diff.by.pixel):
diff<-2*Adeno_bin_current + Adeno_bin_2100_MPI26 
plot(diff)  # diff = 1: présence qu'au future (gain, en orange), diff=2: présence qu'au présent (loss, en vert clair), diff=3: présence aux 2 périodes (stable, en vert foncé), 0 = pixels non-occupés
writeRaster(diff, "6_Cartes_Projections/Adeno_Diff_MPI26.tif")

# Current-MPI85
Shift_current_2100_MPI85 <- BIOMOD_RangeSize(Adeno_bin_current, Adeno_bin_2100_MPI85)
Shift_current_2100_MPI85<-Shift_current_2100_MPI85$Compt.By.Models
Shift_current_2100_MPI85<-as.data.frame(Shift_current_2100_MPI85)
diff<-2*Adeno_bin_current + Adeno_bin_2100_MPI85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Adeno_Diff_MPI85.tif")

# Current-UK26
Shift_current_2100_UK26 <- BIOMOD_RangeSize(Adeno_bin_current, Adeno_bin_2100_UK26)
Shift_current_2100_UK26<-Shift_current_2100_UK26$Compt.By.Models
Shift_current_2100_UK26<-as.data.frame(Shift_current_2100_UK26)
diff<-2*Adeno_bin_current + Adeno_bin_2100_UK26 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Adeno_Diff_UK26.tif")

# Current-UK85
Shift_current_2100_UK85 <- BIOMOD_RangeSize(Adeno_bin_current, Adeno_bin_2100_UK85)
Shift_current_2100_UK85<-Shift_current_2100_UK85$Compt.By.Models
Shift_current_2100_UK85<-as.data.frame(Shift_current_2100_UK85)
diff<-2*Adeno_bin_current + Adeno_bin_2100_UK85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Adeno_Diff_UK85.tif")

# Fusionner les dataframe:
RangeShifts<-bind_rows(Shift_LGM_current,Shift_current_2100_MPI26,Shift_current_2100_MPI85,
                       Shift_current_2100_UK26,Shift_current_2100_UK85,.id=NULL)
write.csv(RangeShifts,"OutputsCluster/Adeno_RangeShifts.csv")


## Rangeshift altitudinal

# load cartes binarisées:
Adeno_bin_current<-raster("Adeno/Adeno_bin_current.tif")
Adeno_bin_2100_MPI26<-raster("Adeno/Adeno_bin_2100_MPI26.tif")
Adeno_bin_2100_MPI85<-raster("Adeno/Adeno_bin_2100_MPI85.tif")
Adeno_bin_2100_UK26<-raster("Adeno/Adeno_bin_2100_UK26.tif")
Adeno_bin_2100_UK85<-raster("Adeno/Adeno_bin_2100_UK85.tif")
Adeno_bin_LGM<-raster("Adeno/Adeno_bin_LGM.tif")

## Shift altitudinal:
altPres <- raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif")
altPast<- raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_dem_-200_V1.0.tif")

# Extraire coordonnées au niveau des présences pour chaque période:
# Current
coordPres <- coordinates(Adeno_bin_current)
valPres <- values(Adeno_bin_current)
coordPres <- coordPres[!is.na(valPres) & valPres == 1,]
altiPres <- extract(altPres,coordPres)
MedAltPres<-median(altiPres) # 1415,375 m
# LGM
coordLGM <- coordinates(Adeno_bin_LGM)
valLGM <- values(Adeno_bin_LGM)
coordLGM <- coordLGM[!is.na(valLGM) & valLGM == 1,]
altiLGM <- extract(altPast,coordLGM)
MedAltLGM<-median(altiLGM) #NULL
# MPI26
coordMPI26 <- coordinates(Adeno_bin_2100_MPI26)
valMPI26 <- values(Adeno_bin_2100_MPI26)
coordMPI26 <- coordMPI26[!is.na(valMPI26) & valMPI26 == 1,]
altiMPI26 <- extract(altPres,coordMPI26)
MedAltMPI26<-median(altiMPI26) # 1535,06 m
# MPI85
coordMPI85 <- coordinates(Adeno_bin_2100_MPI85)
valMPI85 <- values(Adeno_bin_2100_MPI85)
coordMPI85 <- coordMPI85[!is.na(valMPI85) & valMPI85 == 1,]
altiMPI85 <- extract(altPres,coordMPI85)
MedAltMPI85<-median(altiMPI85) # 1030,75 m
# UK26
coordUK26 <- coordinates(Adeno_bin_2100_UK26)
valUK26 <- values(Adeno_bin_2100_UK26)
coordUK26 <- coordUK26[!is.na(valUK26) & valUK26 == 1,]
altiUK26 <- extract(altPres,coordUK26)
MedAltUK26<-median(altiUK26) # 1831,625 m
# UK85
coordUK85 <- coordinates(Adeno_bin_2100_UK85)
valUK85<- values(Adeno_bin_2100_UK85)
coordUK85 <- coordUK85[!is.na(valUK85) & valUK85 == 1,]
altiUK85 <- extract(altPres,coordUK85)
MedAltUK85<-median(altiUK85) # 2273,625 m

# Mettre les médianes en graphique:
pourBoxplotPres <- cbind.data.frame(Altitude = altiPres, Periode = "Present")
#pourBoxplotLGM <- cbind.data.frame(Altitude = altiLGM, Periode = "LGM")
pourBoxplotMPI26 <- cbind.data.frame(Altitude = altiMPI26, Periode = "MPIESM-RCP26")
pourBoxplotMPI85 <- cbind.data.frame(Altitude = altiMPI85, Periode = "MPIESM-RCP85")
pourBoxplotUK26 <- cbind.data.frame(Altitude = altiUK26, Periode = "UKESM-RCP26")
pourBoxplotUK85 <- cbind.data.frame(Altitude = altiUK85, Periode = "UKESM-RCP85")

pourBoxplot <- rbind(pourBoxplotPres,pourBoxplotMPI26,pourBoxplotMPI85,
                     pourBoxplotUK26,pourBoxplotUK85)
# GGPlot2
pourBoxplot$Periode = as.factor(pourBoxplot$Periode)
p <- ggplot(data = pourBoxplot,aes(x=Periode,y=Altitude)) + geom_boxplot() + labs(title="Adenostyles glabra subsp alpina") +theme_classic() 
p


##############################################################################
### Calculer vitesses de migration

# A. Altitudinal

# RangeSchift alt (différences entre médianes; si diff<0: sp est descendue en alt par rapport au présent, si diff>0: sp est montée en alt par rapport au présent):
DiffAlt_Pres_LGM<-MedAltLGM-MedAltPres
DiffAlt_Pres_LGM # NULL
DiffAlt_Pres_MPI26<-MedAltMPI26-MedAltPres
DiffAlt_Pres_MPI26 # 119,69 m
DiffAlt_Pres_MPI85<-MedAltMPI85-MedAltPres
DiffAlt_Pres_MPI85 # -384,625 m
DiffAlt_Pres_UK26<-MedAltUK26-MedAltPres
DiffAlt_Pres_UK26 # 416,25 m
DiffAlt_Pres_UK85<-MedAltUK85-MedAltPres
DiffAlt_Pres_UK85 # 858,25 m

# Diviser les rangeshift altitudinal par deltaT:
# Entre LGM et Current (1995-(-21000)=22995 ans) # vitesse min de migration de l'sp
#MigrationAlt_Pres_LGM<-DiffAlt_Pres_LGM/22995
#MigrationAlt_Pres_LGM #  m/an

# Entre Current et future (2085-1995= 90 ans) selon scénario # vitesse de migration que l'sp doit avoir pour suivre son habitat optimal:
MigAlt_Pres_MPI26<-DiffAlt_Pres_MPI26/90
MigAlt_Pres_MPI26 # 1,33  m/an
MigAlt_Pres_MPI85<-DiffAlt_Pres_MPI85/90
MigAlt_Pres_MPI85 # -4,27  m/an
MigAlt_Pres_UK26<-DiffAlt_Pres_UK26/90
MigAlt_Pres_UK26 # 4,625 m/an
MigAlt_Pres_UK85<-DiffAlt_Pres_UK85/90
MigAlt_Pres_UK85 # 9,54 m/an

# B. latitudinal/longitudinal

# Chargement des coordonnées des centroides (calculé dans QGIS):
Centroides<-read.csv("ShiftLat_Adeno_2.csv",header=F,row.names=1, sep=";",dec=",")
# Calcul des distances entre les centroides des différents scenarios:
Distances<-pointDistance(Centroides, lonlat=T,allpairs=T)
