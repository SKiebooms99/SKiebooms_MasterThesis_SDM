### Modèle Silène Europe ###

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
Silene<-read.table("1_Occurences_data/Occ_Silene_acaulis/SilenePV.csv",header = TRUE, sep = ",", row.names=1)
Silene<-Silene[,c("decimalLon","decimalLat")] #8865
# Supprimer les occurences se trouvant à une altitude aberrante (celles se trouvant en dehors du range altitudinal de l'sp):
altitude<-raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif") # doi:10.1038/sdata.2018.40.
AltPoints<-extract(altitude,Silene)
#View(AltPoints)
Silene <- Silene[AltPoints >= 900,] #On retire les points en dessous de 900m -> reste 4956

# Charger les variables sélectionnées (pour tout l'hémisphère nord)

# variables bioclimatiques (actuelles, à l'échelle de l'Europe, résolution 1km):
bioclim_map <- stack(list.files(path="CHELSA_bioclim/CHELSA_bioclim_1981_2010",pattern = ".tif",full.names = T))
slope<- raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

# Uniformiser la taille des rasters:
Mask<-rast("maskHemisphereNord.tif")
Mask<-terra::project(Mask, "epsg:4326")
#couper le masque à 30° N (réduire surface = réduire nombre de pixels = réduire temps de calculs):
Cut<-as(extent(-180,180,30,90),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask)        
Mask <- crop(Mask, Cut)
Mask <- raster(Mask)

bioclim_map<-crop(bioclim_map, Mask)
slope <-crop(slope, extent(Mask))

ClimTopo<-stack(bioclim_map,slope)

# Sélectionner potentiellememt d'autres prédicteurs? Regarder courbes de densité pour toutes les bio
# -> bio6, bio10, bio12, bio15, slope

ClimTopo_sub<-stack(subset(ClimTopo,c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                                      "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")))

# Générer des background points (Bp) pour l'Europe (et juste dans les biomes comprenant la majorité des points d'occurence)
biome <- raster("shapeEurope/GB_EBR_map_BiomesUE.tif")
crs(biome)
Occ_SileneProj <- SpatialPointsDataFrame(coords = Silene,data = Silene,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
Occ_SileneProj <- spTransform(Occ_SileneProj,CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
Occ_SileneProj <- as.data.frame(Occ_SileneProj)[,-c(1:2)]
biomePoints <- extract(biome, Occ_SileneProj)
BiomeToKeep <- as.numeric(names(table(biomePoints))[table(biomePoints)>= (5*nrow(Occ_SileneProj)/100)]) # on excluera les biomes contenant moins de 5% des occurences
BiomeToKeep

biomeToSample <- rasterToPolygons(biome,fun=function(x){x==1}) #sélectionner le(s) biome(s) à garder
Points <- extract(ClimTopo_sub, Silene)
Bp <- spsample(biomeToSample, 20000, type="random")
Bp <- spTransform(Bp,CRS("+proj=longlat +datum=WGS84 +no_defs "))
Bp <- as.data.frame(Bp)
BpEnv <- raster::extract(ClimTopo_sub,Bp)
Bp <- na.omit(cbind(Bp,BpEnv)) #retirer les NA dans les lignes
#CheckPointInTheWater <- extract(subset(ClimTopo_sub,1),Bp)
#Bp <- Bp[!is.na(CheckPointInTheWater),]
Bp <- Bp[1:10000,1:2]
write.csv(Bp,"BpCoordinates_Silene_EU2.csv")

#Fusionner les données de présences et les background points:
Bp<-read.csv("BpCoordinates_Silene_EU2.csv", row.names = 1)
colnames(Silene) = c("x","y")
Occ_Silene <- gridSample(Silene,ClimTopo_sub) # Ne garder qu'une occurence max par pixel -> 4956
Data_Silene<-rbind(Occ_Silene,Bp)
write.csv(Data_Silene, "Data_Silene_Occ_Bp_ClimTopo_EU2.csv")

###
#Graphique: densité des points (présence et Bp) pour chaque variable pour l'sp considérée
env <- raster::extract(ClimTopo_sub,Data_Silene)
resp <- c(rep("presence", nrow(Occ_Silene)), rep("background",10000))
tab <- cbind.data.frame(resp=resp,env)
tab$resp=as.factor(tab$resp)
tab <- melt(tab,id="resp")
rm(resp)
png("DensityPlotsVariableEnv_Silene_EU_bis.png",width= 5000,height = 8000,res=300)
p <- ggplot(data=tab,aes(x=value,color=resp,fill=resp)) + geom_density(alpha=0.1) +facet_wrap(~variable,scale= "free")
p
dev.off()
###

# Formater les données pour utiliser le package biomod2:
env <- raster::extract(ClimTopo_sub,Data_Silene)
resp <- c(rep(1, nrow(Occ_Silene)), rep(NA,10000))

Data_Silene_EU<-BIOMOD_FormatingData(resp.var = resp,
                                     expl.var = env,
                                     resp.xy = Data_Silene,
                                     resp.name = "SileneModelEU",
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 10000,
                                     PA.strategy = 'random')

Silene_opt<-BIOMOD_ModelingOptions(GBM = list(n.trees = 5000),GAM = list(algo = "GAM_mgcv"))

Silene_models_EU<-BIOMOD_Modeling(data=Data_Silene_EU, models=c("GLM","GBM","GAM"),
                                  models.options = Silene_opt,NbRunEval = 10,DataSplit = 70,
                                  models.eval.meth = c("TSS", "ROC"),
                                  VarImport = 10,SaveObj=TRUE,do.full.models = FALSE,
                                  Prevalence = 0.5)

# Vérifier importance des variables:
Silene_models_var_import_ClimTopo <- get_variables_importance(Silene_models_EU) 
apply(Silene_models_var_import_ClimTopo, c(1, 2), mean) #The higher the score, the more important the variable
# Response curves
Silene_glm <- BIOMOD_LoadModels(Silene_models_EU, models = 'GLM')
Silene_gbm <- BIOMOD_LoadModels(Silene_models_EU, models = 'GBM')
Silene_gam <- BIOMOD_LoadModels(Silene_models_EU, models = 'GAM')

glm_eval_strip <- biomod2::response.plot2(models = Silene_glm,
                                          Data = get_formal_data(Silene_models_EU,
                                                                 'expl.var'),
                                          show.variables = get_formal_data(Silene_models_EU,
                                                                           'expl.var.names'),
                                          do.bivariate = FALSE,
                                          fixed.var.metric = 'mean',
                                          legend = FALSE,
                                          display_title = FALSE,
                                          data_species = get_formal_data(Silene_models_EU,'resp.var'))
gbm_eval_strip <-
  biomod2::response.plot2(models =Silene_gbm,
                          Data = get_formal_data(Silene_models_EU,
                                                 'expl.var'),
                          show.variables = get_formal_data(Silene_models_EU,
                                                           'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(Silene_models_EU,
                                                         'resp.var'))
gam_eval_strip <-
  biomod2::response.plot2(models = Silene_gam,
                          Data = get_formal_data(Silene_models_EU,
                                                 'expl.var'),
                          show.variables= get_formal_data(Silene_models_EU,
                                                          'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(Silene_models_EU,
                                                         'resp.var'))

# Scores des modèles:
Silene_models_scores_EU <- get_evaluations(Silene_models_EU)
Silene_models_scores_EU
TSS_Silene <- as.data.frame(Silene_models_scores_EU["TSS",,,,])
AUC_Silene <- as.data.frame(Silene_models_scores_EU["ROC",,,,])
write.csv(TSS_Silene,"Silene_model_EU_TSS.csv")
write.csv(AUC_Silene,"Silene_model_EU_AUC.csv")

## Ensemble Modelling ##

Silene_ensemble_models_EU <-BIOMOD_EnsembleModeling(modeling.output = Silene_models_EU,
                                                    em.by = 'all',eval.metric = 'TSS',
                                                    eval.metric.quality.threshold = 0,
                                                    models.eval.meth = c("TSS", "ROC"),
                                                    prob.mean = FALSE,
                                                    prob.cv = TRUE,
                                                    committee.averaging = TRUE,
                                                    prob.mean.weight = TRUE,
                                                    VarImport = 0)
(Silene_ensemble_models_scores_EU <- get_evaluations(Silene_ensemble_models_EU)) #fiabilité 
TSS_Silene_EM <- Silene_ensemble_models_EU[[2]]
AUC_Silene_EM <- Silene_ensemble_models_EU[[3]]
write.csv(TSS_Silene_EM,"Silene_Emodel_EU_TSS.csv")
write.csv(AUC_Silene_EM,"Silene_Emodel_EU_AUC.csv")

###################################################################################################
### Projections ### (sur l'Europe)

Silene_models_EU<-get(load("SileneModelEU/SileneModelEU.1649228003.models.out"))
Silene_ensemble_models_EU<-get(load("SileneModelEU/SileneModelEU.1649228003ensemble.models.out"))

# crop des variables sélectionnées au niveau de l'Europe:
bioclim_map <- stack(list.files(path="2_CHELSA_bioclim/CHELSA_bioclim_1981_2010",pattern = ".tif",full.names = T))
slope<- raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")
Mask_EU<-rast("Shapes_Masks/shapeEurope/Mask_EU.tif")
Mask_EU<-terra::project(Mask_EU, "epsg:4326")
#couper le masque à 35° Est (réduire surface = réduire nombre de pixels = réduire temps de calculs):
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
bioclim_map<-crop(bioclim_map, Mask_EU)
slope <-crop(slope, extent(Mask_EU))

Climtopo<-stack(bioclim_map,slope)
ClimTopo_sub_new<-stack(subset(Climtopo,c("CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                                          "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")))

# Remplacer les valeurs dans les mers par des NA (diminue calculs du modèle et donc le temps de génération):
ForWater <- raster("Shapes_Masks/shapeEurope/ForNAWater.tif")
ForWater <- terra::crop(ForWater,Mask_EU)
ClimTopo_sub_new <- mask(ClimTopo_sub_new,ForWater)

## Projections actuelles ##

Silene_models_proj_current_EU <- BIOMOD_Projection(modeling.output = Silene_models_EU,
                                                   new.env = stack(ClimTopo_sub_new),
                                                   proj.name = "Silene_Current_EU",
                                                   binary.meth = "TSS",
                                                   selected.models="all")
#plot(Silene_models_proj_current_EU)

#BIOMOD_EnsembleProjection
Silene_ensemble_models_proj_current_EU <- BIOMOD_EnsembleForecasting(EM.output = Silene_ensemble_models_EU,
                                                                     projection.output = Silene_models_proj_current_EU,
                                                                     binary.meth = "TSS",
                                                                     selected.models="all")
plot(Silene_ensemble_models_proj_current_EU)
# Save:
cartes_Silene <- stack("SileneModelEU/proj_Silene_Current_EU/proj_Silene_Current_EU_SileneModelEU_ensemble.grd")
writeRaster(cartes_Silene$SileneModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"SileneModelEU/Silene_Current_EU_weighted_mean.tif")
carte_SileneEU <- raster("SileneModelEU/Silene_Current_EU_weighted_mean.tif")
plot(carte_SileneEU)

## Projections futures #########################################################
# Revoir les variables
# Télécharger prédicteurs bioclim pour le future (ceux issus de la sélection de variables) :
# 2070-2100
# Pour les scénarios RCP2.6 (optimiste) et RCP8.5 (pessimiste)
# RCP: https://doi.org/10.1016/j.gloenvcha.2015.01.004
# GCM = MPI-ESM
bio6_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio6_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio10_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio10_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio10_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio10_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio12_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio12_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio15_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio15_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")

#Couper les cartes sur la zone d'étude:
bio6_RCP26<-crop(bio6_RCP26, Mask_EU)
bio6_RCP85<-crop(bio6_RCP85, Mask_EU)
bio10_RCP26<-crop(bio10_RCP26, Mask_EU)
bio10_RCP85<-crop(bio10_RCP85, Mask_EU)
bio12_RCP26<-crop(bio12_RCP26, Mask_EU)
bio12_RCP85<-crop(bio12_RCP85, Mask_EU)
bio15_RCP26<-crop(bio15_RCP26, Mask_EU)
bio15_RCP85<-crop(bio15_RCP85, Mask_EU)

# Stack des variables (+ reprendre la pente):
BioClim_RCP26_MPI<-stack(bio6_RCP26,bio10_RCP26,bio12_RCP26,bio15_RCP26,slope)
BioClim_RCP85_MPI<-stack(bio6_RCP85,bio10_RCP85,bio12_RCP85,bio15_RCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26_MPI <- mask(BioClim_RCP26_MPI,ForWater)
BioClim_RCP85_MPI <- mask(BioClim_RCP85_MPI,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_RCP26_MPI)<-c("CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                            "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

# Code à reproduire pour chaque scénarios (GCM x RCP x année):

# GCM = MPIESM

# Projection future optimiste
new_env_26_MPI<-stack(BioClim_RCP26_MPI)
Silene_models_proj_2100_RCP26_MPI <-BIOMOD_Projection(modeling.output = Silene_models_EU,
                                                      new.env = new_env_26_MPI,
                                                      proj.name = "Silene_2100_RCP26_MPI",
                                                      binary.meth = "TSS",
                                                      selected.models = 'all')

Silene_ensemble_models_proj_2100_RCP26_MPI <- BIOMOD_EnsembleForecasting(EM.output = Silene_ensemble_models_EU,
                                                                         projection.output = Silene_models_proj_2100_RCP26_MPI,
                                                                         binary.meth = "TSS",
                                                                         selected.models = 'all')
plot(Silene_ensemble_models_proj_2100_RCP26_MPI)
#Save:
cartes_Silene_RCP26_MPI <- stack("SileneModelEU/proj_Silene_2100_RCP26_MPI/proj_Silene_2100_RCP26_MPI_SileneModelEU_ensemble.grd")
writeRaster(cartes_Silene_RCP26_MPI$SileneModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"SileneModelEU/Silene_2100_RCP26_MPI_weighted_mean.tif")
carte_SileneEU_RCP26_MPI <- raster("SileneModelEU/Silene_2100_RCP26_MPI_weighted_mean.tif")
plot(carte_SileneEU_RCP26_MPI)

# Projection future pessimiste
new_env_85<-stack(BioClim_RCP85_MPI)
names(new_env_85)<-c("CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                     "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

Silene_models_proj_2100_RCP85_MPI <-BIOMOD_Projection(modeling.output = Silene_models_EU,
                                                  new.env = new_env_85,
                                                  proj.name = "Silene_2100_RCP85_MPI",
                                                  binary.meth = "TSS",
                                                  selected.models = 'all')

Silene_ensemble_models_proj_2100_RCP85_MPI <- BIOMOD_EnsembleForecasting(EM.output = Silene_ensemble_models_EU,
                                                                     projection.output = Silene_models_proj_2100_RCP85_MPI,
                                                                     binary.meth = "TSS",
                                                                     selected.models = 'all')
plot(Silene_ensemble_models_proj_2100_RCP85_MPI)
#Save:
cartes_Silene_RCP85_MPI <- stack("SileneModelEU/proj_Silene_2100_RCP85_MPI/proj_Silene_2100_RCP85_MPI_SileneModelEU_ensemble.grd")
writeRaster(cartes_Silene_RCP85_MPI$SileneModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"SileneModelEU/Silene_2100_RCP85_MPI_weighted_mean.tif")
carte_SileneEU_RCP85_MPI <- raster("SileneModelEU/Silene_2100_RCP85_MPI_weighted_mean.tif")
plot(carte_SileneEU_RCP85_MPI)

# GCM = UKESM

bio6_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio6_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio6_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio6_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio10_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio10_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio10_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio10_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio12_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio12_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio12_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio12_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio15_UKRCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio15_UKRCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
slope<-raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

#Couper les cartes sur la zone d'étude:
bio6_UKRCP26<-crop(bio6_UKRCP26, Mask_EU)
bio6_UKRCP85<-crop(bio6_UKRCP85, Mask_EU)
bio10_UKRCP26<-crop(bio10_UKRCP26, Mask_EU)
bio10_UKRCP85<-crop(bio10_UKRCP85, Mask_EU)
bio12_UKRCP26<-crop(bio12_UKRCP26, Mask_EU)
bio12_UKRCP85<-crop(bio12_UKRCP85, Mask_EU)
bio15_UKRCP26<-crop(bio15_UKRCP26, Mask_EU)
bio15_UKRCP85<-crop(bio15_UKRCP85, Mask_EU)
slope<-crop(slope, Mask_EU)

# Stack des variables (+ reprendre la pente):
BioClim_RCP26_UKESM<-stack(bio6_UKRCP26,bio10_UKRCP26,bio12_UKRCP26,bio15_UKRCP26,slope)
BioClim_RCP85_UKESM<-stack(bio6_UKRCP85,bio10_UKRCP85,bio12_UKRCP85,bio15_UKRCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26_UKESM <- mask(BioClim_RCP26_UKESM,ForWater)
BioClim_RCP85_UKESM <- mask(BioClim_RCP85_UKESM,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_RCP26_UKESM)<-c("CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                              "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

# Projection future optimiste
new_env_26_UK<-stack(BioClim_RCP26_UKESM)
Silene_models_proj_2100_RCP26_UKESM <-BIOMOD_Projection(modeling.output = Silene_models_EU,
                                                             new.env = new_env_26_UK,
                                                             proj.name = "Silene_2100_RCP26_UKESM",
                                                             binary.meth = "TSS",
                                                             selected.models = 'all')

Silene_ensemble_models_proj_2100_RCP26_UKESM <- BIOMOD_EnsembleForecasting(EM.output = Silene_ensemble_models_EU,
                                                                                projection.output = Silene_models_proj_2100_RCP26_UKESM,
                                                                                binary.meth = "TSS",
                                                                                selected.models = 'all')
plot(Silene_ensemble_models_proj_2100_RCP26_UKESM)
#Save:
cartes_Silene_RCP26_UKESM <- stack("SileneModelEU/proj_Silene_2100_RCP26_UKESM/proj_Silene_2100_RCP26_UKESM_SileneModelEU_ensemble.grd")
writeRaster(cartes_Silene_RCP26_UKESM$SileneModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"SileneModelEU/Silene_2100_RCP26_UKESM_weighted_mean.tif")
carte_SileneEU_RCP26_UKESM <- raster("SileneModelEU/Silene_2100_RCP26_UKESM_weighted_mean.tif")
plot(carte_SileneEU_RCP26_UKESM)

# Projection future pessimiste
new_env_85_UKESM<-stack(BioClim_RCP85_UKESM)
names(new_env_85_UKESM)<-c("CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                           "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")
Silene_models_proj_2100_RCP85_UKESM <-BIOMOD_Projection(modeling.output = Silene_models_EU,
                                                             new.env = new_env_85_UKESM,
                                                             proj.name = "Silene_2100_RCP85_UKESM",
                                                             binary.meth = "TSS",
                                                             selected.models = 'all')

Silene_ensemble_models_proj_2100_RCP85_UKESM <- BIOMOD_EnsembleForecasting(EM.output = Silene_ensemble_models_EU,
                                                                                projection.output = Silene_models_proj_2100_RCP85_UKESM,
                                                                                binary.meth = "TSS",
                                                                                selected.models = 'all')
plot(Silene_ensemble_models_proj_2100_RCP85_UKESM)
#Save:
cartes_Silene_RCP85_UKESM <- stack("SileneModelEU/proj_Silene_2100_RCP85_UKESM/proj_Silene_2100_RCP85_UKESM_SileneModelEU_ensemble.grd")
writeRaster(cartes_Silene_RCP85_UKESM$SileneModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"SileneModelEU/Silene_2100_RCP85_UKESM_weighted_mean.tif")
carte_SileneEU_RCP85_UKESM <- raster("SileneModelEU/Silene_2100_RCP85_UKESM_weighted_mean.tif")
plot(carte_SileneEU_RCP85_UKESM)
#removeTmpFiles()

## Projections passées #########################################################

# Télécharger variables bioclim pour le LGM (PMIP3, CCSM4):
bio6_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_06.tif")
bio10_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_10.tif")
bio12_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_12.tif")
bio15_LGM<-raster("2_CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_15.tif")

# + recalculer la pente (fait dans SAGAGis):
slope_LGM<-raster("3_VariablesEnvi_horsBioclim/lgm/slope_lgm.tif")

# Couper les cartes sur la zone d'étude:
bio6_LGM<-crop(bio6_LGM, Mask_EU)
bio10_LGM<-crop(bio10_LGM, Mask_EU)
bio12_LGM<-crop(bio12_LGM, Mask_EU)
bio15_LGM<-crop(bio15_LGM, Mask_EU)
slope_LGM<-crop(slope_LGM, Mask_EU)
slope_LGM <- resample(slope_LGM,bio15_LGM)

# Convertir les unités (pour que ce soit les mêmes que pour les proj current et futures):
bio6_LGM<-(bio6_LGM/10)-273.15
bio10_LGM<-(bio10_LGM/10)-273.15
bio12_LGM<-bio12_LGM/10
bio15_LGM<-bio15_LGM/10

# Stack des variables (+ reprendre la pente):
BioClim_LGM<-stack(bio6_LGM,bio10_LGM,bio12_LGM,bio15_LGM,slope_LGM)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_LGM)<-c("CHELSA_bio6_19812010_V21","CHELSA_bio10_19812010_V21","CHELSA_bio12_19812010_V21",
                      "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

# Projection passée (LGM: -21 000 BP):
new_env_LGM<-stack(BioClim_LGM)
Silene_models_proj_LGM <-BIOMOD_Projection(modeling.output = Silene_models_EU,
                                           new.env = new_env_LGM,
                                           proj.name = "Silene_LGM",
                                           binary.meth = "TSS",
                                           selected.models = 'all')

Silene_ensemble_models_proj_LGM <- BIOMOD_EnsembleForecasting(EM.output = Silene_ensemble_models_EU,
                                                              projection.output = Silene_models_proj_LGM,
                                                              binary.meth = "TSS",
                                                              selected.models = 'all')
plot(Silene_ensemble_models_proj_LGM)
#Save:
cartes_Silene_LGM <- stack("SileneModelEU/proj_Silene_LGM/proj_Silene_LGM_SileneModelEU_ensemble.grd")
writeRaster(cartes_Silene_LGM$SileneModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"SileneModelEU/Silene_LGM_weighted_mean.tif")
carte_SileneEU_LGM <- raster("SileneModelEU/Silene_LGM_weighted_mean.tif")
plot(carte_SileneEU_LGM)

##############################################################################################################################################
### Comparer les distributions ###

# Charger les projections (carteWM):
Proj_Current<-raster("SileneModelEU/Silene_Current_EU_weighted_mean.tif")
Proj_LGM_Silene<-raster("SileneModelEU/Silene_LGM_weighted_mean.tif")
Proj_MPI26<-raster("SileneModelEU/Silene_2100_RCP26_MPI_weighted_mean.tif")
Proj_MPI85<-raster("SileneModelEU/Silene_2100_RCP85_MPI_weighted_mean.tif")
Proj_UK26<-raster("SileneModelEU/Silene_2100_RCP26_UKESM_weighted_mean.tif")
Proj_UK85<-raster("SileneModelEU/Silene_2100_RCP85_UKESM_weighted_mean.tif")

# Pour Proj_LGM: retirer les zones sous glace comme étant favorable:
Glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # charger carte présence de glaciers
# Couper carte glacier sur l'EU:
Mask_EU<-rast("Y_Shapes_Masks/shapeEurope/ForNaWater.tif")
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
Glacier<-crop(Glacier,Mask_EU)
glacier<-resample(Glacier,Proj_LGM_Silene,method='ngb')
values(glacier)[is.na(values(glacier))] = 0 #sinon ça vire là où il y a des NA et donc quand il n'y a pas de glacier
glacierAbsence <- 1-glacier # pour mettre la favorabilité à 0 quand présence de glacier
Proj_LGM_HG_Silene<-Proj_LGM_Silene*glacierAbsence
plot(Proj_LGM_HG_Silene) 
writeRaster(Proj_LGM_HG_Silene,"6_Cartes_Projections/Silene_LGM_WM_2.tif")

## Binariser les cartes de projection (en sélectionnant le treshold qui maximise le TSS):
# Treshold ici = 478
Silene_bin_current<-BinaryTransformation(Proj_Current, 478)
Silene_bin_2100_MPI26<-BinaryTransformation(Proj_MPI26, 478)
Silene_bin_2100_MPI85<-BinaryTransformation(Proj_MPI85,478 )
Silene_bin_2100_UK26<-BinaryTransformation(Proj_UK26, 478)
Silene_bin_2100_UK85<-BinaryTransformation(Proj_UK85, 478)
Silene_bin_LGM<-BinaryTransformation(Proj_LGM, 478)
plot(Silene_bin_current)
plot(Silene_bin_2100_MPI26)
plot(Silene_bin_2100_MPI85)
plot(Silene_bin_2100_UK26)
plot(Silene_bin_2100_UK85)
plot(Silene_bin_LGM)
writeRaster(Silene_bin_current,"SileneModelEU/Silene_bin_current.tif")
writeRaster(Silene_bin_2100_MPI26,"SileneModelEU/Silene_bin_2100_MPI26.tif")
writeRaster(Silene_bin_2100_MPI85,"SileneModelEU/Silene_bin_2100_MPI85.tif")
writeRaster(Silene_bin_2100_UK26,"SileneModelEU/Silene_bin_2100_UK26.tif")
writeRaster(Silene_bin_2100_UK85,"SileneModelEU/Silene_bin_2100_UK85.tif")
writeRaster(Silene_bin_LGM,"SileneModelEU/Silene_bin_LGM.tif")

## Rangeshifts:

# Load cartes binarisées:
Silene_bin_current<-raster("SileneModelEU/Silene_bin_current.tif")
Silene_bin_2100_MPI26<-raster("SileneModelEU/Silene_bin_2100_MPI26.tif")
Silene_bin_2100_MPI85<-raster("SileneModelEU/Silene_bin_2100_MPI85.tif")
Silene_bin_2100_UK26<-raster("SileneModelEU/Silene_bin_2100_UK26.tif")
Silene_bin_2100_UK85<-raster("SileneModelEU/Silene_bin_2100_UK85.tif")
Silene_bin_LGM_2<-raster("SileneModelEU/Silene_bin_LGM_2.tif")

# Retirer glacier de Silene_bin_LGM (car surprédiction de sa présence ds et autour du glacier):
glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # carte de présence/absence du glacier
glacier<-crop(glacier, Mask_EU)
values(glacier)[is.na(values(glacier))]=0 #transformer les na en 0
Silene_bin_LGM_HorsGlacier<-Silene_bin_LGM-glacier
values(Silene_bin_LGM_HorsGlacier)[values(Silene_bin_LGM_HorsGlacier)<= 0] = 0 #ne garder que les présences (les valeurs=1, les autres =0)
writeRaster(Silene_bin_LGM_HorsGlacier,"Silene_bin_LGM_2.tif")

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
Shift_LGM_current<-BIOMOD_RangeSize(Silene_bin_current,Silene_bin_LGM_HorsGlacier)
Shift_LGM_current<-Shift_LGM_current$Compt.By.Models
Shift_LGM_current<-as.data.frame(Shift_LGM_current)
#plot(Shift_LGM_current_$Diff.By.Pixel) # 0 = pixels non-occupé, 1 = pixels occupés slmt ds la nouvelle période
# -1 = pixels occupés aux 2 périodes, -2 = pixels perdus entre current et nouvelle période
diff<-2*Silene_bin_current + Silene_bin_LGM_2
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Silene_Diff_LGM.tif")

# Current-MPI26
Shift_current_2100_MPI26 <- BIOMOD_RangeSize(Silene_bin_current, Silene_bin_2100_MPI26)
Shift_current_2100_MPI26<-Shift_current_2100_MPI26$Compt.By.Models
Shift_current_2100_MPI26<-as.data.frame(Shift_current_2100_MPI26)
# Pour le plot (plus lisible que plot diff.by.pixel):
diff<-2*Silene_bin_current + Silene_bin_2100_MPI26 
plot(diff)  # diff = 1: présence qu'au future (gain, en orange), diff=2: présence qu'au présent (loss, en vert clair), diff=3: présence aux 2 périodes (stable, en vert foncé), 0 = pixels non-occupés
writeRaster(diff, "6_Cartes_Projections/Silene_Diff_MPI26.tif")

# Current-MPI85
Shift_current_2100_MPI85 <- BIOMOD_RangeSize(Silene_bin_current, Silene_bin_2100_MPI85)
Shift_current_2100_MPI85<-Shift_current_2100_MPI85$Compt.By.Models
Shift_current_2100_MPI85<-as.data.frame(Shift_current_2100_MPI85)
diff<-2*Silene_bin_current + Silene_bin_2100_MPI85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Silene_Diff_MPI85.tif")

# Current-UK26
Shift_current_2100_UK26 <- BIOMOD_RangeSize(Silene_bin_current, Silene_bin_2100_UK26)
Shift_current_2100_UK26<-Shift_current_2100_UK26$Compt.By.Models
Shift_current_2100_UK26<-as.data.frame(Shift_current_2100_UK26)
diff<-2*Silene_bin_current + Silene_bin_2100_UK26 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Silene_Diff_UK26.tif")

# Current-UK85
Shift_current_2100_UK85 <- BIOMOD_RangeSize(Silene_bin_current, Silene_bin_2100_UK85)
Shift_current_2100_UK85<-Shift_current_2100_UK85$Compt.By.Models
Shift_current_2100_UK85<-as.data.frame(Shift_current_2100_UK85)
diff<-2*Silene_bin_current + Silene_bin_2100_UK85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Silene_Diff_UK85.tif")

# Fusionner les dataframe:
RangeShifts<-bind_rows(Shift_LGM_current,Shift_current_2100_MPI26,Shift_current_2100_MPI85,
                       Shift_current_2100_UK26,Shift_current_2100_UK85,.id=NULL)
write.csv(RangeShifts,"Models_and_Resultats_horsCluster/Silene_RangeShifts_2.csv")

## Shift altitudinal:

altPres <- raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif")
altPast<- raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_dem_-200_V1.0.tif")
# Extraire coordonnées au niveau des présences pour chaque période:
# Current
coordPres <- coordinates(Silene_bin_current)
valPres <- values(Silene_bin_current)
coordPres <- coordPres[!is.na(valPres) & valPres == 1,]
altiPres <- extract(altPres,coordPres)
MedAltPres<-median(altiPres)
# LGM
coordLGM <- coordinates(Silene_bin_LGM_HorsGlacier)
valLGM <- values(Silene_bin_LGM_HorsGlacier)
coordLGM <- coordLGM[!is.na(valLGM) & valLGM == 1,]
altiLGM <- extract(altPast,coordLGM)
MedAltLGM<-median(altiLGM)
# MPI26
coordMPI26 <- coordinates(Silene_bin_2100_MPI26)
valMPI26 <- values(Silene_bin_2100_MPI26)
coordMPI26 <- coordMPI26[!is.na(valMPI26) & valMPI26 == 1,]
altiMPI26 <- extract(altPres,coordMPI26)
MedAltMPI26<-median(altiMPI26)
# MPI85
coordMPI85 <- coordinates(Silene_bin_2100_MPI85)
valMPI85 <- values(Silene_bin_2100_MPI85)
coordMPI85 <- coordMPI85[!is.na(valMPI85) & valMPI85 == 1,]
altiMPI85 <- extract(altPres,coordMPI85)
MedAltMPI85<-median(altiMPI85)
# UK26
coordUK26 <- coordinates(Silene_bin_2100_UK26)
valUK26 <- values(Silene_bin_2100_UK26)
coordUK26 <- coordUK26[!is.na(valUK26) & valUK26 == 1,]
altiUK26 <- extract(altPres,coordUK26)
MedAltUK26<-median(altiUK26)
# UK85
coordUK85 <- coordinates(Silene_bin_2100_UK85)
valUK85<- values(Silene_bin_2100_UK85)
coordUK85 <- coordUK85[!is.na(valUK85) & valUK85 == 1,]
altiUK85 <- extract(altPres,coordUK85)
MedAltUK85<-median(altiUK85)

# Mettre les médianes en graphique:
pourBoxplotPres <- cbind.data.frame(Altitude = altiPres, Periode = "Present")
pourBoxplotLGM <- cbind.data.frame(Altitude = altiLGM, Periode = "LGM")
pourBoxplotMPI26 <- cbind.data.frame(Altitude = altiMPI26, Periode = "MPIESM-RCP26")
pourBoxplotMPI85 <- cbind.data.frame(Altitude = altiMPI85, Periode = "MPIESM-RCP85")
pourBoxplotUK26 <- cbind.data.frame(Altitude = altiUK26, Periode = "UKESM-RCP26")
pourBoxplotUK85 <- cbind.data.frame(Altitude = altiUK85, Periode = "UKESM-RCP85")

pourBoxplot <- rbind(pourBoxplotPres,pourBoxplotLGM,pourBoxplotMPI26,pourBoxplotMPI85,
                     pourBoxplotUK26,pourBoxplotUK85)
# GGPlot2
pourBoxplot$Periode = as.factor(pourBoxplot$Periode)
p <- ggplot(data = pourBoxplot,aes(x=Periode,y=Altitude)) + geom_boxplot() + labs(title="Silene acaulis") +theme_classic() 
p


##############################################################################
### Calculer vitesses de migration

# A. Altitudinal

# RangeSchift alt (différences entre médianes; si diff<0: sp est descendue en alt par rapport au présent, si diff>0: sp est montée en alt par rapport au présent):
DiffAlt_Pres_LGM<-MedAltLGM-MedAltPres
DiffAlt_Pres_LGM # -470,7 m
DiffAlt_Pres_MPI26<-MedAltMPI26-MedAltPres
DiffAlt_Pres_MPI26 # 77,6 m
DiffAlt_Pres_MPI85<-MedAltMPI85-MedAltPres
DiffAlt_Pres_MPI85 # -33 m
DiffAlt_Pres_UK26<-MedAltUK26-MedAltPres
DiffAlt_Pres_UK26 # -12 m
DiffAlt_Pres_UK85<-MedAltUK85-MedAltPres
DiffAlt_Pres_UK85 # 117,6 m

# Diviser les rangeshift altitudinal par deltaT:
# Entre LGM et Current (1995-(-21000)=22995 ans)
MigrationAlt_Pres_LGM<-DiffAlt_Pres_LGM/22995
MigrationAlt_Pres_LGM # 0,02 m/an 

# Entre Current et future (2085-1995= 90 ans) selon scénario:
MigAlt_Pres_MPI26<-DiffAlt_Pres_MPI26/90
MigAlt_Pres_MPI26 # 0,86 m/an
MigAlt_Pres_MPI85<-DiffAlt_Pres_MPI85/90
MigAlt_Pres_MPI85 # -0,37 m/an
MigAlt_Pres_UK26<-DiffAlt_Pres_UK26/90
MigAlt_Pres_UK26 # -0,13 m/an
MigAlt_Pres_UK85<-DiffAlt_Pres_UK85/90
MigAlt_Pres_UK85 # 1,3 m/an

# B. latitudinal/longitudinal

# Chargement des coordonnées des centroides (calculé dans QGIS):
Centroides<-read.csv("ShiftLat_Silene_2.csv",header=F,row.names=1, sep=";",dec=",")
# Calcul des distances entre les centroides des différents scenarios:
Distances_Silene<-pointDistance(Centroides, lonlat=T,allpairs=T)
