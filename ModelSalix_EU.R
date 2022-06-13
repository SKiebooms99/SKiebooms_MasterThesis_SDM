### Modèle Salix Europe ###

setwd("Z:/projects-unil/GEN4MIG/Solene")

library(biomod2)
library(raster)
library(foreign)
library(dismo)
library(ggplot2)
library(reshape2)
library(terra)

## télécharger les données (pour l'Europe) ##

# Espèce:
SalixPV<-read.table("Occurences_data/Occ_Salix_retusa/SalixPV.csv",header = TRUE, sep = ",", row.names=1)
SalixPV<-SalixPV[,c("decimalLon","decimalLat")] # 1228 points
# Supprimer les occurences se trouvant à une altitude aberrante (celles se trouvant en dehors du range altitudinal de l'sp):
altitude<-raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif") # doi:10.1038/sdata.2018.40.
AltPoints<-extract(altitude,SalixPV)
#View(AltPoints)
SalixPV <- SalixPV[AltPoints >= 1100,] #On retire les points en dessous de 1100m -> restent 1220

# Variables environnementales:

# variables bioclimatiques (actuelles, à l'échelle de l'Europe, résolution 1km):
bioclim_map <- stack(list.files(path="CHELSA_bioclim/CHELSA_bioclim_1981_2010",pattern = ".tif",full.names = T))

# Autres variables bio:

# Srad
srad <- stack(paste0("VariablesEnvi/CHELSA_srad/CHELSA_srad_",1:12,".tif"))
Srad<-sum(srad,na.rm = T)
plot(Srad)

# Degrés-jour: gdd0: heat sum of all days above the treshold temperature (0°C) accumulated over 1 year.
gdd0<-raster("VariablesEnvi/CHELSA_gdd0_1981-2010_V.2.1.tif")

# Indice d'humidité
#rh (humidité relative en %): daily mean near surface relative humidity averaged over 1 month
# rh<-stack(paste0("VariablesEnvi/CHELSA_rh/CHELSA_rh_",1:12,"_1979-2013_V1.0.tif"))
rh1<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_1_1979-2013_V1.0.tif")
rh2<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_2_1979-2013_V1.0.tif")
rh3<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_3_1979-2013_V1.0.tif")
rh4<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_4_1979-2013_V1.0.tif")
rh5<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_5_1979-2013_V1.0.tif")
rh6<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_6_1979-2013_V1.0.tif")
rh7<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_7_1979-2013_V1.0.tif")
rh8<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_8_1979-2013_V1.0.tif")
rh9<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_9_1979-2013_V1.0.tif")
rh10<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_10_1979-2013_V1.0.tif")
rh11<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_11_1979-2013_V1.0.tif")
rh12<-raster("VariablesEnvi/CHELSA_rh/CHELSA_rh_12_1979-2013_V1.0.tif")
rh<-mean(rh1,rh2,rh3,rh4,rh5,rh6,rh7,rh8,rh9,rh10,rh11,rh12)

#swb (mean summer (juin, juillet et aout) water balance): was calculated by subtracting mean monthly potential evapotranspiration (PET) from mean monthly precipitation
swb<-raster("VariablesEnvi/CHELSA_swb_1981-2010_V.2.1.tif")

# variables topographiques (à l'échelle européenne, résolution 1km: Amatulli et al., 2018: DOI: doi:10.1038/sdata.2018.40)
slope<- raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")
topo<- raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/tpi_1KMmd_GMTEDmd.tif") #TPI
aspectSin<- raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/aspectsine_1KMmd_GMTEDmd.tif") #sin: axe nord-sud
aspectCos<-raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/aspectcosine_1KMmd_GMTEDmd.tif") # cos: axe ouest-est
#northness<- raster("VariablesTopo_EU_Amatulli2018/aspect_northness_1KMmd_GMTEDmd.tif")
#eastness<- raster("VariablesTopo_EU_Amatulli2018/aspect_eastness_1KMmn_GMTEDmd.tif")
#courbeP <- raster("VariablesTopo_EU_Amatulli2018/curvP_1KMmd_GMTEDmd.tif")
#courbeT<- raster("VariablesTopo_EU_Amatulli2018/curvT_1KMmn_GMTEDmd.tif")

# Uniformiser la taille des rasters:
  
Mask_EU<-rast("shapeEurope/Mask_EU.tif")
Mask_EU<-terra::project(Mask_EU, "epsg:4326")
#couper le masque à 35° Est (réduire surface = réduire nombre de pixels = réduire temps de calculs):
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)

bioclim_map<-crop(bioclim_map, Mask_EU)
slope <-crop(slope, extent(Mask_EU))
topo<-crop(topo, Mask_EU)
#northness<-crop(northness, Mask_EU)
#eastness<-terra::crop(eastness, Mask_EU)
#courbeP<-terra::crop(courbeP, Mask_EU)
#courbeT<-terra::crop(courbeT, Mask_EU)
#soil<-terra::crop(soil, Mask_EU)
gdd0<-crop(gdd0, Mask_EU)
Srad<-crop(Srad, Mask_EU)
rh<-crop(rh, Mask_EU)
swb<-crop(swb, Mask_EU)
aspectSin<-raster::crop(aspectSin, Mask_EU)
aspectCos<-raster::crop(aspectCos, Mask_EU)
  
# Stack des variables

Climtopo<-stack(slope,topo,aspectSin,aspectCos,bioclim_map,rh,Srad,gdd0,swb)

#############################################
## Sélection des prédicteurs environnementaux (variables explicatives parmi Climtopo)

# Test de corrélation de Pearson:
corr<-layerStats(Climtopo,stat="pearson",na.rm=T)
write.csv(corr$`pearson correlation coefficient`, file="matriceCorr_ClimBioTopo2.csv")
MatriceCorr_ClimTopo<-read.csv("matriceCorr_ClimBioTopo2.csv")
# Visualiser via un dendrogramme (clustering):
corr<-corr$`pearson correlation coefficient`
dendro <- hclust(as.dist(1-abs(corr)),method = "ward.D2")
plot(dendro)
abline(h=0.35,col="red",lty=2) # Treshold de corrélation < 0,65
NC <- abs(corr)
NC[NC<0.65]=NA # variables non-corrélées mises en évidence dans la matrice de corrélation

# Continuer avec le subset de variables sélectionnées (slope, srad, bio2,bio6, bio12, bio15, swb)
ClimTopo_sub<-stack(subset(Climtopo,c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                                      "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")))

#################################################################################
## Modélisation ##

## Sélection des prédicteurs environnementaux (variables explicatives parmis Climtopo)
# Partir de la sélection par test de corrélation de Pearson 
# et éventuellement ajuster en fonction des courbes de densité Bp - Presences

# Remplacer les valeurs dans les mers par des NA (diminue calculs du modèle et donc le temps de génération):
ForWater <- raster("shapeEurope/ForNAWater.tif")
ForWater <- crop(ForWater,Mask_EU)
ClimTopo_sub <- mask(ClimTopo_sub,ForWater) # -> saved as "ClimTopo_sub_Current"

# Générer des background points (Bp) pour l'Europe (et juste dans les biomes comprenant la majorité des points d'occurence)
biome <- raster("shapeEurope/GB_EBR_map_BiomesUE.tif")
crs(biome)
Occ_SalixProj <- SpatialPointsDataFrame(coords = SalixPV,data = SalixPV,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
Occ_SalixProj <- spTransform(Occ_SalixProj,CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
Occ_SalixProj <- as.data.frame(Occ_SalixProj)[,-c(1:2)]
biomePoints <- extract(biome, Occ_SalixProj)
BiomeToKeep <- as.numeric(names(table(biomePoints))[table(biomePoints)>= (5*nrow(Occ_SalixProj)/100)]) # on excluera les biomes contenant moins de 5% des occurences
BiomeToKeep

biomeToSample <- rasterToPolygons(biome,fun=function(x){x==1}) #sélectionner le(s) biome(s) à garder
Points <- extract(ClimTopo_sub, SalixPV)
Bp <- spsample(biomeToSample, 20000, type="random")
Bp <- spTransform(Bp,CRS("+proj=longlat +datum=WGS84 +no_defs "))
Bp <- as.data.frame(Bp)
BpEnv <- raster::extract(ClimTopo_sub,Bp)
Bp <- na.omit(cbind(Bp,BpEnv)) #retirer les NA dans les lignes
#CheckPointInTheWater <- extract(subset(ClimTopo_sub,1),Bp)
#Bp <- Bp[!is.na(CheckPointInTheWater),]
Bp <- Bp[1:10000,1:2]
write.csv(Bp,"BpCoordinates_Salix_EU2.csv")


#Fusionner les données de présences et les background points:
Bp<-read.csv("BpCoordinates_Salix_EU2.csv", row.names = 1)
colnames(SalixPV) = c("x","y")
Occ_Salix <- gridSample(SalixPV,ClimTopo_sub) # Ne garder qu'une occurence max par pixel -> 1208
Data_Salix<-rbind(Occ_Salix,Bp)
write.csv(Data_Salix, "Data_Salix_Occ_Bp_ClimTopo_EU2.csv")

###
#Graphique: densité des points (présence et Bp) pour chaque variable pour l'sp considérée
env <- raster::extract(ClimTopo_sub,Data_Salix)
resp <- c(rep("presence", nrow(Occ_Salix)), rep("background",10000))
tab <- cbind.data.frame(resp=resp,env)
tab$resp=as.factor(tab$resp)
tab <- melt(tab,id="resp")
rm(resp)
png("DensityPlotsVariableEnv_Salix_EU2.png",width= 5000,height = 8000,res=300)
p <- ggplot(data=tab,aes(x=value,color=resp,fill=resp)) + geom_density(alpha=0.1) +facet_wrap(~variable,scale= "free")
p
dev.off()
###

write.table(resp,"resp_Salix.txt") # pour le cluster de calibration

# Formater les données pour utiliser le package biomod2:
env <- raster::extract(ClimTopo_sub,Data_Salix)
resp <- c(rep(1, nrow(Occ_Salix)), rep(NA,10000))

Data_Salix_EU<-BIOMOD_FormatingData(resp.var = resp,
                                    expl.var = env,
                                    resp.xy = Data_Salix,
                                    resp.name = "SalixModelEU",
                                    PA.nb.rep = 1,
                                    PA.nb.absences = 10000,
                                    PA.strategy = 'random')

Salix_opt<-BIOMOD_ModelingOptions(GBM = list(n.trees = 5000),GAM = list(algo = "GAM_mgcv"))

Salix_models_EU<-BIOMOD_Modeling(data=Data_Salix_EU, models=c("GLM","GBM","GAM"),
                                 models.options = Salix_opt,NbRunEval = 10,DataSplit = 70,
                                 models.eval.meth = c("TSS", "ROC"),
                                 VarImport = 10,SaveObj=TRUE,do.full.models = FALSE,
                                 Prevalence = 0.5)

# Vérifier importance des variables:
Salix_models_var_import_ClimTopo <- get_variables_importance(Salix_models_EU) 
apply(Salix_models_var_import_ClimTopo, c(1, 2), mean) #The higher the score, the more important the variable
# Response curves
Salix_glm <- BIOMOD_LoadModels(Salix_models_EU, models = 'GLM')
Salix_gbm <- BIOMOD_LoadModels(Salix_models_EU, models = 'GBM')
Salix_gam <- BIOMOD_LoadModels(Salix_models_EU, models = 'GAM')

glm_eval_strip <- biomod2::response.plot2(models = Salix_glm,
                                          Data = get_formal_data(Salix_models_EU,
                                                                 'expl.var'),
                                          show.variables = get_formal_data(Salix_models_EU,
                                                                           'expl.var.names'),
                                          do.bivariate = FALSE,
                                          fixed.var.metric = 'mean',
                                          legend = FALSE,
                                          display_title = FALSE,
                                          data_species = get_formal_data(Salix_models_EU,'resp.var'))
gbm_eval_strip <-
  biomod2::response.plot2(models =Salix_gbm,
                          Data = get_formal_data(Salix_models_EU,
                                                 'expl.var'),
                          show.variables = get_formal_data(Salix_models_EU,
                                                           'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(Salix_models_EU,
                                                         'resp.var'))
gam_eval_strip <-
  biomod2::response.plot2(models = Salix_gam,
                          Data = get_formal_data(Salix_models_EU,
                                                 'expl.var'),
                          show.variables= get_formal_data(Salix_models_EU,
                                                          'expl.var.names'),
                          do.bivariate = FALSE,
                          fixed.var.metric = 'mean',
                          legend = FALSE,
                          display_title = FALSE,
                          data_species = get_formal_data(Salix_models_EU,
                                                         'resp.var'))

# Scores des modèles:
Salix_models_scores_EU <- get_evaluations(Salix_models_EU)
Salix_models_scores_EU
TSS_Salix <- as.data.frame(Salix_models_scores_EU["TSS",,,,])
AUC_Salix <- as.data.frame(Salix_models_scores_EU["ROC",,,,])
write.csv(TSS_Salix,"Salix_model_EU_TSS.csv")
write.csv(AUC_Salix,"Salix_model_EU_AUC.csv")

## Ensemble Modelling ##

Salix_ensemble_models_EU <-BIOMOD_EnsembleModeling(modeling.output = Salix_models_EU,
                                                   em.by = 'all',eval.metric = 'TSS',
                                                   eval.metric.quality.threshold = 0,
                                                   models.eval.meth = c("TSS", "ROC"),
                                                   prob.mean = FALSE,
                                                   prob.cv = TRUE,
                                                   committee.averaging = TRUE,
                                                   prob.mean.weight = TRUE,
                                                   VarImport = 0)
(Salix_ensemble_models_scores_EU <- get_evaluations(Salix_ensemble_models_EU)) #fiabilité 
TSS_Salix_EM <- Salix_ensemble_models_EU[[2]]
AUC_Salix_EM <- Salix_ensemble_models_EU[[3]]
write.csv(TSS_Salix_EM,"Salix_Emodel_EU_TSS.csv")
write.csv(AUC_Salix_EM,"Salix_Emodel_EU_AUC.csv")
##############################################################################################
### Projections ###

## Projections actuelles ##

Salix_models_proj_current_EU <- BIOMOD_Projection(modeling.output = Salix_models_EU,
                                                  new.env = stack(ClimTopo_sub),
                                                  proj.name = "Salix_Current_EU",
                                                  binary.meth = "TSS",
                                                  selected.models="all")
plot(Salix_models_proj_current_EU)

#BIOMOD_EnsembleProjection
Salix_ensemble_models_proj_current_EU <- BIOMOD_EnsembleForecasting(EM.output = Salix_ensemble_models_EU,
                                                                    projection.output = Salix_models_proj_current_EU,
                                                                    binary.meth = "TSS",
                                                                    selected.models="all")
plot(Salix_ensemble_models_proj_current_EU)
# Save:
cartes_Salix <- stack("OutputsCluster/SalixModelEU/proj_SalixCurrent/proj_SalixCurrent_SalixModelEU_ensemble.grd")
writeRaster(cartes_Salix$SalixModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"OutputsCluster/SalixModelEU/SalixCurrent_weighted_mean.tif")
carte_SalixEU <- raster("OutputsCluster/SalixModelEU/SalixCurrent_weighted_mean.tif")
plot(carte_SalixEU)

## Projections futures #########################################################

# Télécharger prédicteurs bioclim pour le future (ceux issus de la sélection de variables) :
# 2070-2100
# Pour les scénarios RCP2.6 (optimiste) et RCP8.5 (pessimiste)
# RCP: https://doi.org/10.1016/j.gloenvcha.2015.01.004
# GCM = MPI-ESM
bio2_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio2_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio2_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio2_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio6_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio6_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio12_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio12_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio15_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio15_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")

#Couper les cartes sur la zone d'étude:
bio2_RCP26<-crop(bio2_RCP26, Mask_EU)
bio2_RCP85<-crop(bio2_RCP85, Mask_EU)
bio6_RCP26<-crop(bio6_RCP26, Mask_EU)
bio6_RCP85<-crop(bio6_RCP85, Mask_EU)
bio12_RCP26<-crop(bio12_RCP26, Mask_EU)
bio12_RCP85<-crop(bio12_RCP85, Mask_EU)
bio15_RCP26<-crop(bio15_RCP26, Mask_EU)
bio15_RCP85<-crop(bio15_RCP85, Mask_EU)

# Stack des variables (+ reprendre la pente):
BioClim_RCP26<-stack(bio2_RCP26,bio6_RCP26,bio12_RCP26,bio15_RCP26,slope)
BioClim_RCP85<-stack(bio2_RCP85,bio6_RCP85,bio12_RCP85,bio15_RCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26 <- mask(BioClim_RCP26,ForWater)
BioClim_RCP85 <- mask(BioClim_RCP85,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_RCP26)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                        "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

# Code à reproduire pour chaque scénarios (GCM x RCP x année):

# Projection future optimiste
new_env<-stack(BioClim_RCP26)
Salix_models_proj_2100_RCP26 <-BIOMOD_Projection(modeling.output = Salix_models_EU,
                                                       new.env = new_env,
                                                       proj.name = "Salix_2100_RCP26",
                                                       binary.meth = "TSS",
                                                       selected.models = 'all')

Salix_ensemble_models_proj_2100_RCP26 <- BIOMOD_EnsembleForecasting(EM.output = Salix_ensemble_models_EU,
                                                                          projection.output = Salix_models_proj_2100_RCP26,
                                                                          binary.meth = "TSS",
                                                                          selected.models = 'all')
plot(Salix_ensemble_models_proj_2100_RCP26)
#Save:
cartes_Salix <- stack("OutputsCluster/SalixModelEU/proj_SalixMPIESM26/proj_SalixMPIESM26_SalixModelEU_ensemble.grd")
writeRaster(cartes_Salix$SalixModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"OutputsCluster/SalixModelEU/SalixMPIESM26_weighted_mean.tif")
carte_SalixEU <- raster("OutputsCluster/SalixModelEU/SalixMPIESM26_weighted_mean.tif")
plot(carte_SalixEU)

# Projection future pessimiste
new_env_85<-stack(BioClim_RCP85)
names(new_env_85)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                     "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")
Salix_models_proj_2100_RCP85 <-BIOMOD_Projection(modeling.output = Salix_models_EU,
                                                 new.env = new_env,
                                                 proj.name = "Salix_2100_RCP85",
                                                 binary.meth = "TSS",
                                                 selected.models = 'all')

Salix_ensemble_models_proj_2100_RCP85 <- BIOMOD_EnsembleForecasting(EM.output = Salix_ensemble_models_EU,
                                                                    projection.output = Salix_models_proj_2100_RCP85,
                                                                    binary.meth = "TSS",
                                                                    selected.models = 'all')
plot(Salix_ensemble_models_proj_2100_RCP85)
#Save:
cartes_Salix <- stack("OutputsCluster/SalixModelEU/proj_SalixMPIESM85/proj_SalixMPIESM85_SalixModelEU_ensemble.grd")
writeRaster(cartes_Salix$SalixModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"OutputsCluster/SalixModelEU/SalixMPIESM85_weighted_mean.tif")
carte_SalixEU <- raster("OutputsCluster/SalixModelEU/SalixMPIESM85_weighted_mean.tif")
plot(carte_SalixEU)

## Projections passées #########################################################

# BioClim_LGM<-raster("ClimTopo_sub_LGM.tif")

# Télécharger variables bioclim pour le LGM (PMIP3, CCSM4):
bio2_LGM<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_02.tif")
bio6_LGM<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_06.tif")
bio12_LGM<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_12.tif")
bio15_LGM<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_15.tif")
#BioClim_LGM<-stack(bio2_LGM,bio6_LGM,bio12_LGM,bio15_LGM)

# + recalculer la pente (fait dans SAGAGis):
slope_LGM<-raster("VariablesEnvi/lgm/slope_lgm.tif")

# Couper les cartes sur la zone d'étude:
bio2_LGM<-crop(bio2_LGM, Mask_EU)
bio6_LGM<-crop(bio6_LGM, Mask_EU)
bio12_LGM<-crop(bio12_LGM, Mask_EU)
bio15_LGM<-crop(bio15_LGM, Mask_EU)
slope_LGM<-crop(slope_LGM, Mask_EU)
slope_LGM <- resample(slope_LGM,bio15_LGM)

# Stack des variables (+ reprendre la pente):
BioClim_LGM<-stack(bio2_LGM,bio6_LGM,bio12_LGM,bio15_LGM,slope_LGM)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_LGM)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                      "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

# Projection passée (LGM: -21 000 BP):
new_env<-stack(BioClim_LGM)
Salix_models_proj_LGM <-BIOMOD_Projection(modeling.output = Salix_models_EU,
                                                 new.env = new_env,
                                                 proj.name = "Salix_LGM",
                                                 binary.meth = "TSS",
                                                 selected.models = 'all')

Salix_ensemble_models_proj_LGM <- BIOMOD_EnsembleForecasting(EM.output = Salix_ensemble_models_EU,
                                                                    projection.output = Salix_models_proj_LGM,
                                                                    binary.meth = "TSS",
                                                                    selected.models = 'all')
plot(Salix_ensemble_models_proj_LGM)
#Save:
cartes_Salix <- stack("OutputsCluster/SalixModelEU/proj_SalixLGM/proj_SalixLGM_SalixModelEU_ensemble.grd")
writeRaster(cartes_Salix$SalixModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"OutputsCluster/SalixModelEU/SalixLGM_weighted_mean.tif")
carte_SalixEU <- raster("OutputsCluster/SalixModelEU/SalixLGM_weighted_mean.tif")
plot(carte_SalixEU)

## GCM= UKESM

cartes_Salix <- stack("OutputsCluster/SalixModelEU/proj_SalixUKESM85/proj_SalixUKESM85_SalixModelEU_ensemble.grd")
writeRaster(cartes_Salix$SalixModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"OutputsCluster/SalixModelEU/SalixUKESM85_weighted_mean.tif")
carte_SalixEU <- raster("OutputsCluster/SalixModelEU/SalixUKESM85_weighted_mean.tif")
plot(carte_SalixEU)

cartes_Salix <- stack("OutputsCluster/SalixModelEU/proj_SalixUKESM26/proj_SalixUKESM26_SalixModelEU_ensemble.grd")
writeRaster(cartes_Salix$SalixModelEU_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData,"OutputsCluster/SalixModelEU/SalixUKESM26_weighted_mean.tif")
carte_SalixEU <- raster("OutputsCluster/SalixModelEU/SalixUKESM26_weighted_mean.tif")
plot(carte_SalixEU)

##############################################################################################################################################
### Comparer les distributions ###

# Charger les projections (carteWM):
Proj_Current<-raster("4_OutputsCluster/SalixModelEU/SalixCurrent_weighted_mean.tif")
Proj_LGM_Salix<-raster("4_OutputsCluster/SalixModelEU/SalixLGM_weighted_mean.tif")
Proj_MPI26<-raster("4_OutputsCluster/SalixModelEU/SalixMPIESM26_weighted_mean.tif")
Proj_MPI85<-raster("4_OutputsCluster/SalixModelEU/SalixMPIESM85_weighted_mean.tif")
Proj_UK26<-raster("4_OutputsCluster/SalixModelEU/SalixUKESM26_weighted_mean.tif")
Proj_UK85<-raster("4_OutputsCluster/SalixModelEU/SalixUKESM85_weighted_mean.tif")

# Pour Proj_LGM: retirer les zones sous glace comme étant favorable:
Glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # charger carte présence de glaciers
# Couper carte glacier sur l'EU:
Mask_EU<-rast("Y_Shapes_Masks/shapeEurope/ForNaWater.tif")
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
Glacier<-crop(Glacier,Mask_EU)
glacier<-resample(Glacier,Proj_LGM_Salix,method='ngb')
values(glacier)[is.na(values(glacier))] = 0 #sinon ça vire là où il y a des NA et donc quand il n'y a pas de glacier
glacierAbsence <- 1-glacier # pour mettre la favorabilité à 0 quand présence de glacier
Proj_LGM_HG_Salix<-Proj_LGM_Salix*glacierAbsence
plot(Proj_LGM_HG_Salix) 
writeRaster(Proj_LGM_HG_Salix,"6_Cartes_Projections/Salix_LGM_WM_2.tif")

## Binariser les cartes de projection (en sélectionnant le treshold qui maximise le TSS):
# Treshold ici = 459
Salix_bin_current<-BinaryTransformation(Proj_Current, 459)
Salix_bin_2100_MPI26<-BinaryTransformation(Proj_MPI26, 459)
Salix_bin_2100_MPI85<-BinaryTransformation(Proj_MPI85, 459)
Salix_bin_2100_UK26<-BinaryTransformation(Proj_UK26, 459)
Salix_bin_2100_UK85<-BinaryTransformation(Proj_UK85, 459)
Salix_bin_LGM<-BinaryTransformation(Proj_LGM, 459)
plot(Salix_bin_current)
plot(Salix_bin_2100_MPI26)
plot(Salix_bin_2100_MPI85)
plot(Salix_bin_2100_UK26)
plot(Salix_bin_2100_UK85)
plot(Salix_bin_LGM)
writeRaster(Salix_bin_current,"Cartes_Projections/Salix_bin_current.tif")
writeRaster(Salix_bin_2100_MPI26,"Cartes_Projections/Salix_bin_2100_MPI26.tif")
writeRaster(Salix_bin_2100_MPI85,"Cartes_Projections/Salix_bin_2100_MPI85.tif")
writeRaster(Salix_bin_2100_UK26,"Cartes_Projections/Salix_bin_2100_UK26.tif")
writeRaster(Salix_bin_2100_UK85,"Cartes_Projections/Salix_bin_2100_UK85.tif")
writeRaster(Salix_bin_LGM,"Cartes_Projections/Salix_bin_LGM.tif")

## Rangeshifts:

# Load cartes binarisées:
Salix_bin_current<-raster("6_Cartes_Projections/Salix_bin_current.tif")
Salix_bin_2100_MPI26<-raster("6_Cartes_Projections/Salix_bin_2100_MPI26.tif")
Salix_bin_2100_MPI85<-raster("6_Cartes_Projections/Salix_bin_2100_MPI85.tif")
Salix_bin_2100_UK26<-raster("6_Cartes_Projections/Salix_bin_2100_UK26.tif")
Salix_bin_2100_UK85<-raster("6_Cartes_Projections/Salix_bin_2100_UK85.tif")
Salix_bin_LGM<-raster("6_Cartes_Projections/Salix_bin_LGM.tif")

# Pour Proj_LGM: retirer les zones sous glace comme étant favorable:
Glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # charger carte présence de glaciers
# Couper carte glacier sur l'EU:
Mask_EU<-rast("Y_Shapes_Masks/shapeEurope/ForNaWater.tif")
Mask_EU<-terra::project(Mask_EU, "epsg:4326")
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
Glacier<-crop(Glacier,Mask_EU)
Glacier<-resample(Glacier,Salix_bin_LGM,method='ngb')
values(Glacier)[is.na(values(Glacier))]=0 #transformer les na en 0
Salix_LGM_HorsGlacier<-Salix_bin_LGM-Glacier
values(Salix_LGM_HorsGlacier)[values(Salix_LGM_HorsGlacier)<= 0] = 0 #ne garder que les présences (les valeurs=1, les autres =0)
plot(Salix_LGM_HorsGlacier)
writeRaster(Salix_LGM_HorsGlacier,"6_Cartes_Projections/Salix_bin_LGM_2.tif")
Salix_bin_LGM<-raster("6_Cartes_Projections/Salix_bin_LGM_2.tif")

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
Shift_LGM_current<-BIOMOD_RangeSize(Salix_bin_current,Salix_bin_LGM)
Shift_LGM_current<-Shift_LGM_current$Compt.By.Models
Shift_LGM_current<-as.data.frame(Shift_LGM_current)
#plot(Shift_LGM_current$Diff.By.Pixel) # 0 = pixels non-occupé, 1 = pixels occupés slmt ds la nouvelle période
# -1 = pixels occupés aux 2 périodes, -2 = pixels perdus entre current et nouvelle période
diff<-2*Salix_bin_current + Salix_bin_LGM 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Salix_Diff_LGM_2.tif")

# Current-MPI26
Shift_current_2100_MPI26 <- BIOMOD_RangeSize(Salix_bin_current, Salix_bin_2100_MPI26)
Shift_current_2100_MPI26<-Shift_current_2100_MPI26$Compt.By.Models
Shift_current_2100_MPI26<-as.data.frame(Shift_current_2100_MPI26)
# Pour le plot (plus lisible que plot diff.by.pixel):
diff<-2*Salix_bin_current + Salix_bin_2100_MPI26 
plot(diff)  # diff = 1: présence qu'au future (gain, en orange), diff=2: présence qu'au présent (loss, en vert clair), diff=3: présence aux 2 périodes (stable, en vert foncé), 0 = pixels non-occupés
writeRaster(diff, "6_Cartes_Projections/Salix_Diff_MPI26.tif")

# Current-MPI85
Shift_current_2100_MPI85 <- BIOMOD_RangeSize(Salix_bin_current, Salix_bin_2100_MPI85)
Shift_current_2100_MPI85<-Shift_current_2100_MPI85$Compt.By.Models
Shift_current_2100_MPI85<-as.data.frame(Shift_current_2100_MPI85)
diff<-2*Salix_bin_current + Salix_bin_2100_MPI85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Salix_Diff_MPI85.tif")

# Current-UK26
Shift_current_2100_UK26 <- BIOMOD_RangeSize(Salix_bin_current, Salix_bin_2100_UK26)
Shift_current_2100_UK26<-Shift_current_2100_UK26$Compt.By.Models
Shift_current_2100_UK26<-as.data.frame(Shift_current_2100_UK26)
diff<-2*Salix_bin_current + Salix_bin_2100_UK26 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Salix_Diff_UK26.tif")

# Current-UK85
Shift_current_2100_UK85 <- BIOMOD_RangeSize(Salix_bin_current, Salix_bin_2100_UK85)
Shift_current_2100_UK85<-Shift_current_2100_UK85$Compt.By.Models
Shift_current_2100_UK85<-as.data.frame(Shift_current_2100_UK85)
diff<-2*Salix_bin_current + Salix_bin_2100_UK85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Salix_Diff_UK85.tif")

# Fusionner les dataframe:
RangeShifts<-bind_rows(Shift_LGM_current,Shift_current_2100_MPI26,Shift_current_2100_MPI85,
                       Shift_current_2100_UK26,Shift_current_2100_UK85,.id=NULL)
write.csv(RangeShifts,"OutputsCluster/Salix_RangeShifts.csv")


## Shift altitudinal:
altPres <- raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif")
altPast<- raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_dem_-200_V1.0.tif")

# Extraire coordonnées au niveau des présences pour chaque période:
# Current
coordPres <- coordinates(Salix_bin_current)
valPres <- values(Salix_bin_current)
coordPres <- coordPres[!is.na(valPres) & valPres == 1,]
altiPres <- extract(altPres,coordPres)
MedAltPres<-median(altiPres)
# LGM
coordLGM <- coordinates(Salix_bin_LGM)
valLGM <- values(Salix_bin_LGM)
coordLGM <- coordLGM[!is.na(valLGM) & valLGM == 1,]
altiLGM <- extract(altPast,coordLGM)
MedAltLGM<-median(altiLGM)
# MPI26
coordMPI26 <- coordinates(Salix_bin_2100_MPI26)
valMPI26 <- values(Salix_bin_2100_MPI26)
coordMPI26 <- coordMPI26[!is.na(valMPI26) & valMPI26 == 1,]
altiMPI26 <- extract(altPres,coordMPI26)
MedAltMPI26<-median(altiMPI26)
# MPI85
coordMPI85 <- coordinates(Salix_bin_2100_MPI85)
valMPI85 <- values(Salix_bin_2100_MPI85)
coordMPI85 <- coordMPI85[!is.na(valMPI85) & valMPI85 == 1,]
altiMPI85 <- extract(altPres,coordMPI85)
MedAltMPI85<-median(altiMPI85)
# UK26
coordUK26 <- coordinates(Salix_bin_2100_UK26)
valUK26 <- values(Salix_bin_2100_UK26)
coordUK26 <- coordUK26[!is.na(valUK26) & valUK26 == 1,]
altiUK26 <- extract(altPres,coordUK26)
MedAltUK26<-median(altiUK26)
# UK85
coordUK85 <- coordinates(Salix_bin_2100_UK85)
valUK85<- values(Salix_bin_2100_UK85)
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
p <- ggplot(data = pourBoxplot,aes(x=Periode,y=Altitude)) + geom_boxplot() + labs(title="Salix retusa") +theme_classic() 
p


##############################################################################
### Calculer vitesses de migration

# A. Altitudinal

# RangeSchift alt (différences entre médianes; si diff<0: sp est descendue en alt par rapport au présent, si diff>0: sp est montée en alt par rapport au présent):
DiffAlt_Pres_LGM<-MedAltLGM-MedAltPres
DiffAlt_Pres_LGM # -975,1875m
DiffAlt_Pres_MPI26<-MedAltMPI26-MedAltPres
DiffAlt_Pres_MPI26 # 76,6m
DiffAlt_Pres_MPI85<-MedAltMPI85-MedAltPres
DiffAlt_Pres_MPI85 # -978,3m
DiffAlt_Pres_UK26<-MedAltUK26-MedAltPres
DiffAlt_Pres_UK26 # 393,0m
DiffAlt_Pres_UK85<-MedAltUK85-MedAltPres
DiffAlt_Pres_UK85 # 792,0m

# Diviser les rangeshift altitudinal par deltaT:
# Entre LGM et Current (1995-(-21000)=22995 ans) # vitesse min de migration de l'sp
MigrationAlt_Pres_LGM<-DiffAlt_Pres_LGM/22995
MigrationAlt_Pres_LGM # -0,04 m/an soit 4m par décennie

# Entre Current et future (2085-1995= 90 ans) selon scénario # vitesse de migration que l'sp doit avoir pour suivre son habitat optimal:
MigAlt_Pres_MPI26<-DiffAlt_Pres_MPI26/90
MigAlt_Pres_MPI26 # 0,85m/an
MigAlt_Pres_MPI85<-DiffAlt_Pres_MPI85/90
MigAlt_Pres_MPI85 # -10,9m/an
MigAlt_Pres_UK26<-DiffAlt_Pres_UK26/90
MigAlt_Pres_UK26 # 4,4m/an
MigAlt_Pres_UK85<-DiffAlt_Pres_UK85/90
MigAlt_Pres_UK85 # 8,8m/an

# B. latitudinal/longitudinal

# Chargement des coordonnées des centroides (calculé dans QGIS):
Centroides<-read.csv("ShiftLat_Salix_2.csv",header=F,row.names=1, sep=";",dec=",")
# Calcul des distances entre les centroides des différents scenarios:
Distances<-pointDistance(Centroides, lonlat=T,allpairs=T)
