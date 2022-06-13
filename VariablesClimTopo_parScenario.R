### stack des variables environnementales par scenario

## LGM:
bio2<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_02.tif")
bio6<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_06.tif")
bio12<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_12.tif")
bio15<-raster("CHELSA_bioclim/CHELSA_bioclim_LGM/CHELSA_PMIP_CCSM4_BIO_15.tif")
slope_LGM<-raster("VariablesEnvi/lgm/slope_lgm.tif")

# Couper les cartes sur la zone d'étude:
Mask_EU<-rast("shapeEurope/Mask_EU.tif")
Mask_EU<-terra::project(Mask_EU, "epsg:4326")
#couper le masque à 35° Est (réduire surface = réduire nombre de pixels = réduire temps de calculs):
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)

bio2_LGM<-crop(bio2, Mask_EU)
bio6_LGM<-crop(bio6, Mask_EU)
bio12_LGM<-crop(bio12, Mask_EU)
bio15_LGM<-crop(bio15, Mask_EU)
slope_LGM<-crop(slope_LGM, Mask_EU)
slope_LGM <- resample(slope_LGM,bio15_LGM)

# Stack des variables (+ reprendre la pente):
BioClim_LGM<-stack(bio2_LGM,bio6_LGM,bio12_LGM,bio15_LGM,slope_LGM)

# # Remplacer les valeurs en mer par des NA:
# ForWater <- raster("shapeEurope/ForNAWater.tif")
# ForWater <- crop(ForWater,Mask_EU)
# BioClim_LGM <- mask(BioClim_LGM,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre):
names(BioClim_LGM)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                        "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

writeRaster(BioClim_LGM,"ClimTopo_sub_LGM.tif")

#############################

#Selection de 2 GCM éloigné (https://www.isimip.org/documents/413/ISIMIP3b_bias_adjustment_fact_sheet_Gnsz7CO.pdf) et de deux RCP/SSP

## MPI-ESM, RCP26 et RCP85

bio2_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio2_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio2_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio2_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio6_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio6_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio6_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio12_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio12_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio12_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
bio15_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp126/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1.tif")
bio15_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/MPI-ESM1-2-HR/ssp585/bio/CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1.tif")
slope<-raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

#Couper les cartes sur la zone d'étude:
bio2_RCP26<-crop(bio2_RCP26, Mask_EU)
bio2_RCP85<-crop(bio2_RCP85, Mask_EU)
bio6_RCP26<-crop(bio6_RCP26, Mask_EU)
bio6_RCP85<-crop(bio6_RCP85, Mask_EU)
bio12_RCP26<-crop(bio12_RCP26, Mask_EU)
bio12_RCP85<-crop(bio12_RCP85, Mask_EU)
bio15_RCP26<-crop(bio15_RCP26, Mask_EU)
bio15_RCP85<-crop(bio15_RCP85, Mask_EU)
slope<-crop(slope, Mask_EU)

# Stack des variables (+ reprendre la pente):
BioClim_RCP26<-stack(bio2_RCP26,bio6_RCP26,bio12_RCP26,bio15_RCP26,slope)
BioClim_RCP85<-stack(bio2_RCP85,bio6_RCP85,bio12_RCP85,bio15_RCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26 <- mask(BioClim_RCP26,ForWater)
BioClim_RCP85 <- mask(BioClim_RCP85,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre) et sauvegarder:
names(BioClim_RCP26)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                        "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

writeRaster(BioClim_RCP26,"ClimTopo_sub_MPIESM_RCP26.tif")

names(BioClim_RCP85)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                        "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

writeRaster(BioClim_RCP85,"ClimTopo_sub_MPIESM_RCP85.tif")


## UKSEM: RCP26 et RCP25

bio2_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio2_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio2_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio2_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio6_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio6_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio6_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio6_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio12_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio12_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio12_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio12_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
bio15_RCP26<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp126/bio/CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp126_V.2.1.tif")
bio15_RCP85<-raster("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/2071-2100/UKESM1-0-LL/ssp585/bio/CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1.tif")
slope<-raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/slope_1KMmd_GMTEDmd.tif")

#Couper les cartes sur la zone d'étude:
bio2_RCP26<-crop(bio2_RCP26, Mask_EU)
bio2_RCP85<-crop(bio2_RCP85, Mask_EU)
bio6_RCP26<-crop(bio6_RCP26, Mask_EU)
bio6_RCP85<-crop(bio6_RCP85, Mask_EU)
bio12_RCP26<-crop(bio12_RCP26, Mask_EU)
bio12_RCP85<-crop(bio12_RCP85, Mask_EU)
bio15_RCP26<-crop(bio15_RCP26, Mask_EU)
bio15_RCP85<-crop(bio15_RCP85, Mask_EU)
slope<-crop(slope, Mask_EU)

# Stack des variables (+ reprendre la pente):
BioClim_RCP26<-stack(bio2_RCP26,bio6_RCP26,bio12_RCP26,bio15_RCP26,slope)
BioClim_RCP85<-stack(bio2_RCP85,bio6_RCP85,bio12_RCP85,bio15_RCP85,slope)

# Remplacer les valeurs en mer par des NA:
BioClim_RCP26 <- mask(BioClim_RCP26,ForWater)
BioClim_RCP85 <- mask(BioClim_RCP85,ForWater)

# Donner les mêmes noms aux variables que ceux utilisés pour la construction du modèle (et dans le même ordre) et sauvegarder:
names(BioClim_RCP26)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                        "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

writeRaster(BioClim_RCP26,"ClimTopo_sub_UKESM_RCP26.tif")

names(BioClim_RCP85)<-c("CHELSA_bio2_19812010_V21","CHELSA_bio6_19812010_V21","CHELSA_bio12_19812010_V21",
                        "CHELSA_bio15_19812010_V21","slope_1KMmd_GMTEDmd")

writeRaster(BioClim_RCP85,"ClimTopo_sub_UKESM_RCP85.tif")