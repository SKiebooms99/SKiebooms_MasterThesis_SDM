### Mod�le Dryopteris EU ###

setwd("Z:/projects-unil/GEN4MIG/Solene")

library(biomod2)
library(raster)
library(foreign)
library(dismo)
library(ggplot2)
library(reshape2)
library(terra)
library(dplyr)

## t�l�charger les donn�es (pour l'Europe) ##

# Esp�ce:
Dryopteris<-read.table("Occurences_data/Occ_Dryopteris_villarii/DryopterisPV.csv",header = TRUE, sep = ",", row.names=1)
Dryopteris<-Dryopteris[,c("decimalLon","decimalLat")] # 444 points
# Supprimer les occurences se trouvant � une altitude aberrante (celles se trouvant en dehors du range altitudinal de l'sp):
altitude<-raster("VariablesEnvi/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif") # doi:10.1038/sdata.2018.40.
AltPoints<-extract(altitude,Dryopteris)
#View(AltPoints)
Dryopteris <- Dryopteris[AltPoints >= 900,] #On retire les points en dessous de 900m -> restent 435

# Charger les variables s�lectionn�es (ClimTopo_sub):
ClimTopo_sub<-stack("ClimTopo_sub.tif")

# G�n�rer des background points (Bp) pour l'Europe (et juste dans les biomes comprenant la majorit� des points d'occurence)
biome <- raster("shapeEurope/GB_EBR_map_BiomesUE.tif")
crs(biome)
Occ_DryopterisProj <- SpatialPointsDataFrame(coords = Dryopteris,data = Dryopteris,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
Occ_DryopterisProj <- spTransform(Occ_DryopterisProj,CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
Occ_DryopterisProj <- as.data.frame(Occ_DryopterisProj)[,-c(1:2)]
biomePoints <- extract(biome, Occ_DryopterisProj)
BiomeToKeep <- as.numeric(names(table(biomePoints))[table(biomePoints)>= (5*nrow(Occ_DryopterisProj)/100)]) # on excluera les biomes contenant moins de 5% des occurences
BiomeToKeep

biomeToSample <- rasterToPolygons(biome,fun=function(x){x==1 | x==7}) #s�lectionner le(s) biome(s) � garder
Points <- extract(ClimTopo_sub, Dryopteris)
Bp <- spsample(biomeToSample, 20000, type="random")
Bp <- spTransform(Bp,CRS("+proj=longlat +datum=WGS84 +no_defs "))
Bp <- as.data.frame(Bp)
BpEnv <- raster::extract(ClimTopo_sub,Bp)
Bp <- na.omit(cbind(Bp,BpEnv)) #retirer les NA dans les lignes
#CheckPointInTheWater <- extract(subset(ClimTopo_sub,1),Bp)
#Bp <- Bp[!is.na(CheckPointInTheWater),]
Bp <- Bp[1:10000,1:2]
write.csv(Bp,"BpCoordinates_Dryopteris_EU2.csv")

#Fusionner les donn�es de pr�sences et les background points:
Bp<-read.csv("BpCoordinates_Dryopteris_EU2.csv", row.names = 1)
colnames(Dryopteris) = c("x","y")
Occ_Dryopteris <- gridSample(Dryopteris,ClimTopo_sub) # Ne garder qu'une occurence max par pixel -> 435
Data_Dryopteris<-rbind(Occ_Dryopteris,Bp)
write.csv(Data_Dryopteris, "Data_Dryopteris_Occ_Bp_ClimTopo_EU2.csv")

###
#Graphique: densit� des points (pr�sence et Bp) pour chaque variable pour l'sp consid�r�e
env <- raster::extract(ClimTopo_sub,Data_Dryopteris)
resp <- c(rep("presence", nrow(Occ_Dryopteris)), rep("background",10000))
tab <- cbind.data.frame(resp=resp,env)
tab$resp=as.factor(tab$resp)
tab <- melt(tab,id="resp")
rm(resp)
png("DensityPlotsVariableEnv_Dryopteris_EU.png",width= 5000,height = 8000,res=300)
p <- ggplot(data=tab,aes(x=value,color=resp,fill=resp)) + geom_density(alpha=0.1) +facet_wrap(~variable,scale= "free")
p
dev.off()
###

## -> Direction cluster pour la calibration et la g�n�ration du mod�le (Code "Calibration.R")
## -> Cluster pour les projections
############################################################################################################################
### Comparer les distributions ###

# Charger les projections (cartes par weighted mean):
Proj_Current<-raster("6_Cartes_Projections/Dryopteris_Current_weighted_mean.tif")
Proj_LGM_Dryopteris<-raster("6_Cartes_Projections/Dryopteris_LGM_weighted_mean.tif")
Proj_MPI26<-raster("6_Cartes_Projections/Dryopteris_MPIESM26_weighted_mean.tif")
Proj_MPI85<-raster("6_Cartes_Projections/Dryopteris_MPIESM85_weighted_mean.tif")
Proj_UK26<-raster("6_Cartes_Projections/Dryopteris_UKESM26_weighted_mean.tif")
Proj_UK85<-raster("6_Cartes_Projections/Dryopteris_UKESM85_weighted_mean.tif")

# Pour Proj_LGM: retirer les zones sous glace comme �tant favorable:
Glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # charger carte pr�sence de glaciers
# Couper carte glacier sur l'EU:
Mask_EU<-rast("Y_Shapes_Masks/shapeEurope/ForNaWater.tif")
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonn�es: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
Glacier<-crop(Glacier,Mask_EU)
glacier<-resample(Glacier,Proj_LGM_Dryopteris,method='ngb')
values(glacier)[is.na(values(glacier))] = 0 #sinon �a vire l� o� il y a des NA et donc quand il n'y a pas de glacier
glacierAbsence <- 1-glacier # pour mettre la favorabilit� � 0 quand pr�sence de glacier
Proj_LGM_HG_Dryopteris<-Proj_LGM_Dryopteris*glacierAbsence
plot(Proj_LGM_HG_Dryopteris) 
writeRaster(Proj_LGM_HG_Dryopteris,"6_Cartes_Projections/Dryopteris_LGM_WM_2.tif")

## Binariser les cartes de projection (en s�lectionnant le treshold qui maximise le TSS):
# Treshold ici = 403
Dryopteris_bin_current<-BinaryTransformation(Proj_Current, 403)
Dryopteris_bin_2100_MPI26<-BinaryTransformation(Proj_MPI26, 403)
Dryopteris_bin_2100_MPI85<-BinaryTransformation(Proj_MPI85, 403)
Dryopteris_bin_2100_UK26<-BinaryTransformation(Proj_UK26, 403)
Dryopteris_bin_2100_UK85<-BinaryTransformation(Proj_UK85, 403)
Dryopteris_bin_LGM<-BinaryTransformation(Proj_LGM, 403)
plot(Dryopteris_bin_current)
plot(Dryopteris_bin_2100_MPI26)
plot(Dryopteris_bin_2100_MPI85)
plot(Dryopteris_bin_2100_UK26)
plot(Dryopteris_bin_2100_UK85)
plot(Dryopteris_bin_LGM)
writeRaster(Dryopteris_bin_current,"Cartes_Projections/Dryopteris_bin_current.tif")
writeRaster(Dryopteris_bin_2100_MPI26,"Cartes_Projections/Dryopteris_bin_2100_MPI26.tif")
writeRaster(Dryopteris_bin_2100_MPI85,"Cartes_Projections/Dryopteris_bin_2100_MPI85.tif")
writeRaster(Dryopteris_bin_2100_UK26,"Cartes_Projections/Dryopteris_bin_2100_UK26.tif")
writeRaster(Dryopteris_bin_2100_UK85,"Cartes_Projections/Dryopteris_bin_2100_UK85.tif")
writeRaster(Dryopteris_bin_LGM,"Cartes_Projections/Dryopteris_bin_LGM_2.tif")

## Rangeshifts:

# load cartes binaris�es:
Dryopteris_bin_current<-raster("6_Cartes_Projections/Dryopteris_bin_current.tif")
Dryopteris_bin_2100_MPI26<-raster("6_Cartes_Projections/Dryopteris_bin_2100_MPI26.tif")
Dryopteris_bin_2100_MPI85<-raster("6_Cartes_Projections/Dryopteris_bin_2100_MPI85.tif")
Dryopteris_bin_2100_UK26<-raster("6_Cartes_Projections/Dryopteris_bin_2100_UK26.tif")
Dryopteris_bin_2100_UK85<-raster("6_Cartes_Projections/Dryopteris_bin_2100_UK85.tif")
Dryopteris_bin_LGM<-raster("6_Cartes_Projections/Dryopteris_bin_LGM.tif")

# Pour Proj_LGM: retirer les zones sous glace comme �tant favorable:
Glacier<-raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_gle_-200_V1.0.tif") # charger carte pr�sence de glaciers
# Couper carte glacier sur l'EU:
Mask_EU<-rast("Y_Shapes_Masks/shapeEurope/ForNaWater.tif")
Mask_EU<-terra::project(Mask_EU, "epsg:4326")
Cut<-as(extent(-25,35,35,75),'SpatialPolygons') #Coordonn�es: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
crs(Cut)<-crs(Mask_EU)        
Mask_EU <- crop(Mask_EU, Cut)
Mask_EU <- raster(Mask_EU)
Glacier<-crop(Glacier,Mask_EU)
Glacier<-resample(Glacier,Dryopteris_bin_LGM,method='ngb')
values(Glacier)[is.na(values(Glacier))]=0 #transformer les na en 0
Dryopteris_LGM_HorsGlacier<-Dryopteris_bin_LGM-Glacier
values(Dryopteris_LGM_HorsGlacier)[values(Dryopteris_LGM_HorsGlacier)<= 0] = 0 #ne garder que les pr�sences (les valeurs=1, les autres =0)
plot(Dryopteris_LGM_HorsGlacier)
writeRaster(Dryopteris_LGM_HorsGlacier,"6_Cartes_Projections/Dryopteris_bin_LGM_2.tif")
Dryopteris_bin_LGM<-raster("6_Cartes_Projections/Dryopteris_bin_LGM_2.tif")

# Current-LGM (par �tapes pour comprendre)
#diff<-2*Pedicularis_bin_current + Pedicularis_bin_LGM 
#plot(diff)  # diff = 1: pr�sence qu'au LGM, diff=2: pr�sence qu'au pr�sent, diff=3: pr�sence aux 2 p�riodes
#valPixel<-freq(diff)
#valPixel # affiche nb de pixels de pr�sence dans chaque cat�gorie (1, 2 ou 3)
#valPixel <- valPixel[-5,] # retirer les NA
# % de diff:
#PixelConserve <- valPixel[valPixel[,1]==3,"count"]/ (valPixel[valPixel[,1]==3,"count"] + valPixel[valPixel[,1]==2,"count"]) # % de l'aire de distribution conserv�
#PixelPresenceOnly <- valPixel[valPixel[,1]==2,"count"]/ (valPixel[valPixel[,1]==3,"count"] + valPixel[valPixel[,1]==2,"count"]) # ou, en compl�mentaire, % de l'aire perdue
#DistPastParRapPres <- (valPixel[valPixel[,1]==1,"count"] +valPixel[valPixel[,1]==3,"count"])/ (valPixel[valPixel[,1]==3,"count"] + valPixel[valPixel[,1]==2,"count"]) # % de l'aire actuelle occup�e au pass�

# Comparer chaque proj au current (utilisation de la fonction BIOMOD_RangeSize, plus rapide):

# Current-LGM
Shift_LGM_current<-BIOMOD_RangeSize(Dryopteris_bin_current,Dryopteris_bin_LGM)
Shift_LGM_current<-Shift_LGM_current$Compt.By.Models
Shift_LGM_current<-as.data.frame(Shift_LGM_current)
#plot(Shift_LGM_current$Diff.By.Pixel) # 0 = pixels non-occup�, 1 = pixels occup�s slmt ds la nouvelle p�riode
# -1 = pixels occup�s aux 2 p�riodes, -2 = pixels perdus entre current et nouvelle p�riode
diff<-2*Dryopteris_bin_current + Dryopteris_bin_LGM 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Dryopteris_Diff_LGM_2.tif")

# Current-MPI26
Shift_current_2100_MPI26 <- BIOMOD_RangeSize(Dryopteris_bin_current, Dryopteris_bin_2100_MPI26)
Shift_current_2100_MPI26<-Shift_current_2100_MPI26$Compt.By.Models
Shift_current_2100_MPI26<-as.data.frame(Shift_current_2100_MPI26)
# Pour le plot (plus lisible que plot diff.by.pixel):
diff<-2*Dryopteris_bin_current + Dryopteris_bin_2100_MPI26 
plot(diff)  # diff = 1: pr�sence qu'au future (gain, en orange), diff=2: pr�sence qu'au pr�sent (loss, en vert clair), diff=3: pr�sence aux 2 p�riodes (stable, en vert fonc�), 0 = pixels non-occup�s
writeRaster(diff, "6_Cartes_Projections/Dryopteris_Diff_MPI26.tif")

# Current-MPI85
Shift_current_2100_MPI85 <- BIOMOD_RangeSize(Dryopteris_bin_current, Dryopteris_bin_2100_MPI85)
Shift_current_2100_MPI85<-Shift_current_2100_MPI85$Compt.By.Models
Shift_current_2100_MPI85<-as.data.frame(Shift_current_2100_MPI85)
diff<-2*Dryopteris_bin_current + Dryopteris_bin_2100_MPI85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Dryopteris_Diff_MPI85.tif")

# Current-UK26
Shift_current_2100_UK26 <- BIOMOD_RangeSize(Dryopteris_bin_current, Dryopteris_bin_2100_UK26)
Shift_current_2100_UK26<-Shift_current_2100_UK26$Compt.By.Models
Shift_current_2100_UK26<-as.data.frame(Shift_current_2100_UK26)
diff<-2*Dryopteris_bin_current + Dryopteris_bin_2100_UK26 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Dryopteris_Diff_UK26.tif")

# Current-UK85
Shift_current_2100_UK85 <- BIOMOD_RangeSize(Dryopteris_bin_current, Dryopteris_bin_2100_UK85)
Shift_current_2100_UK85<-Shift_current_2100_UK85$Compt.By.Models
Shift_current_2100_UK85<-as.data.frame(Shift_current_2100_UK85)
diff<-2*Dryopteris_bin_current + Dryopteris_bin_2100_UK85 
plot(diff)
writeRaster(diff, "6_Cartes_Projections/Dryopteris_Diff_UK85.tif")

# Fusionner les dataframe:
RangeShifts<-bind_rows(Shift_LGM_current,Shift_current_2100_MPI26,Shift_current_2100_MPI85,
                       Shift_current_2100_UK26,Shift_current_2100_UK85,.id=NULL)
write.csv(RangeShifts,"Models_and_Resultats_horsCluster/Dryopteris_RangeShifts.csv")


## Shift altitudinal:
altPres <- raster("3_VariablesEnvi_horsBioclim/VariablesTopo_EU_Amatulli2018/elevation_1KMmn_GMTEDmn.tif")
altPast<- raster("3_VariablesEnvi_horsBioclim/lgm/CHELSA_TraCE21k_dem_-200_V1.0.tif")

# Extraire coordonn�es au niveau des pr�sences pour chaque p�riode:
# Current
coordPres <- coordinates(Dryopteris_bin_current)
valPres <- values(Dryopteris_bin_current)
coordPres <- coordPres[!is.na(valPres) & valPres == 1,]
altiPres <- extract(altPres,coordPres)
MedAltPres<-median(altiPres)
# LGM
coordLGM <- coordinates(Dryopteris_bin_LGM)
valLGM <- values(Dryopteris_bin_LGM)
coordLGM <- coordLGM[!is.na(valLGM) & valLGM == 1,]
altiLGM <- extract(altPast,coordLGM)
MedAltLGM<-median(altiLGM)
# MPI26
coordMPI26 <- coordinates(Dryopteris_bin_2100_MPI26)
valMPI26 <- values(Dryopteris_bin_2100_MPI26)
coordMPI26 <- coordMPI26[!is.na(valMPI26) & valMPI26 == 1,]
altiMPI26 <- extract(altPres,coordMPI26)
MedAltMPI26<-median(altiMPI26)
# MPI85
coordMPI85 <- coordinates(Dryopteris_bin_2100_MPI85)
valMPI85 <- values(Dryopteris_bin_2100_MPI85)
coordMPI85 <- coordMPI85[!is.na(valMPI85) & valMPI85 == 1,]
altiMPI85 <- extract(altPres,coordMPI85)
MedAltMPI85<-median(altiMPI85)
# UK26
coordUK26 <- coordinates(Dryopteris_bin_2100_UK26)
valUK26 <- values(Dryopteris_bin_2100_UK26)
coordUK26 <- coordUK26[!is.na(valUK26) & valUK26 == 1,]
altiUK26 <- extract(altPres,coordUK26)
MedAltUK26<-median(altiUK26)
# UK85
coordUK85 <- coordinates(Dryopteris_bin_2100_UK85)
valUK85<- values(Dryopteris_bin_2100_UK85)
coordUK85 <- coordUK85[!is.na(valUK85) & valUK85 == 1,]
altiUK85 <- extract(altPres,coordUK85)
MedAltUK85<-median(altiUK85)

# Mettre les m�dianes en graphique:
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
p <- ggplot(data = pourBoxplot,aes(x=Periode,y=Altitude)) + geom_boxplot() + labs(title="Dryopteris villarii") +theme_classic() 
p


##############################################################################
### Calculer vitesses de migration

# A. Altitudinal

# RangeSchift alt (diff�rences entre m�dianes; si diff<0: sp est descendue en alt par rapport au pr�sent, si diff>0: sp est mont�e en alt par rapport au pr�sent):
DiffAlt_Pres_LGM<-MedAltLGM-MedAltPres
DiffAlt_Pres_LGM # -471,1875m
DiffAlt_Pres_MPI26<-MedAltMPI26-MedAltPres
DiffAlt_Pres_MPI26 # 119,9m
DiffAlt_Pres_MPI85<-MedAltMPI85-MedAltPres
DiffAlt_Pres_MPI85 # -646,25 m
DiffAlt_Pres_UK26<-MedAltUK26-MedAltPres
DiffAlt_Pres_UK26 # 448,7 m
DiffAlt_Pres_UK85<-MedAltUK85-MedAltPres
DiffAlt_Pres_UK85 # 808,8 m

# Diviser les rangeshift altitudinal par deltaT:
# Entre LGM et Current (1995-(-21000)=22995 ans) # vitesse min de migration de l'sp
MigrationAlt_Pres_LGM<-DiffAlt_Pres_LGM/22995
MigrationAlt_Pres_LGM # -0,02 m/an soit 2m par d�cennie

# Entre Current et future (2085-1995= 90 ans) selon sc�nario # vitesse de migration que l'sp doit avoir pour suivre son habitat optimal:
MigAlt_Pres_MPI26<-DiffAlt_Pres_MPI26/90
MigAlt_Pres_MPI26 # 1,3 m/an
MigAlt_Pres_MPI85<-DiffAlt_Pres_MPI85/90
MigAlt_Pres_MPI85 # -7,2 m/an
MigAlt_Pres_UK26<-DiffAlt_Pres_UK26/90
MigAlt_Pres_UK26 # 4,9 m/an
MigAlt_Pres_UK85<-DiffAlt_Pres_UK85/90
MigAlt_Pres_UK85 # 8,9 m/an

# B. latitudinal/longitudinal

# Chargement des coordonn�es des centroides (calcul�es dans QGIS):
Centroides<-read.csv("ShiftLat_Dryopteris_2.csv",header=F,row.names=1, sep=";",dec=",")
# Calcul des distances entre les centroides des diff�rents scenarios:
Distances<-pointDistance(Centroides, lonlat=T,allpairs=T)
