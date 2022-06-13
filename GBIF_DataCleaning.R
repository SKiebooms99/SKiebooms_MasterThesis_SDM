setwd("Z:/projects-unil/GEN4MIG/Solene")

library(rgbif)
library(CoordinateCleaner)
library(readr)

### Cleaning the GBIF data
GBIF_Cirsium <- read_delim("Z:/projects-unil/GEN4MIG/Solene/Occ_Cirsium-spinosissimum/GBIF_Cirsium.csv",escape_double = FALSE, col_names=TRUE, trim_ws=TRUE)
cleaned <- clean_coordinates(x=GBIF_Cirsium,
                             lon = "decimalLongitude",
                             lat = "decimalLatitude",
                             species = "species",
                             countries = NULL,
                             tests = c("capitals", "centroids", "equal", "gbif", "institutions","seas", "zeros"),
                             capitals_rad = 10000,
                             centroids_rad = 1000,
                             centroids_detail = "both",
                             inst_rad = 100,
                             range_rad = 0,
                             zeros_rad = 0.5,
                             capitals_ref = NULL,
                             centroids_ref = NULL,
                             country_ref = NULL,
                             country_refcol = "iso_a3",
                             inst_ref = NULL,
                             range_ref = NULL,
                             seas_ref = NULL,
                             seas_scale = 50,
                             urban_ref = NULL,
                             value = "spatialvalid",
                             verbose = TRUE,
                             report = FALSE)
cleaned <- cleaned[cleaned$.summary,]
View(cleaned)
cleaned <- cleaned[,c("gbifID","species","decimalLongitude","decimalLatitude")]
write.csv(cleaned,"GBIF_Cirsium_cleaned.csv") #sauvegarde du tableau

######
### Merge the GBIF data with InfoFlora
library(dismo)
InfoFlo <- read.csv("InfofloraData3CoordinateSystems.csv")
GBIF_A_alpina<-read.csv("GBIF_Adenostyle-glabra-alpina/A_glabra_alpina_cleaned.csv",row.names = 1)
GBIF_Cirsium<-read.csv("GBIF_Cirsium_cleaned.csv",row.names = 1)
GBIF_P_verticillata<-read.csv("GBIF_Pedicularis-verticillata/Pedicularis-verticillata.csv",row.names = 1)
GBIF_P_alpina<-read.csv("GBIF_Pulsatilla-alpina/Pulsatilla_alpina.csv",row.names = 1)
GBIF_Salix<-read.csv("GBIF_Salix-retusa/Salix_retusa.csv",row.names = 1)
GBIF_Senecio<-read.csv("GBIF_Senecio doronicum/Senecio_doronicum.csv",row.names = 1)
GBIF_Silene<-read.csv("GBIF_Silene-acaulis/Silene_acaulis.csv",row.names = 1)
GBIF_Dryopteris<-read.csv("GBIF_Dryopteris-villarii/Dryopteris_villarii.csv",row.names = 1)

#Pour chaque espèce:
Cirsium<-InfoFlo[InfoFlo$Taxon=="Cirsium spinosissimum (L.) Scop.",c("Taxon","LonWGS84","LatWGS84")] #Selection sp dans infoflora
colnames(Cirsium)<- c("species","decimalLongitude","decimalLatitude")
Cirsium <- rbind.data.frame(Cirsium,GBIF_Cirsium[,-1])
write.csv(Cirsium, "Cirsium_InfoFloGBIF.csv")

#Merge with the GEN4MIG data
GEN4MIG<-read.csv("Gen4MigCleaned.csv")
InfoFloGBIF_Dryopteris<-read.csv("Occ_Dryopteris-villarii/Dryopteris_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_A_alpina<-read.csv("Occ_Adenostyle-glabra-alpina/A_alpina_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_Cirsium<-read.csv("Cirsium_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_Verticillata<-read.csv("Occ_Pedicularis-verticillata/P_verticillata_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_P_alpina<-read.csv("Occ_Pulsatilla-alpina/Pulsatile_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_Salix<-read.csv("Occ_Salix-retusa/Salix_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_Senecio<-read.csv("Occ_Senecio doronicum/Senecio_InfoFloGBIF.csv",row.names = 1)
InfoFloGBIF_Silene<-read.csv("Occ_Silene-acaulis/Silene_InfoFloGBIF.csv",row.names = 1)

#Pour chaque espèce:
Cirsium<-GEN4MIG[GEN4MIG$Species.Name=="Cirsium spinosissimum",c("Species.Name","X","Y")] #Selection sp dans GEN4MIG et des colonnes d'intéret
colnames(Cirsium) <- c("species","decimalLongitude","decimalLatitude")
Cirsium <- rbind.data.frame(Cirsium,InfoFloGBIF_Cirsium)
write.csv(Cirsium,"Cirsium2_full.csv")

#Homogénéiser les noms d'espèce (à faire pour chaque sp):
dryopteris_full<-read.csv("Dryopteris_full.csv",row.names = 1)
dryopteris_full$species = "Dryopteris villarii"

Adenostyle_full<-read.csv("Adenostyle_full.csv",row.names = 1)
Adenostyle_full$species = "Adenostyles alpina"

Cirsium2_full<-read.csv("Cirsium2_full.csv",row.names = 1)
Cirsium2_full$species = "Cirsium spinosissimum"

Pulsatile_full<-read.csv("Pulsatile_full.csv",row.names = 1)
Pulsatile_full$species = "Pulsatila alpina"

Salix_full<-read.csv("Salix_full.csv",row.names = 1)
Salix_full$species = "Salix retusa"

Senecio_full<-read.csv("Senecio_full.csv",row.names = 1)
Senecio_full$species = "Senecio doronicum"

Silene_full<-read.csv("Silene_full.csv",row.names = 1)
Silene_full$species = "Silene acaulis"

Verticillata_full<-read.csv("Verticillata_full.csv",row.names = 1)
Verticillata_full$species = "Pedicularis verticilata"

## Ne garder qu'une coordonnée par pixel (n=1), pour chaque espèce:
xyDryopteris<-dryopteris_full[,c("decimalLongitude","decimalLatitude")] #dataframe à 2 colonnes avec les coordonnées
xyA_alpina<-Adenostyle_full[,c("decimalLongitude","decimalLatitude")]
xyCirsium<-Cirsium2_full[,c("decimalLongitude","decimalLatitude")]
xyPulsatile<-Pulsatile_full[,c("decimalLongitude","decimalLatitude")]
xySalix<-Salix_full[,c("decimalLongitude","decimalLatitude")]
xySenecio<-Senecio_full[,c("decimalLongitude","decimalLatitude")]
xySilene<-Silene_full[,c("decimalLongitude","decimalLatitude")]
xyVerticilata<-Verticillata_full[,c("decimalLongitude","decimalLatitude")]

#r<- bioclim_map #à rechercher dans 'Memoire_SelectVariables':
#Load worldClim data (depuis le nas) de 1981 à 2010:
#bioclim_world <- raster::stack(list.files("Z:/common/50_data/GeoData_other/EnvData/10_World/ChelsaV2.1/1981-2010/bio",pattern = ".tif", full.names = TRUE))
#Couper les cartes pour n'avoir que l'hémisphère Nord (les occurences se trouvant uniquememt dans cet hémispère; diminuera les temps de calcul)
#crop raster with spatial object
#Coordonnées: X=long, Y=lat, extent(xmin,xmax,ymin,ymax)
#Nord<-as(extent(-180,180,0,90),'SpatialPolygons')
#crs(Nord)<-crs(bioclim_world)        
#bioclim_map <- crop(bioclim_world, Nord)

#for(i in 2:20){
#print(i)
#writeRaster(subset(bioclim_map,i), paste0("Z:/projects-unil/GEN4MIG/Solene/",names(bioclim_map)[i],".tif"))
#}

Occ_Verticilata<-gridSample(xy=xyVerticilata,r="CHELSA_bio1_1981.2010_V.2.1.tif", n=1) #A faire pour chaque espèce
write.csv(Occ_Verticilata,"Occ_Verticilata.csv")
