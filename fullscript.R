###this is for the full script of our GEB paper: Native and Alien Ungulate Expansion Across North America: A Path to Ecosystem Functionality Restoration
##load the packages
library(raster)
library(terra)
library(data.table)
library(tidyverse)
library(flexsdm)
library(sp)
library(viridis)
library(sdm)
library(mapview)
library(rasterVis)
library(Rmisc)

setwd("/Users/nihao/Library/CloudStorage/Dropbox/Work/NAungulate/forpublic")
###define spaital reference system (crs): 
wgs = '+proj=longlat +datum=WGS84 +no_defs +type=crs'
laea = " +proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs"  #the equal area projection good for North America, which is the major projection in ths study
##the north american map, and raster map 
naboundary = vect('naboundary_no.shp')
naboundary_wgs = project(naboundary, wgs)  ##we also need the WGS84 projection version 
##make the 50km raster 
na_raster50 = rast(naboundary, resolution = 50000)  ##we create 50km raster map 
na_raster50 = rasterize(naboundary, na_raster50, touches = TRUE)
na_raster50 = mask (na_raster50, naboundary)
### the climate variables, from worldclim https://www.worldclim.org/data/bioclim.html .  and we keep the five variables we need 
env_all_raster  = rast('climate_na.tif')

#-------------------------------------------------------------------------------
#now for the current distributions of specific species, there are 25 species in total:
##about the species presence data, which can be downloaded from https://www.inaturalist.org/observations; and here is an instruction: https://help.inaturalist.org/en/support/solutions/articles/151000170342-how-can-i-download-data-from-inaturalist-
na_allherbivores = readRDS('na_allherbivores.rds') #here we compile the 25 species together:
na_allherbivores = na_allherbivores[na_allherbivores$positional_accuracy < 10000, ] ##remove the records with high uncertainty 
na_allherbivores = na_allherbivores[!is.na(na_allherbivores$latitude), ] 
##get the record density, which can be sued to sample background points
recordnum_na = rast('recordnum_NA.tif')
recordnum_na = project(recordnum_na, laea)
# we start with moose: 
moose    = na_allherbivores[na_allherbivores$scientific_name == 'Alces alces', ]
moose.sp = vect(moose, geom = c('longitude','latitude'), crs = wgs )  ##spatialize it 
moose.sp <- project(moose.sp, laea)
dim(moose.sp)
moose.cor  = as.data.frame(crds(moose.sp)); names(moose.cor) = c('longitude','latitude')
na_raster5 = rast('na_raster5.tif') ##we use a 5*5km raster map to sample the background points
### we use the sample_background function from flexsdm 
moose.bg <- sample_background(
  data = moose.cor,
  x='longitude',
  y='latitude',
  n = dim(moose.cor)[1],
  method = "biased",
  rlayer = na_raster5,
  rbias = recordnum_na)
moose.bg.sp = vect(moose.bg, geom = c('longitude','latitude'), crs = laea )
plot(moose.bg.sp)
###extract the environmental variables: 
moose.pre = terra::extract(env_all_raster, moose.sp)
moose.pre$pr_ab = 1
moose.ab  =  terra::extract(env_all_raster, moose.bg.sp)
moose.ab$pr_ab = 0
moose.all =  rbind (moose.pre, moose.ab)
#moose.all$pr_ab = as.factor(moose.all$pr_ab)
###train the model:
moose.all = na.omit(moose.all)
moose_modeldata <- sdmData(pr_ab ~ mat + map  + dryprep + tempseason + prepseason, train=moose.all)
moose_model     <-  sdm(pr_ab~ ., data = moose_modeldata, methods=c('rf')) ##training the model 
moose_modeltest <-  sdm(as.factor(pr_ab)~ ., data=moose_modeldata,  methods=c('rf'), test.percent=20,n=5) #test the model performance
##predict species current distributions 
moose_prediciton = predict(moose_model, newdata=env_all_raster)
plot(moose_prediciton)
## we use this the buffer plot to restrict current distributions 
moose.sp.buffer   = terra::buffer(moose.sp, 50000)
moose_prediciton2 = mask(moose_prediciton, moose.sp.buffer)
values(moose_prediciton2)[is.na(values(moose_prediciton2))] = 0 
moose_prediciton2 = mask(moose_prediciton2, moose_prediciton)
writeRaster(moose_pre_restrict2, './currentdistribution/moose_pre_restrict.tif')
levelplot(moose_prediciton2, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
###we use the same procedure to predict the current distributions of other 25 species, the resulting maps can be found in 'ungulate.allraster.tif'

#------------------------------------------------------------------------------------------------------------------------
###we use a similar way, but global records, to estimate the potential distributions of alien herbivores in North America:
##first we get the global ungulate records from iNaturalized and clean as ungulates_unique:
ungulates_unique = readRDS('ungulates_all.rds')
ungulates_unique$scientific_name = word(ungulates_unique$scientific_name, 1,2, sep=" ") #get the binary name
ungulates_unique = ungulates_unique[!is.na(ungulates_unique$latitude), ] #remove the missing values 
#change to spatial data: 
ungulates_loc = vect    (ungulates_unique, geom = c('longitude','latitude'), crs = wgs )
ungulates_loc = project (ungulates_loc, laea)
mapview(ungulates_loc) ##this is the global records
recordnum_global = rast('recordnum_global.tif') ##we also calculate the global record density
globalenv_all_raster = rast('globalenv_all_raster.tif')
globalenv_all_raster = project(globalenv_all_raster, laea)
######we start from Axis axis:
Axis.axis    = ungulates_loc[ungulates_loc$scientific_name=="Axis axis", ]
axis.cor     = as.data.frame(crds(Axis.axis)); names(axis.cor) = c('longitude','latitude')
###because Axis native to eurasia, so we only sample the background points from Euroasia: 
eurasiamap_raster5 = rast('eurasiamap_raster5.tif')
axis.global.bg1 <- sample_background(
  data = axis.cor,
  x ='longitude',
  y ='latitude',
  n = dim(axis.cor)[1],
  method = "biased",
  rlayer = eurasiamap_raster5,
  rbias =  recordnum_global )
axis.global.bgsp    = vect(axis.global.bg1, geom = c('longitude','latitude'), crs = laea ) 
axisdeer.global.pre =  terra::extract(globalenv_all_raster, Axis.axis)
axisdeer.global.pre$pr_ab = 1
axisdeer.global.ab  =  terra::extract(globalenv_all_raster, axis.global.bgsp)
axisdeer.global.ab$pr_ab = 0
#combine it:
axisdeer.global.all =  rbind(axisdeer.global.pre, axisdeer.global.ab)
axisdeer.global.all = na.omit(axisdeer.global.all)

axisdeer_globalmodeldata <- sdmData(as.factor(pr_ab)~  mat + map  + dryprep + tempseason + prepseason, 
                                     train = axisdeer.global.all)
axisdeer_globalmodel      <-  sdm(as.factor(pr_ab)~ . , data=axisdeer_globalmodeldata, methods=c('rf'))
axisdeer_globalmodel.test <-  sdm(as.factor(pr_ab)~ .,  data=axisdeer_globalmodeldata, methods=c('rf'), test.percent=20,n=5)
##make the prediction: 
axis.potential.prediction = predict(axisdeer_globalmodel, newdata=env_all_raster)
levelplot(axis.potential.prediction, margin = FALSE, scales=list(draw=FALSE), par.settings = RdBuTheme)
writeRaster(axis.potential.prediction, 'axis.potential.prediction.tif', overwrite=TRUE)
###we use the same procedure for the potential distributions of remaining alien species 


######--------------------------------------------------------------------------------------------------------
######then we give analyses for richness pattern and functional composition: start from current distributions 
##first we transform the predicted occurrence probability into binary presence/absence
ungulate.allraster = rast('ungulate.allraster.tif') #this is current distributions
names(ungulate.allraster) = c('Odocoileus virginianus', 'Odocoileus hemionus', 'Rangifer tarandus', 'Cervus canadensis', 'Alces alces','Oreamnos americanus',
                              'Bison bison', 'Ovis canadensis' ,'Ovis dalli' , 'Ovibos moschatus',  'Pecari tajacu',   'Antilocapra americana',
                              'Axis axis' , 'Dama dama','Oryx gazella' ,'Cervus nippon','Ovis orientalis aries', 'Bos taurus',
                              'Capra aegagrus hircus', 'Boselaphus tragocamelus', 'Antilope cervicapra', 'Ammotragus lervia' , 'Equus ferus caballus', 'Equus africanus asinus', 'Sus scrofa')
model_current      = readRDS('model_current.rds') # the sdm model for each species 
ungulate.allraster2= ungulate.allraster
for (i in 1:25) {ungulate.allraster2[[i]] = pa(ungulate.allraster[[i]], model_current[[i]], opt = 2 )} ##pa is from 'sdm' package

###map the richness
native.raster = ungulate.allraster2[[1:12]]
alien.raster = ungulate.allraster2[[13:25]]
###summing  up the richness 
ungulate.richness = sum(ungulate.allraster2) 
plot(ungulate.richness)
mean(values(ungulate.richness), na.rm=TRUE) ##1.273421
sd(values(ungulate.richness), na.rm=TRUE) #1.7684
native.richness = sum(native.raster) 
mean(values(native.richness), na.rm=TRUE) #0.795
sd(values(native.richness), na.rm=TRUE) ##1.04
alien.richness = sum(alien.raster)
mean(values(alien.richness), na.rm=TRUE) ##0.48
sd(values(alien.richness), na.rm=TRUE)   ##1.21
#map the richness pattern:
all.richness = c(ungulate.richness, native.richness, alien.richness)
names(all.richness) = c('All', 'Native', 'Alien')
levelplot(all.richness, margin = FALSE, 
          par.settings = list(axis.line = list(col = "transparent"), 
                              strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c("blue", "white", "red"))(100),
          scales=list(draw=FALSE))

###show the ecoregion differences: 
ecoregion  =  vect('./Ecoregions2017/Ecoregions2017.shp')
ecoregion  =  project(ecoregion, laea)
ecoregioncurrent =  rasterize (ecoregion, all.richness, "BIOME_NAME")
allrichness_eco  =  c(all.richness, ecoregioncurrent )
plot(allrichness_eco)
allrichness_ecodata = values(allrichness_eco, dataframe=TRUE)
allrichness_ecodata = allrichness_ecodata[!is.na(allrichness_ecodata$All),]
table(allrichness_ecodata$BIOME_NAME)
###we only select part of ecoregions with sufficient grids
select_reg = c('Boreal Forests/Taiga', 'Deserts & Xeric Shrublands', 
               'Temperate Broadleaf & Mixed Forests', 
               'Temperate Conifer Forests', 'Temperate Grasslands, Savannas & Shrublands', 'Tundra')
allrichness_eco = allrichness_ecodata[allrichness_ecodata$BIOME_NAME %in% select_reg, ]
names(allrichness_eco) =  c('All', 'Native', 'Alien','Ecoregion')
#####rename the data
allrichness_eco2 = gather(allrichness_eco, 'group','richness',1:3)
allrichness_eco2 $Ecoregion = as.character(allrichness_eco2$Ecoregion)
allrichness_eco2$Ecoregion[allrichness_eco2$Ecoregion=='Temperate Broadleaf & Mixed Forests'] = 'Broadleaf forests'
allrichness_eco2$Ecoregion[allrichness_eco2$Ecoregion=='Boreal Forests/Taiga'] = 'Boreal forests'
allrichness_eco2$Ecoregion[allrichness_eco2$Ecoregion=='Temperate Grasslands, Savannas & Shrublands'] = 'Grasslands'
allrichness_eco2$Ecoregion[allrichness_eco2$Ecoregion=='Deserts & Xeric Shrublands'] = 'Deserts'
allrichness_eco2$Ecoregion[allrichness_eco2$Ecoregion=='Temperate Conifer Forests'] = 'Temperate conifers'
allrichness_eco2$Ecoregion  = factor(allrichness_eco2$Ecoregion , levels = c('Deserts','Temperate conifers','Grasslands','Broadleaf forests',  'Boreal forests', 'Tundra'))
allrichness_eco2$group      = factor(allrichness_eco2$group , levels = c('Alien', 'Native','All'))
allrichness_eco2    = allrichness_eco2[!is.na(allrichness_eco2$richness),]
allrichness_eco_sum = summarySE(allrichness_eco2, measurevar="richness", groupvars=c("Ecoregion","group"))  ##summarySE calculate the mean and SE, we provide the data lower site;
cbPalette <- c("#cb181d",   "#f7f7f7", "#525252") #the colour plan
ggplot(allrichness_eco_sum, aes(x=Ecoregion, y=richness, fill=group))  + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=richness-se, ymax=richness+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  xlab("") + theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylab("Herbivore richness") +  theme_bw() + 
  theme(axis.text.x = element_text(color="black", 
                                   size=8)) +
  scale_fill_manual(values=cbPalette)

##----------------------------------------------------------------------------------------------------
##now for the trait composition: 
##change the raster to matrix： 
##change the raster data to matrix: 
ungulate.matrix = as.data.frame(values(ungulate.allraster))
rownames(ungulate.matrix) = 1:21222
native.matrix = as.data.frame(values(native.raster))
rownames(native.matrix) = 1:21222
alien.matrix = as.data.frame(values(alien.raster))
rownames(alien.matrix) = 1:21222
#get the HerbiTraits data 
HerbiTraits_1.2 <- read.csv("./HerbiTraits_1.2/HerbiTraits_1.2.csv")
##we only choose the functional traits useful for us:
herbtrait = HerbiTraits_1.2 [ , c(1, 6, 7, 8, 13, 14, 16, 22, 23, 28,33) ]
our_sp = data.frame(species = names(ungulate.matrix) )
#transform the data and fix the name
herbtrait_use = herbtrait[, -c(2, 3)]
rownames(herbtrait_use) = herbtrait_use$Binomial
herbtrait_use$Diet.Graminoids = as.numeric(herbtrait_use$Diet.Graminoids + 1, ordered = TRUE)
herbtrait_use$Diet.Browse.Fruit = as.numeric(herbtrait_use$Diet.Browse.Fruit + 1, ordered = TRUE)
herbtrait_use$Mass.g = as.numeric(herbtrait_use$Mass.g) 
herbtrait_use$Binomial[herbtrait_use$Binomial=='Bos primigenius taurus'] = 'Bos taurus' ##fix the name 

#now we calculate the community-level functional features:  
ungulate.matrix.all = ungulate.matrix
herbtrait.current    = herbtrait_use[rownames(herbtrait_use)%in% names(ungulate.matrix.all), ]
ungulate.matrix.all  = ungulate.matrix.all [,rownames(herbtrait.current)]
###this is sum of body mass: 
site.bodymass   = as.matrix(ungulate.matrix.all) %*% herbtrait.current$Mass.g
bodymass.raster = ungulate.richness
values(bodymass.raster)  = site.bodymass [, 1]
plot(bodymass.raster)
bodymass.raster = bodymass.raster/1000
levelplot(bodymass.raster, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
###sum of grazing 
site.grazing.total           =  as.matrix(ungulate.matrix.all) %*% herbtrait.current$Diet.Graminoids
grazing.raster.total          =  ungulate.richness;   values(grazing.raster.total)  =  site.grazing.total
levelplot(grazing.raster.total, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
###sum of browsing 
site.browsing.total            =  as.matrix(ungulate.matrix.all) %*% herbtrait.current$Diet.Browse.Fruit
browsing.raster.total          =  ungulate.richness;   values(browsing.raster.total)  =  site.browsing.total
levelplot(browsing.raster.total, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

##now for the native and invasive species separately: 
#for native species 
herbtrait.currentnative    = herbtrait_use[rownames(herbtrait_use)%in% names(native.matrix), ]
nativesite.bodymass   = as.matrix(native.matrix) %*% herbtrait.currentnative$Mass.g
nativebodymass.raster = ungulate.richness
values(nativebodymass.raster)  = nativesite.bodymass [, 1]
nativebodymass.raster = nativebodymass.raster/1000
levelplot(nativebodymass.raster, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
site.grazing.native            =  as.matrix(native.matrix)%*% herbtrait.currentnative$Diet.Graminoids  
grazing.raster.native          =  ungulate.richness;   values(grazing.raster.native)  =  site.grazing.native
levelplot(grazing.raster.native, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
site.browsing.native            =  as.matrix(native.matrix) %*% herbtrait.currentnative$Diet.Browse.Fruit 
browsing.raster.native          =  ungulate.richness;   values(browsing.raster.native)  =  site.browsing.native
levelplot(browsing.raster.native, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
#for alien species 
herbtrait.currentalien  =  herbtrait_use[rownames(herbtrait_use)%in% names(alien.matrix), ]
aliensite.bodymass   = as.matrix(alien.matrix) %*% herbtrait.currentalien$Mass.g
alienbodymass.raster = ungulate.richness
values(alienbodymass.raster)  = aliensite.bodymass [, 1]
alienbodymass.raster = alienbodymass.raster/1000
levelplot(alienbodymass.raster, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
site.grazing.alien            =  as.matrix(alien.matrix)%*% herbtrait.currentalien$Diet.Graminoids  
grazing.raster.alien          =  ungulate.richness;   values(grazing.raster.alien)  =  site.grazing.alien
levelplot(grazing.raster.alien, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
site.browsing.alien            =  as.matrix(alien.matrix) %*% herbtrait.currentalien$Diet.Browse.Fruit 
browsing.raster.alien          =  ungulate.richness;   values(browsing.raster.alien)  =  site.browsing.alien
levelplot(browsing.raster.alien, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
####draw the functional composition: 
bodymass.raster = nativebodymass.raster+alienbodymass.raster
bodymass.all = c(bodymass.raster, nativebodymass.raster, alienbodymass.raster)
names (bodymass.all) = c('All', 'Native', 'Alien')
levelplot(bodymass.all, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
grazing.all = stack(raster(grazing.raster.total), raster(grazing.raster.native), raster(grazing.raster.alien))
names (grazing.all) = c('All', 'Native', 'Alien')
levelplot(grazing.all, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))
browsing.all = stack(raster(browsing.raster.total), raster(browsing.raster.native), raster(browsing.raster.alien))
names (browsing.all) = c('All', 'Native', 'Alien')
levelplot(browsing.all, margin = FALSE, par.settings=rasterTheme(viridis_pal(option = "D")(255)), scales=list(draw=FALSE))



#---------------------------------------------------------------------------------------------------------------
#now for the alien species potential richness and composition: 
all_ungupotential = rast('all_ungupotential.tif')
names(all_ungupotential) = c('Axis axis',  'Rucervus duvaucelii', 'Ammotragus lervia', 'Antilope cervicapra', 'Bos taurus', 
                             'Dama dama', 'Equus africanus asinus', 'Capra aegagrus hircus','Hemitragus jemlahicus', 'Equus ferus caballus', 
                             'Aepyceros melampus', 'Kobus leche', 'Boselaphus tragocamelus','Tragelaphus angasii', 
                             'Oryx dammah', 'Oryx gazella', 'Sus scrofa', 'Rusa unicolor', 'Ovis orientalis aries', 
                             'Cervus nippon', 'Antidorcas marsupialis', 'Phacochoerus africanus', 'Kobus ellipsiprymnus' )
##change to binary occurrence: 
potential_model = readRDS('./potential_model.rds')
all_ungupotential2 =  all_ungupotential
for (i in 1:23) {all_ungupotential2[[i]] = pa(all_ungupotential2[[i]],  potential_model[[i]], opt = 2 )}
#first, give the richness: 
potentialrichness.raster = sum(all_ungupotential2)
mean(values(potentialrichness.raster),na.rm=TRUE)  #2.19
sd(values(potentialrichness.raster),na.rm=TRUE)  ##2.68
levelplot(potentialrichness.raster, 
          margin = list(FUN = 'mean'),
          par.settings = list(
            strip.background = list(col = 'transparent'), 
            strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c("blue", "white", "red"))(100),
          scales=list(draw=FALSE))
##get the traits:
ungulatepotential.matrix        = as.data.frame (values(all_ungupotential2))
herbtrait.potential = herbtrait_use[herbtrait_use$Binomial %in% names(ungulatepotential.matrix), ]
ungulatepotential.matrix = ungulatepotential.matrix[, herbtrait.potential$Binomial]
##now mean body mass, change community to binary data: 
potentialsite.bodymasstotal     =  (as.matrix(ungulatepotential.matrix) %*% herbtrait.potential$Mass.g )
potentialbodymasstotal.raster   =  potentialrichness.raster ; values(potentialbodymasstotal.raster)  =  potentialsite.bodymasstotal
potentialbodymasstotal.raster = potentialbodymasstotal.raster/1000
levelplot(potentialbodymasstotal.raster, 
          margin = list(FUN = 'mean'),
          par.settings = list(
            strip.background = list(col = 'transparent'), 
            strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c("blue", "white", "red"))(100),
          scales=list(draw=FALSE))
##the grazing : 
site.grazingtotal  =  as.matrix(ungulatepotential.matrix) %*% herbtrait.potential$Diet.Graminoids 
potentialgrazingtotal.raster          =  potentialrichness.raster ; values(potentialgrazingtotal.raster)  =  site.grazingtotal
levelplot(potentialgrazingtotal.raster, 
          margin = list(FUN = 'mean'),
          par.settings = list(
            strip.background = list(col = 'transparent'), 
            strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c("blue", "white", "red"))(100),
          scales=list(draw=FALSE))
##the browsing behavior: 
site.browsingtotal       =  as.matrix(ungulatepotential.matrix) %*% herbtrait.potential$Diet.Browse.Fruit 
potentialbrowsingtotal.raster          =  potentialrichness.raster ; values(potentialbrowsingtotal.raster)  =  site.browsingtotal
levelplot(potentialbrowsingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
levelplot(potentialbrowsingtotal.raster, 
          margin = list(FUN = 'mean'),
          par.settings = list(
            strip.background = list(col = 'transparent'), 
            strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c("blue", "white", "red"))(100),
          scales=list(draw=FALSE))
###give a ecoregion-level analsis: 
potentialrichness_eco = c   (potentialrichness.raster,    potentialbodymasstotal.raster,
                             potentialgrazingtotal.raster,potentialbrowsingtotal.raster,
                             ecoregioncurrent)
plot(potentialrichness_eco)
names(potentialrichness_eco) = c('richness','bodymass','grazing','browsing','BIOME_NAME')
potentialrichness_ecodata = values(potentialrichness_eco, dataframe=TRUE)
potentialrichness_ecodata = potentialrichness_ecodata[!is.na(potentialrichness_ecodata$richness),]
select_reg = c('Boreal Forests/Taiga', 'Deserts & Xeric Shrublands', 
               'Temperate Broadleaf & Mixed Forests', 
               'Temperate Conifer Forests', 'Temperate Grasslands, Savannas & Shrublands', 'Tundra')
potentialrichness_eco2 = potentialrichness_ecodata[potentialrichness_ecodata$BIOME_NAME %in% select_reg, ]
potentialrichness_eco2$BIOME_NAME = as.character(potentialrichness_eco2$BIOME_NAME)
potentialrichness_eco2$BIOME_NAME[potentialrichness_eco2$BIOME_NAME ==  'Temperate Broadleaf & Mixed Forests'] = 'Broadleaf Forests'
potentialrichness_eco2$BIOME_NAME[potentialrichness_eco2$BIOME_NAME ==  'Boreal Forests/Taiga'] = 'Boreal Forests'
potentialrichness_eco2$BIOME_NAME[potentialrichness_eco2$BIOME_NAME ==  'Temperate Grasslands, Savannas & Shrublands'] = 'Grasslands'
potentialrichness_eco2$BIOME_NAME[potentialrichness_eco2$BIOME_NAME ==  'Deserts & Xeric Shrublands'] = 'Deserts'
potentialrichness_eco2$BIOME_NAME[potentialrichness_eco2$BIOME_NAME ==  'Boreal Forests/Taiga'] = 'Boreal Forests'
potentialrichness_eco2$BIOME_NAME[potentialrichness_eco2$BIOME_NAME ==  'Temperate Conifer Forests'] = 'Temperate Conifers'
potentialrichness_eco2$BIOME_NAME  = factor(potentialrichness_eco2$BIOME_NAME, 
                                            levels = c('Deserts','Grasslands',
                                                       'Broadleaf Forests', 'Temperate Conifers', 
                                                       'Boreal Forests', 'Tundra'))
names(potentialrichness_eco2) = c('richness','biomass','grazing','browsing','ecoregion')
potentialrichness_eco2_richness = summarySE(potentialrichness_eco2, measurevar= "richness", groupvars=c("ecoregion"))
potentialrichness_eco2_biomass = summarySE(potentialrichness_eco2, measurevar = "biomass", groupvars=c("ecoregion"))
potentialrichness_eco2_grazing = summarySE(potentialrichness_eco2, measurevar = "grazing", groupvars=c("ecoregion"))
potentialrichness_eco2_browsing = summarySE(potentialrichness_eco2, measurevar= "browsing", groupvars=c("ecoregion"))

#------------------------------------------------------------------------------------------------
#about the present-natural species before extinctions: 
allherbivore = rast("allherbivore_natual.tif") ###download the data from PHYLACINE 
names(allherbivore) = c( "Alces_alces","Antilocapra_americana","Bison_bison","Bootherium_bombifrons",
"Camelops_hesternus","Capromeryx_minor","Castoroides_ohioensis","Cervalces_scotti", "Cervus_canadensis","Cuvieronius_hyodon","Dasypus_bellus","Equus_ferus",
"Equus_francisci","Eremotherium_laurillardi",  "Euceratherium_collinum",    "Glyptotherium_cylindricum", "Glyptotherium_floridanum", "Hemiauchenia_macrocephala", "Holmesina_septentrionalis", "Mammut_americanum",        
"Mammuthus_columbi",         "Mammuthus_primigenius",    "Megalonyx_jeffersonii",     "Mixotoxodon_larensis",   
"Muknalia_minima" ,          "Mylohyus_nasutus" ,         "Navahoceros_fricki"   ,     "Neochoerus_aesopi"   ,     
"Nothrotheriops_shastensis" ,"Nothrotherium_maquinense" , "Odocoileus_hemionus"  ,     "Odocoileus_virginianus",  
"Oreamnos_americanus"   ,    "Oreamnos_harringtoni"  ,    "Ovibos_moschatus"   ,       "Ovis_canadensis"  ,        
"Ovis_dalli"     ,           "Palaeolama_mirifica"  ,     "Paramylodon_harlani"  ,     "Pecari_tajacu",            
"Platygonus_compressus"   ,  "Rangifer_tarandus"  ,       "Saiga_tatarica"     ,       "Sangamona_fugitiva" ,     
"Stockoceros_conklingi"   ,  "Tapirus_merriami"    ,      "Tapirus_veroensis"     ,    "Tetrameryx_shuleri"  ,    
"Tremarctos_floridanus")
#give some clean:
allherbivore2   = project(allherbivore, laea)
allherbivore2   = mask(allherbivore2, naboundary)
allherbivore2   = crop(allherbivore2, ext(ungulate.richness))
ext(allherbivore2) = ext(ungulate.richness)
#change to matrix: 
naturalungulate.matrix = as.data.frame(values(allherbivore2))
names(naturalungulate.matrix) = stri_replace_all_fixed(names(allherbivore), '_', ' ')
#first, give the richness: 
naturalrichness = sum(allherbivore2)
mean(values(naturalrichness), na.rm=TRUE) ##13.98
sd(values(naturalrichness), na.rm=TRUE)   ##9.41
levelplot(naturalrichness, margin = FALSE, 
          par.settings = list(axis.line = list(col = "transparent"), 
                              strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c("blue", "white", "red"))(100),
          scales=list(draw=FALSE))

##get the traits:
herbtrait.natural    = herbtrait_use[rownames(herbtrait_use) %in% names(naturalungulate.matrix), ]
naturalungulate.matrix      =    naturalungulate.matrix[, herbtrait.natural$Binomial]
naturalsite.bodymasstotal           =  (as.matrix(naturalungulate.matrix) %*% herbtrait.natural$Mass.g )
naturalbodymasstotal.raster         =  naturalrichness ; values(naturalbodymasstotal.raster)  =  naturalsite.bodymasstotal
naturalbodymasstotal.raster = naturalbodymasstotal.raster/1000
levelplot(naturalbodymasstotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the grazing behavior: 
site.grazingtotal       =  as.matrix(naturalungulate.matrix) %*% herbtrait.natural$Diet.Graminoids 
naturalgrazingtotal.raster          =  naturalrichness ; values(naturalgrazingtotal.raster)  =  site.grazingtotal
levelplot(naturalgrazingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the browsing behavior: 
site.browsingtotal       =  as.matrix(naturalungulate.matrix) %*% herbtrait.natural$Diet.Browse.Fruit 
naturalbrowsingtotal.raster          =  naturalrichness ; values(naturalbrowsingtotal.raster)  =  site.browsingtotal
levelplot(naturalbrowsingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
natural.all = c(naturalrichness, naturalbodymasstotal.raster, naturalgrazingtotal.raster, naturalbrowsingtotal.raster )
names(natural.all) = c('Richness', 'Total body mass', 'Grazing intensity', 'Browzing intensity')
plot(natural.all)


##--------------------------------------------------------------------------------------------------------
##now we compare the present-natural baseline with current native, current all richness, current all + alien potential
##resample current distributions to 100 km scale 
native.sp2 = resample(native.raster,  naturalrichness,  method = 'max')
native.richness2 = sum(native.sp2)
plot(native.richness2)
##resample alien species: 
alien.sp2 = resample(alien.raster, naturalrichness, method = 'max')
alien.richness2 = sum(alien.sp2)
plot(alien.richness2)
###resample the alien potential distribution 
all_ungupotential3 = resample(all_ungupotential2, naturalrichness, method = 'max')
alien.potentialrichness = sum (all_ungupotential3)
plot(alien.potentialrichness)
##get the native potential distributions: 
native_potenital  =  c(allherbivore2$Odocoileus_virginianus, allherbivore2$Odocoileus_hemionus, 
                       allherbivore2$Rangifer_tarandus, allherbivore2$Cervus_canadensis, allherbivore2$Alces_alces,
                       allherbivore2$Oreamnos_americanus, allherbivore2$Bison_bison, allherbivore2$Ovis_canadensis,
                       allherbivore2$Ovis_dalli, allherbivore2$Ovibos_moschatus, allherbivore2$Pecari_tajacu, 
                       allherbivore2$Antilocapra_americana)
native_potenitalrichness = sum(native_potenital)
plot(native_potenitalrichness)
all_potenital = c(naturalrichness, native_potenitalrichness, alien.potentialrichness)
names(all_potenital) = c('Natural','CurrentNative','AlienPotential')
levelplot(all_potenital, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

#####----------------------------------------------------------------------------------
##compare with the natural baseline:
current_allrichness2  = native.richness2 + alien.richness2
potenital_allrichness = native_potenitalrichness + alien.potentialrichness
richness.all = c(naturalrichness,          native.richness2, 
                 current_allrichness2,     native_potenitalrichness, 
                 potenital_allrichness)
names(richness.all) = c('Natural_richness', 'Current_native_richness',
                        'Current_native_alien_richness',
                        'Native_potential_richness', 
                        'Native_alien_potential_richness')
levelplot(richness.all, margin = FALSE, scales=list(draw=FALSE))
###calculate the differences: 
diff_currentnative  = naturalrichness - native.richness2   ##current native deficit
diff_currentnative2 = diff_currentnative / naturalrichness
diff_currentnative2[is.infinite(diff_currentnative2)]=NA
mean(values(diff_currentnative2), na.rm=TRUE) ##0.842
diff_currentnative[diff_currentnative<0] = NA  #we exclude the grids with negative value 
levelplot(diff_currentnative, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

current_allrichness2 = native.richness2 + alien.richness2  ## current all species 
diff_currentall = naturalrichness - current_allrichness2
diff_currentall = mask (diff_currentall, diff_currentnative)
levelplot(diff_currentall, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

nativepotential = native_potenitalrichness + alien.richness2
diff_potentialnative = naturalrichness - nativepotential
diff_potentialnative = mask (diff_potentialnative, diff_currentnative)
levelplot(diff_potentialnative, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

potenital_allrichness = native_potenitalrichness + alien.potentialrichness
diff_potentialall = naturalrichness - potenital_allrichness
diff_potentialall = mask(diff_potentialall, diff_currentnative)
levelplot(native_potenitalrichness, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
levelplot(diff_potentialall, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
deficit_all = c(diff_currentall, diff_potentialnative, diff_potentialall)
names(deficit_all) = c('Current_Alien', 'Potential_Native', 'Potential_Native_Alien')
levelplot(deficit_all, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##get the values of deficit: 
mean(values(diff_currentnative),na.rm=TRUE)  #13.09
mean(values(deficit_all$Current_Alien),na.rm=TRUE) #12.23
mean(values(deficit_all$Potential_Native),na.rm=TRUE)  #10.17
mean(values(deficit_all$Potential_Native_Alien),na.rm=TRUE)  #8.35
##get the proportion of changes:
deficit_percent  =  (c(values(diff_currentnative)) - values(deficit_all)) / c(values(diff_currentnative))
deficit_percent [is.infinite(deficit_percent)] <- NA
deficit_percent [deficit_percent < 0] = 0
deficit_percent [deficit_percent > 1 ] = 1
deficit_percent_raster          =  deficit_all 
values(deficit_percent_raster)  =  deficit_percent*100
plot(deficit_percent_raster)
names(deficit_percent_raster) = c('Current_alien','Potential_Native','Potential_Native_Alien')
levelplot(deficit_percent_raster*100, margin = FALSE,  par.strip.text = list(cex = 0),
          par.settings = list(axis.line = list(col = "transparent"), 
                              strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c('blue', "#fed976", "red"))(100),
          scales=list(draw=FALSE))
##save in 825*275
mean(deficit_percent[,1], na.rm=TRUE) ##0.047
mean(deficit_percent[,2], na.rm=TRUE) ##0.3026
mean(deficit_percent[,3], na.rm=TRUE)  ##0.4206
####we also compare between ecoregions: 
ecoregion2  =  rasterize (ecoregion, naturalrichness, "BIOME_NAME")
deficit_eco =  c(deficit_percent_raster, ecoregion2 )
plot(deficit_eco)
deficit_ecodata = values(deficit_eco, dataframe=TRUE)
deficit_ecodata = deficit_ecodata[!is.na(deficit_ecodata$Current),]
table(deficit_ecodata$BIOME_NAME)
###we may only analyse some biomes, like:
select_reg = c('Boreal Forests/Taiga', 'Deserts & Xeric Shrublands', 
               'Temperate Broadleaf & Mixed Forests', 
               'Temperate Conifer Forests', 'Temperate Grasslands, Savannas & Shrublands', 'Tundra')
deficit_eco = deficit_ecodata[deficit_ecodata$BIOME_NAME %in% select_reg, ]
deficit_eco2 = gather(deficit_eco, 'DeficitRestoration','Percent',1:3)
ggplot(deficit_eco2, aes(x=BIOME_NAME, y=Percent, fill=DeficitRestoration)) + geom_boxplot() 
class(deficit_eco2$BIOME_NAME)
deficit_eco2$BIOME_NAME = as.character(deficit_eco2$BIOME_NAME)
deficit_eco2$BIOME_NAME [deficit_eco2$BIOME_NAME =='Temperate Broadleaf & Mixed Forests'] = 'Temperate Broadleaf'
deficit_eco2$BIOME_NAME [deficit_eco2$BIOME_NAME =='Boreal Forests/Taiga'] = 'Boreal Forests'
deficit_eco2$BIOME_NAME [deficit_eco2$BIOME_NAME =='Temperate Grasslands, Savannas & Shrublands'] = 'Grasslands'
deficit_eco2$BIOME_NAME [deficit_eco2$BIOME_NAME =='Deserts & Xeric Shrublands'] = 'Deserts'
deficit_eco2$BIOME_NAME [deficit_eco2$BIOME_NAME =='Temperate Conifer Forests'] = 'Temperate Conifers'
deficit_eco2$BIOME_NAME = factor(deficit_eco2$BIOME_NAME, 
                                 levels = c('Deserts','Grasslands',
                                            'Temperate Broadleaf', 'Temperate Conifers', 
                                            'Boreal Forests', 'Tundra'))
deficit_eco2    = deficit_eco2[!is.na(deficit_eco2$Percent),]
deficit_eco_sum = summarySE(deficit_eco2, measurevar="Percent", groupvars=c("BIOME_NAME","DeficitRestoration"))
names(deficit_eco_sum)[2] = 'group'
ggplot(deficit_eco_sum, aes(x=BIOME_NAME, y=Percent, fill=group))  + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=Percent-se, ymax=Percent+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  xlab("") + theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylab("Deficit restoration percentage (%)") +  theme_bw() + 
  theme(axis.text.x = element_text(color="black", 
                                   size=8)) +
  scale_fill_manual(values=cbPalette)  ##save in 900*400

###calculate the functional composition: 
ungu_alienpoten_matrix = as.data.frame(values(all_ungupotential3))  ##this time is in 100km resolution
ungu_alienpoten_matrix = ungu_alienpoten_matrix[ , herbtrait.potential$Binomial]
#####get the traits:ungu_alienpoten_matrix
alienpotensite.bodymasstotal   =  (as.matrix(ungu_alienpoten_matrix) %*% herbtrait.potential$Mass.g )
alienpotenbodymasstotal.raster =  native_potenitalrichness ; values(alienpotenbodymasstotal.raster)  =  alienpotensite.bodymasstotal
alienpotenbodymasstotal.raster = alienpotenbodymasstotal.raster/1000
levelplot(alienpotenbodymasstotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the grazing behavior: 
alienpotensite.grazingtotal    =  as.matrix(ungu_alienpoten_matrix)    %*% herbtrait.potential$Diet.Graminoids 
alienpotengrazingtotal.raster  =  native_potenitalrichness ; values(alienpotengrazingtotal.raster)  =  alienpotensite.grazingtotal
levelplot(alienpotengrazingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the browsing behavior: 
alienpotensite.browsingtotal    = as.matrix(ungu_alienpoten_matrix) %*% herbtrait.potential$Diet.Browse.Fruit 
alienpotenbrowsingtotal.raster  =  native_potenitalrichness ; values(alienpotenbrowsingtotal.raster)  =  alienpotensite.browsingtotal
levelplot(alienpotenbrowsingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

#current alien species : 
ungu_currentalien_matrix = as.data.frame(values(alien.sp2))
names(ungu_currentalien_matrix) = c( 'Axis axis',  'Dama dama', 'Ovis orientalis aries', 
                                     'Bos taurus', 'Capra aegagrus hircus', 
                                     'Boselaphus tragocamelus', 'Antilope cervicapra', 'Oryx gazella',
                                     'Ammotragus lervia' , 
                                     'Equus ferus caballus', 'Equus africanus asinus', 'Sus scrofa',
                                     "Cervus nippon")
ungu_currentalien_matrix = ungu_currentalien_matrix[,c(herbtrait.currentalien$Binomial)]
currentaliensite.bodymasstotal   =  (as.matrix(ungu_currentalien_matrix) %*% herbtrait.currentalien$Mass.g )
currentalienbodymasstotal.raster =  native_potenitalrichness ; values(currentalienbodymasstotal.raster)  =  currentaliensite.bodymasstotal
currentalienbodymasstotal.raster =  currentalienbodymasstotal.raster/1000
levelplot(currentalienbodymasstotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the grazing behavior: 
currentaliensite.grazingtotal    =  as.matrix(ungu_currentalien_matrix)    %*% herbtrait.currentalien$Diet.Graminoids 
currentaliengrazingtotal.raster  =  native_potenitalrichness ; values(currentaliengrazingtotal.raster)  =  currentaliensite.grazingtotal
levelplot(currentaliengrazingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the browsing behavior: 
currentaliensite.browsingtotal    = as.matrix(ungu_currentalien_matrix) %*% herbtrait.currentalien$Diet.Browse.Fruit 
currentalienbrowsingtotal.raster  =  native_potenitalrichness ; values(currentalienbrowsingtotal.raster)  =  currentaliensite.browsingtotal
levelplot(currentalienbrowsingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

# current native
names(native.sp2) = c('Odocoileus virginianus', 'Odocoileus hemionus', 'Rangifer tarandus', 
                      'Cervus canadensis', 'Alces alces','Oreamnos americanus',
                      'Bison bison', 'Ovis canadensis' , 'Ovis dalli',
                      'Ovibos moschatus','Pecari tajacu', 'Antilocapra americana' )
ungu_currentnative_matrix = as.data.frame(values(native.sp2))
herbtrait.currentnative = herbtrait_use[rownames(herbtrait_use)%in% names(ungu_currentnative_matrix), ]
ungu_currentnative_matrix = ungu_currentnative_matrix[,herbtrait.currentnative$Binomial ]
currentnativesite.bodymasstotal   =  (as.matrix(ungu_currentnative_matrix) %*% herbtrait.currentnative$Mass.g )
currentnativebodymasstotal.raster =  native_potenitalrichness ; values(currentnativebodymasstotal.raster)  =  currentnativesite.bodymasstotal
currentnativebodymasstotal.raster = currentnativebodymasstotal.raster/1000
levelplot(currentnativebodymasstotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the grazing behavior: 
currentnativesite.grazingtotal    =  as.matrix(ungu_currentnative_matrix)    %*% herbtrait.currentnative$Diet.Graminoids 
currentnativegrazingtotal.raster  =  native_potenitalrichness ; values(currentnativegrazingtotal.raster)  =  currentnativesite.grazingtotal
levelplot(currentnativegrazingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the browsing behavior: 
currentnativesite.browsingtotal    = as.matrix(ungu_currentnative_matrix) %*% herbtrait.currentnative$Diet.Browse.Fruit 
currentnativebrowsingtotal.raster  =  native_potenitalrichness ; values(currentnativebrowsingtotal.raster)  =  currentnativesite.browsingtotal
levelplot(currentnativebrowsingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
currentallbodymasstotal.raster =  currentnativebodymasstotal.raster + currentalienbodymasstotal.raster
currentallgrazingtotal.raster  =  currentnativegrazingtotal.raster  + currentaliengrazingtotal.raster
currentallbrowsingtotal.raster =  currentnativebrowsingtotal.raster + currentalienbrowsingtotal.raster

###############potential native 
names(native_potenital) = stri_replace_all_fixed(names(native_potenital),  '_',  ' ')
ungu_nativepotential_matrix =  as.data.frame(values(native_potenital))
herbtrait.nativepotential   =  herbtrait_use[rownames(herbtrait_use)%in% names(ungu_nativepotential_matrix), ]
ungu_nativepotential_matrix = ungu_nativepotential_matrix[, herbtrait.nativepotential$Binomial]
nativepotentialsite.bodymasstotal   =  (as.matrix(ungu_nativepotential_matrix) %*% herbtrait.nativepotential$Mass.g )
nativepotentialbodymasstotal.raster =  native_potenitalrichness ; values(nativepotentialbodymasstotal.raster)  =  nativepotentialsite.bodymasstotal
nativepotentialbodymasstotal.raster = nativepotentialbodymasstotal.raster/1000
levelplot(nativepotentialbodymasstotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the grazing behavior: 
nativepotentialsite.grazingtotal    =  as.matrix(ungu_nativepotential_matrix)    %*% herbtrait.nativepotential$Diet.Graminoids 
nativepotentialgrazingtotal.raster  =  native_potenitalrichness ; values(nativepotentialgrazingtotal.raster)  =  nativepotentialsite.grazingtotal
levelplot(nativepotentialgrazingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))
##the browsing behavior: 
nativepotentialsite.browsingtotal    = as.matrix(ungu_nativepotential_matrix) %*% herbtrait.nativepotential$Diet.Browse.Fruit 
nativepotentialbrowsingtotal.raster  =  native_potenitalrichness ; values(nativepotentialbrowsingtotal.raster)  =  nativepotentialsite.browsingtotal
levelplot(nativepotentialbrowsingtotal.raster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

nativepotentialbodymass.raster =  nativepotentialbodymasstotal.raster + currentalienbodymasstotal.raster
nativepotentialgrazing.raster  =  nativepotentialgrazingtotal.raster  + currentaliengrazingtotal.raster
nativepotentialbrowsing.raster =  nativepotentialbrowsingtotal.raster + currentalienbrowsingtotal.raster

fullrangebodymasstotal.raster = nativepotentialbodymasstotal.raster + alienpotenbodymasstotal.raster
fullrangegrazingtotal.raster  = nativepotentialgrazingtotal.raster  + alienpotengrazingtotal.raster
fullrangebrowsingtotal.raster = nativepotentialbrowsingtotal.raster + alienpotenbrowsingtotal.raster

####we compare the following scenarios of trait composition: 
##  1. natural;   2.current native species;    3.current alien and native 4. current alien and potential native species; 5,  potential alien and native species 
bodaymass.allraster = c(naturalbodymasstotal.raster,     currentnativebodymasstotal.raster, 
                        currentallbodymasstotal.raster,  nativepotentialbodymass.raster,
                        fullrangebodymasstotal.raster)
names(bodaymass.allraster) = c('natural','current_native','current_all','native_potential','all_potential')
levelplot(bodaymass.allraster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

grazing.allraster  =  c(naturalgrazingtotal.raster,   currentnativegrazingtotal.raster,
                        currentallgrazingtotal.raster, nativepotentialgrazing.raster,
                        fullrangegrazingtotal.raster)
names(grazing.allraster) = c('natural','current_native','current_all','native_potential','all_potential')
levelplot(grazing.allraster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

browsing.allraster =  c(naturalbrowsingtotal.raster,  currentnativebrowsingtotal.raster,
                        currentallbrowsingtotal.raster, nativepotentialbrowsing.raster,
                        fullrangebrowsingtotal.raster)
names(browsing.allraster) = c('natural','current_native','current_all','native_potential','all_potential')
levelplot(browsing.allraster, margin = FALSE, par.settings = BuRdTheme, scales=list(draw=FALSE))

##give the deficit map, biomass: 
deficit_biomass_currentnative  = naturalbodymasstotal.raster - currentnativebodymasstotal.raster  ##current native deficit 
deficit_biomass_currentnative[deficit_biomass_currentnative<0] = NA
deficit_biomass_currentnative2 = deficit_biomass_currentnative/naturalbodymasstotal.raster
deficit_biomass_currentnative2[is.infinite(deficit_biomass_currentnative2)] = NA
plot(deficit_biomass_currentnative)
mean(values(deficit_biomass_currentnative2),na.rm=TRUE) ##0.9431

deficit_biomass_currentall    = naturalbodymasstotal.raster - currentallbodymasstotal.raster
biomass_restoration_currentall = (currentallbodymasstotal.raster - currentnativebodymasstotal.raster)/deficit_biomass_currentnative
biomass_restoration_currentall[biomass_restoration_currentall<0] = 0
plot(biomass_restoration_currentall)
mean(values(biomass_restoration_currentall),na.rm=TRUE) #0.00432

nativepotenital_deficit = naturalbodymasstotal.raster - nativepotentialbodymass.raster
plot(nativepotenital_deficit)
biomass_restoration_nativepoten = (deficit_biomass_currentnative - nativepotenital_deficit)/deficit_biomass_currentnative
biomass_restoration_nativepoten[biomass_restoration_nativepoten<0] = 0
biomass_restoration_nativepoten[biomass_restoration_nativepoten>1] = 1
plot (biomass_restoration_nativepoten)
mean (values(biomass_restoration_nativepoten),na.rm=TRUE) ##0.157

biomass_restoration_alienpoten = (fullrangebodymasstotal.raster - currentnativebodymasstotal.raster)/deficit_biomass_currentnative
biomass_restoration_alienpoten[is.infinite(biomass_restoration_alienpoten)] = NA
biomass_restoration_alienpoten[biomass_restoration_alienpoten<0] = 0
biomass_restoration_alienpoten[biomass_restoration_alienpoten>1] = 1
plot(biomass_restoration_alienpoten)
mean (values(biomass_restoration_alienpoten),na.rm=TRUE)  ##0.186

biomass_restoration = c (biomass_restoration_currentall, 
                         biomass_restoration_nativepoten, 
                         biomass_restoration_alienpoten)
names(biomass_restoration ) =  c( 'Current_Alien', 'Potential_Native','Potential_Native_Alien')
levelplot(biomass_restoration*100, margin = FALSE,  par.strip.text = list(cex = 0),
          par.settings = list(axis.line = list(col = "transparent"), 
                              strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c('blue', "#fed976", "red"))(100),
          scales=list(draw=FALSE))

####we also compare between ecoregions: 
deficit_mass_eco =  c(biomass_restoration*100, ecoregion2 )
plot(deficit_mass_eco)
deficit_mass_ecodata = values(deficit_mass_eco, dataframe=TRUE)
deficit_mass_ecodata = deficit_mass_ecodata[!is.na(deficit_mass_ecodata$Current_Alien),]
table(deficit_mass_ecodata$BIOME_NAME)
###we may only analyse some biomes, like:
select_reg = c('Boreal Forests/Taiga', 'Deserts & Xeric Shrublands', 
               'Temperate Broadleaf & Mixed Forests', 
               'Temperate Conifer Forests', 'Temperate Grasslands, Savannas & Shrublands', 'Tundra')
deficit_mass_eco = deficit_mass_ecodata[deficit_mass_ecodata$BIOME_NAME %in% select_reg, ]
deficit_mass_eco2 = gather(deficit_mass_eco, 'DeficitRestoration','Percent',1:3)
ggplot(deficit_mass_eco2, aes(x=BIOME_NAME, y=Percent, fill=DeficitRestoration)) + geom_boxplot() 
class(deficit_mass_eco2$BIOME_NAME)
deficit_mass_eco2$BIOME_NAME = as.character(deficit_mass_eco2$BIOME_NAME)
deficit_mass_eco2$BIOME_NAME [deficit_mass_eco2$BIOME_NAME =='Temperate Broadleaf & Mixed Forests'] = 'Temperate Broadleaf'
deficit_mass_eco2$BIOME_NAME [deficit_mass_eco2$BIOME_NAME =='Boreal Forests/Taiga'] = 'Boreal Forests'
deficit_mass_eco2$BIOME_NAME [deficit_mass_eco2$BIOME_NAME =='Temperate Grasslands, Savannas & Shrublands'] = 'Grasslands'
deficit_mass_eco2$BIOME_NAME [deficit_mass_eco2$BIOME_NAME =='Deserts & Xeric Shrublands'] = 'Deserts'
deficit_mass_eco2$BIOME_NAME [deficit_mass_eco2$BIOME_NAME =='Temperate Conifer Forests'] = 'Temperate Conifers'

deficit_mass_eco2$BIOME_NAME = factor(deficit_mass_eco2$BIOME_NAME, 
                                      levels = c('Deserts','Grasslands',
                                                 'Temperate Broadleaf', 'Temperate Conifers', 
                                                 'Boreal Forests', 'Tundra'))
deficit_mass_eco2    = deficit_mass_eco2[!is.na(deficit_mass_eco2$Percent),]
deficit_mass_eco_sum = summarySE(deficit_mass_eco2, measurevar="Percent", groupvars=c("BIOME_NAME","DeficitRestoration"))
names(deficit_mass_eco_sum)[2] = 'group'
ggplot(deficit_mass_eco_sum, aes(x=BIOME_NAME, y=Percent, fill=group))  + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=Percent-se, ymax=Percent+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  xlab("") + theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylab("Deficit restoration percentage (%)") +  theme_bw() + 
  theme(axis.text.x = element_text(color="black", size=8)) +
  scale_fill_manual(values=cbPalette)
##save in 900*400


#---------------------------------------------------------------------
##give the deficit map, grazing: 
deficit_grazing_currentnative = naturalgrazingtotal.raster  - currentnativegrazingtotal.raster
plot(deficit_grazing_currentnative)
deficit_grazing_currentnative[deficit_grazing_currentnative<0]  = NA
deficit_grazing_currentnative2 = deficit_grazing_currentnative/naturalgrazingtotal.raster
mean(values(deficit_grazing_currentnative2), na.rm= TRUE)  #0.879

deficit_grazing_currentall     = naturalgrazingtotal.raster  - currentallgrazingtotal.raster
grazing_restoration_currentall = (currentallgrazingtotal.raster - currentnativegrazingtotal.raster)/deficit_grazing_currentnative
grazing_restoration_currentall[grazing_restoration_currentall<0] = 0
grazing_restoration_currentall[grazing_restoration_currentall>1] = 1
plot(grazing_restoration_currentall)
mean(values(grazing_restoration_currentall), na.rm= TRUE)  #0.063

grazing_restoration_nativepoten = (nativepotentialgrazing.raster - currentnativegrazingtotal.raster)/deficit_grazing_currentnative
grazing_restoration_nativepoten[grazing_restoration_nativepoten<0] = 0
grazing_restoration_nativepoten[grazing_restoration_nativepoten>1] = 1
plot(grazing_restoration_nativepoten)
mean(values(grazing_restoration_nativepoten), na.rm= TRUE)  #0.300

grazing_restoration_alienpoten = (fullrangegrazingtotal.raster - currentnativegrazingtotal.raster)/deficit_grazing_currentnative
grazing_restoration_alienpoten[grazing_restoration_alienpoten<0] = 0
grazing_restoration_alienpoten[grazing_restoration_alienpoten>1] = 1
plot(grazing_restoration_alienpoten)
mean(values(grazing_restoration_alienpoten), na.rm= TRUE)  #0.461


grazing_restoration = c(grazing_restoration_currentall, grazing_restoration_nativepoten, 
                        grazing_restoration_alienpoten)
names(grazing_restoration ) =  c( 'Current_Alien', 'Potential_Native','Potential_Native_Alien')
levelplot(grazing_restoration*100, margin = FALSE,  par.strip.text = list(cex = 0),
          par.settings = list(axis.line = list(col = "transparent"), 
                              strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c('blue', "#fed976", "red"))(100),
          scales=list(draw=FALSE))  ###save in plot 600*200

####we also compare between ecoregions: 
deficit_grazing_eco =  c(grazing_restoration*100, ecoregion2 )
plot(deficit_grazing_eco)
deficit_grazing_ecodata = values(deficit_grazing_eco, dataframe=TRUE)
deficit_grazing_ecodata = deficit_grazing_ecodata[!is.na(deficit_grazing_ecodata$Current_Alien),]
table(deficit_grazing_ecodata$BIOME_NAME)
###we may only analyse some biomes, like:
select_reg = c('Boreal Forests/Taiga', 'Deserts & Xeric Shrublands', 
               'Temperate Broadleaf & Mixed Forests', 
               'Temperate Conifer Forests', 'Temperate Grasslands, Savannas & Shrublands', 'Tundra')
deficit_grazing_eco = deficit_grazing_ecodata[deficit_grazing_ecodata$BIOME_NAME %in% select_reg, ]
table(deficit_grazing_ecodata$BIOME_NAME)
deficit_grazing_eco2 = gather(deficit_grazing_eco, 'DeficitRestoration','Percent',1:3)
ggplot(deficit_grazing_eco2, aes(x=BIOME_NAME, y=Percent, fill=DeficitRestoration)) + geom_boxplot() 
class(deficit_grazing_eco2$BIOME_NAME)
deficit_grazing_eco2$BIOME_NAME = as.character(deficit_grazing_eco2$BIOME_NAME)
deficit_grazing_eco2$BIOME_NAME [deficit_grazing_eco2$BIOME_NAME =='Temperate Broadleaf & Mixed Forests'] = 'Temperate Broadleaf'
deficit_grazing_eco2$BIOME_NAME [deficit_grazing_eco2$BIOME_NAME =='Boreal Forests/Taiga'] = 'Boreal Forests'
deficit_grazing_eco2$BIOME_NAME [deficit_grazing_eco2$BIOME_NAME =='Temperate Grasslands, Savannas & Shrublands'] = 'Grasslands'
deficit_grazing_eco2$BIOME_NAME [deficit_grazing_eco2$BIOME_NAME =='Deserts & Xeric Shrublands'] = 'Deserts'
deficit_grazing_eco2$BIOME_NAME [deficit_grazing_eco2$BIOME_NAME =='Temperate Conifer Forests'] = 'Temperate Conifers'

deficit_grazing_eco2$BIOME_NAME = factor(deficit_grazing_eco2$BIOME_NAME, 
                                         levels = c('Deserts','Grasslands',
                                                    'Temperate Broadleaf', 'Temperate Conifers', 
                                                    'Boreal Forests', 'Tundra'))
deficit_grazing_eco2    = deficit_grazing_eco2[!is.na(deficit_grazing_eco2$Percent),]
deficit_grazing_eco_sum = summarySE(deficit_grazing_eco2, measurevar="Percent", groupvars=c("BIOME_NAME","DeficitRestoration"))
names(deficit_grazing_eco_sum)[2] = 'group'
ggplot(deficit_grazing_eco_sum, aes(x=BIOME_NAME, y=Percent, fill=group))  + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=Percent-se, ymax=Percent+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  xlab("") + theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylab("Deficit restoration percentage (%)") +  theme_bw() + 
  theme(axis.text.x = element_text(color="black", size=8)) +
  scale_fill_manual(values=cbPalette)


##give the deficit map, browsing: 
deficit_browsing_currentnative = naturalbrowsingtotal.raster  - currentnativebrowsingtotal.raster
deficit_browsing_currentnative[deficit_browsing_currentnative<0] =NA
deficit_browsing_currentnative2 = deficit_browsing_currentnative/ naturalbrowsingtotal.raster
#deficit_browsing_currentnative2[deficit_browsing_currentnative2<0] = 0
mean(values(deficit_browsing_currentnative2), na.rm=TRUE); #0.856

deficit_browsing_currentall     = naturalbrowsingtotal.raster  - currentallbrowsingtotal.raster
browsing_restoration_currentall = (currentallbrowsingtotal.raster - currentnativebrowsingtotal.raster)/deficit_browsing_currentnative
browsing_restoration_currentall[browsing_restoration_currentall<0] = 0
browsing_restoration_currentall[browsing_restoration_currentall>1] = 1
plot(browsing_restoration_currentall)
mean(values(browsing_restoration_currentall), na.rm=TRUE)  #0.047

browsing_restoration_nativepoten = (nativepotentialbrowsing.raster - currentnativebrowsingtotal.raster)/deficit_browsing_currentnative
browsing_restoration_nativepoten[browsing_restoration_nativepoten<0] = 0
browsing_restoration_nativepoten[browsing_restoration_nativepoten>1] = 1
plot(browsing_restoration_nativepoten)
mean(values(browsing_restoration_nativepoten), na.rm=TRUE)  #0.3209

browsing_restoration_alienpoten = (fullrangebrowsingtotal.raster - currentnativebrowsingtotal.raster)/deficit_browsing_currentnative
browsing_restoration_alienpoten[is.infinite(browsing_restoration_alienpoten)] = 0
browsing_restoration_alienpoten[browsing_restoration_alienpoten<0] = 0
browsing_restoration_alienpoten[browsing_restoration_alienpoten>1] = 1
plot(browsing_restoration_alienpoten)
mean(values(browsing_restoration_alienpoten), na.rm=TRUE)###41.5%

browsing_restoration = c(browsing_restoration_currentall, browsing_restoration_nativepoten, 
                         browsing_restoration_alienpoten)
names(browsing_restoration ) =  c( 'Current_Alien', 'Potential_Native','Potential_Native_Alien')

levelplot(browsing_restoration*100, margin = FALSE,  par.strip.text = list(cex = 0),
          par.settings = list(axis.line = list(col = "transparent"), 
                              strip.background = list(col = 'transparent'), 
                              strip.border = list(col = 'transparent')), 
          col.regions = colorRampPalette(c('blue', "#fed976", "red"))(100),
          scales=list(draw=FALSE))


####we also compare between ecoregions: 
deficit_browsing_eco = c(browsing_restoration*100, ecoregion2 )
plot(deficit_browsing_eco)
deficit_browsing_ecodata = values(deficit_browsing_eco, dataframe=TRUE)
deficit_browsing_ecodata = deficit_browsing_ecodata[!is.na(deficit_browsing_ecodata$Current_Alien),]
table(deficit_browsing_ecodata$BIOME_NAME)
###we may only analyse some biomes, like:
select_reg = c('Boreal Forests/Taiga', 'Deserts & Xeric Shrublands', 
               'Temperate Broadleaf & Mixed Forests', 
               'Temperate Conifer Forests', 'Temperate Grasslands, Savannas & Shrublands', 'Tundra')
deficit_browsing_eco = deficit_browsing_ecodata[deficit_browsing_ecodata$BIOME_NAME %in% select_reg, ]
table(deficit_browsing_ecodata$BIOME_NAME)
deficit_browsing_eco2 = gather(deficit_browsing_eco, 'DeficitRestoration','Percent',1:3)
ggplot(deficit_browsing_eco2, aes(x=BIOME_NAME, y=Percent, fill=DeficitRestoration)) + geom_boxplot() 
class(deficit_browsing_eco2$BIOME_NAME)
deficit_browsing_eco2$BIOME_NAME = as.character(deficit_browsing_eco2$BIOME_NAME)
deficit_browsing_eco2$BIOME_NAME [deficit_browsing_eco2$BIOME_NAME =='Temperate Broadleaf & Mixed Forests'] = 'Temperate Broadleaf'
deficit_browsing_eco2$BIOME_NAME [deficit_browsing_eco2$BIOME_NAME =='Boreal Forests/Taiga'] = 'Boreal Forests'
deficit_browsing_eco2$BIOME_NAME [deficit_browsing_eco2$BIOME_NAME =='Temperate Grasslands, Savannas & Shrublands'] = 'Grasslands'
deficit_browsing_eco2$BIOME_NAME [deficit_browsing_eco2$BIOME_NAME =='Deserts & Xeric Shrublands'] = 'Deserts'
deficit_browsing_eco2$BIOME_NAME [deficit_browsing_eco2$BIOME_NAME =='Temperate Conifer Forests'] = 'Temperate Conifers'

deficit_browsing_eco2$BIOME_NAME = factor(deficit_browsing_eco2$BIOME_NAME, 
                                          levels = c('Deserts','Grasslands',
                                                     'Temperate Broadleaf', 'Temperate Conifers', 
                                                     'Boreal Forests', 'Tundra'))
deficit_browsing_eco2    = deficit_browsing_eco2[!is.na(deficit_browsing_eco2$Percent),]
deficit_browsing_eco_sum = summarySE(deficit_browsing_eco2, measurevar="Percent", groupvars=c("BIOME_NAME","DeficitRestoration"))

names(deficit_browsing_eco_sum)[2] = 'group'
ggplot(deficit_browsing_eco_sum, aes(x=BIOME_NAME, y=Percent, fill=group))  + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=Percent-se, ymax=Percent+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  xlab("") + theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylab("Deficit restoration percentage (%)") +  theme_bw() + 
  theme(axis.text.x = element_text(color="black", size=8)) +
  scale_fill_manual(values=cbPalette) 
