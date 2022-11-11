
library(ape)
library(phytools)
library(picante)
library(ape)
library(picante)
library(ggplot2)
library(ggtree)
library(dplyr)
library(raster)
library(data.table)
library(foreach)
library(doParallel)
library(rgdal)
library(terra)
require(CAST)

##### future PD calculation ####
dis.future<-get(load("./distribution/future_dis_stable/RCP85/proj_2070_rcp8.5.RData"))
tree.100<-get(load("Random100.RData"))
names_match<-read.csv("names_match.csv",row.names=1)
colnames(dis.future)<-names_match$name
tree.exp<-tree.100[[1]]
spnames<-colnames(dis.future)[which(colnames(dis.future) %in% tree.exp$tip.label==T)]
dis.future.1<-dis.future[,spnames]


##### make cluster ###
library(parallel)
pd.fun<-function(x) pd1<-pd(spdis,x,include.root=TRUE)
spdis<-dis.future.1
cl<-makeCluster(100)
clusterExport(cl = cl, varlist = c("spdis"))
clusterEvalQ(cl = cl, library(picante))
pd.obs<-parLapply(cl,tree.100,pd.fun)
stopCluster(cl)
### save(pd.obs,file="PD/ensemble/pd.obs.2070.rcp8.5.RData")

#### branch analyses ###
source("phylo_branch_matrix.R")

# Prepare basic shapefile  
basic_shp <- shapefile('./grid20km/grid20kmChina.shp')
basic_shp$GRIDCODE <- as.numeric(basic_shp$GRIDCODE)

# Prepare point to extract values from shapefiles
basic_raster <- raster('./raster/China_r.tif')
coord <- coordinates(basic_raster) %>% data.frame()
coord <- coord[!is.na(values(basic_raster)),]
coordinates(coord) <- ~x+y
crs(coord) <- crs(basic_shp)

#######Processing phylogenetic data for many trees #######
####
setwd("D:/Shijia/Tree_of_life/Zonation/Spatial_analysis")
dir.create('./full/GeoPhylo_rcp85') #Create a directory where will be stored geographical expression of phylogenetic branches for each one of the 100 phylogenetic trees.
##### Get data.frame with distribution of branches for each row (i.e. cell)
for(i in 1:length(tree)) {
  message(i)
  dirs <- file.path('./full/GeoPhylo_rcp85', paste0('Phylotree_', i))
  dir.create(dirs)
  m <-phylo_branch_matrix(site_sp_matrix = data_dis_1, tree_file = tree[[i]])
  data.table::fwrite(m, file.path(dirs, paste0("phylo_", i, '.gz')))
}

## Transform each matrix to raster 
mlist <- list.dirs('./full/GeoPhylo_rcp85') %>% grep('Phylotree_',., value = T) %>% list.files(pattern = '.gz$', full.names = TRUE)

for (i in 1:length(mlist)){
  message(i)
  dirs <- file.path('./full/GeoPhylo_rcp85', paste0('Phylotree_', i))
  m <- data.table::fread(mlist[i]) %>% data.frame()
  basic_shp2 <- basic_shp["GRIDCODE"]
  basic_shp2@data <-left_join(basic_shp2@data, m, by = c("GRIDCODE" = 'ncell'))
  basic_shp2$"GRIDCODE" <- NULL
  
  ####Extract values for each raster cell from shapefile
  
  layers <- names(basic_shp2)
  layers <- gsub('X', '', layers)
  cl <- parallel::makeCluster(30)
  registerDoParallel(cl)
  foreach(ii = 1:ncol(basic_shp2),
          .packages = c('raster')) %dopar% {
            basic_raster2 <- basic_raster
            values(basic_raster2)[!is.na(values(basic_raster))] <-
              over(coord, basic_shp2[ii])[, 1]
            raster::writeRaster(basic_raster2,filename=paste(layers[ii],".tif",seq=""))
          }
  
  stopCluster(cl)
}

#### delete useless files
allfile = dir() 
txtfile <- grep("*.tif.aux.xml", allfile)
file.remove(allfile[txtfile])

#####Branches weights#####

branches<-cbind(tree_exp$edge, data.frame(tree_exp$edge.length))
colnames(branches)[3]<-"edge.length"
path<-"./GeoPhylo/Phylotree_1"
dat <- data.frame(list.files(path, pattern=".tif"))
colnames(dat)<-"edge"
dat$edge <- gsub(".tif", "", dat$edge, fixed = T)
p<-merge(branches, dat, by.x = '2', by.y = 'edge', all.x = T)[,c(1,3)]
write.table(p,file = file.path(d, "branch_weight.spp"),col.names = FALSE,row.names = FALSE,sep = "\t",quote = FALSE)

###### Gap_analysis ######################################################

#### Species gap ###
dis_current <- get(load("./distribution/current/proj_current_total.RData")) 

### RCP85 as an example ###
dis_full <- get(load("./distribution/future_dis_fully/RCP85/proj_RCP85_total.RData"))
dis_buffer <- get(load("./distribution/future_dis_20km_buffer/RCP85/proj_RCP85_total.RData"))
dis_stable <- get(load("./distribution/future_dis_stable/RCP85/proj_RCP85_total.RData"))

##################
dis_current <- data.frame(dis_current)
spp <- colnames(dis_current)
dis_current_1 <- dis_current
dis_current_1 <- dis_current_1 %>% dplyr::mutate(GRIDCODE = as.numeric(rownames(dis_current_1))) %>% dplyr::relocate(GRIDCODE) %>% tibble()

china_grid <- terra::vect("./Gap_analysis/grid20kmChina_with_PAs/grid20kmChina_with_PAs.shp")
china_grid2 <- as.data.frame(china_grid) %>% dplyr::select("ID", "GRIDCODE", "AREA", "LON", "LAT", "PAcovarage")
names(china_grid2)
dis_current_1 <- left_join(dis_current_1, china_grid2) %>% dplyr::relocate(colnames(china_grid2))

sp_pas <- tibble(sp = spp, total_area = 0, within_pa = 0)
cell_area <- dis_current_1$AREA
PAcovarage <- dis_current_1$PAcovarage
for (i in 1:length(spp)) {
  sp_name <- spp[i]
  dist <- dis_current_1 %>% pull(sp_name)
  sp_pas[sp_pas$sp == sp_name, "total_area"] <- sum(dist * cell_area)
  sp_pas[sp_pas$sp == sp_name, "within_pa"] <- sum(dist * cell_area * PAcovarage)
}

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

target <- na.omit(sp_pas$total_area) %>% unique() %>% sort()
target <- tibble(total_area = c(1, 1000, target, 250000) %>% sort(), target = NA)
target$target[target$total_area <= 1000] <- 1
target$target[target$total_area >= 250000] <- 0.10

da_model <- target %>% dplyr::filter(total_area %in% c(1000,250000))
inter <-approx(da_model$total_area, da_model$target, target$total_area[is.na(target$target)])
target$target[is.na(target$target)] <- inter$y

sp_pas <- left_join(sp_pas, target)
sp_pas$target_achieved <- (sp_pas$proportion) / (sp_pas$target)
sp_pas$category <-cut(sp_pas$target_achieved,breaks = c(-Inf, 0, 0.2, .9, Inf),labels = c("No protected", "Gap", "Partial gap", "Protected"))


##################
dis_future <- data.frame(dis_stable)
spp <- colnames(dis_future)
dis_future_1 <- dis_future
dis_future_1 <- dis_future_1 %>% dplyr::mutate(GRIDCODE = as.numeric(rownames(dis_future_1))) %>% dplyr::relocate(GRIDCODE) %>% tibble()

china_grid <- terra::vect("./Gap_analysis/grid20kmChina_with_PAs/grid20kmChina_with_PAs.shp")
china_grid2 <- as.data.frame(china_grid) %>% dplyr::select("ID", "GRIDCODE", "AREA", "LON", "LAT", "PAcovarage")
names(china_grid2)
dis_future_1 <- left_join(dis_future_1, china_grid2) %>% dplyr::relocate(colnames(china_grid2))

sp_pas <- tibble(sp = spp, total_area = 0, within_pa = 0)
cell_area <- dis_future_1$AREA
PAcovarage <- dis_future_1$PAcovarage
for (i in 1:length(spp)) {
  sp_name <- spp[i]
  dist <- dis_future_1 %>% pull(sp_name)
  sp_pas[sp_pas$sp == sp_name, "total_area"] <- sum(dist * cell_area)
  sp_pas[sp_pas$sp == sp_name, "within_pa"] <- sum(dist * cell_area * PAcovarage)
}

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

target <- na.omit(sp_pas$total_area) %>% unique() %>% sort()
target <- tibble(total_area = c(1, 1000, target, 250000) %>% sort(), target = NA)
target$target[target$total_area <= 1000] <- 1
target$target[target$total_area >= 250000] <- 0.10

da_model <- target %>% dplyr::filter(total_area %in% c(1000,250000))
inter <-approx(da_model$total_area, da_model$target, target$total_area[is.na(target$target)])
target$target[is.na(target$target)] <- inter$y

sp_pas <- left_join(sp_pas, target)
sp_pas$target_achieved <- (sp_pas$proportion) / (sp_pas$target)
sp_pas$category <-cut(sp_pas$target_achieved,breaks = c(-Inf, 0, 0.2, .9, Inf),labels = c("No protected", "Gap", "Partial gap", "Protected"))
write.csv(sp_pas,file="./Gap_analysis/Species/RCP26/Gap_stable.csv")


#### branch gap ###################################
dis_current <- get(load("./SR_additional/branch/current_branches.RData")) 

### RCP85 as an example ###
dis_full <- get(load("./SR_additional/branch/RCP85/full_branch.RData"))
dis_buffer <- get(load("./SR_additional/branch/RCP85/buffer_branch.RData"))
dis_stable <- get(load("./SR_additional/branch/RCP85/stable_branch.RData"))

##################
dis_current <- data.frame(dis_current)
spp <- colnames(dis_current)
dis_current_1 <- dis_current
dis_current_1 <- dis_current_1 %>% dplyr::mutate(GRIDCODE = as.numeric(rownames(dis_current_1))) %>% dplyr::relocate(GRIDCODE) %>% tibble()

china_grid <- terra::vect("./Gap_analysis/grid20kmChina_with_PAs/grid20kmChina_with_PAs.shp")
china_grid2 <- as.data.frame(china_grid) %>% dplyr::select("ID", "GRIDCODE", "AREA", "LON", "LAT", "PAcovarage")
names(china_grid2)
dis_current_1 <- left_join(dis_current_1, china_grid2) %>% dplyr::relocate(colnames(china_grid2))

sp_pas <- tibble(sp = spp, total_area = 0, within_pa = 0)
cell_area <- dis_current_1$AREA
PAcovarage <- dis_current_1$PAcovarage
for (i in 1:length(spp)) {
  sp_name <- spp[i]
  dist <- dis_current_1 %>% pull(sp_name)
  sp_pas[sp_pas$sp == sp_name, "total_area"] <- sum(dist * cell_area)
  sp_pas[sp_pas$sp == sp_name, "within_pa"] <- sum(dist * cell_area * PAcovarage)
}

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

target <- na.omit(sp_pas$total_area) %>% unique() %>% sort()
target <- tibble(total_area = c(1, 1000, target, 250000) %>% sort(), target = NA)
target$target[target$total_area <= 1000] <- 1
target$target[target$total_area >= 250000] <- 0.10

da_model <- target %>% dplyr::filter(total_area %in% c(1000,250000))
inter <-approx(da_model$total_area, da_model$target, target$total_area[is.na(target$target)])
target$target[is.na(target$target)] <- inter$y

sp_pas <- left_join(sp_pas, target)
sp_pas$target_achieved <- (sp_pas$proportion) / (sp_pas$target)
sp_pas$category <-cut(sp_pas$target_achieved,breaks = c(-Inf, 0, 0.2, .9, Inf),labels = c("No protected", "Gap", "Partial gap", "Protected"))

##################
dis_future <- data.frame(dis_full)
spp <- colnames(dis_future)
dis_future_1 <- dis_future
dis_future_1 <- dis_future_1 %>% dplyr::mutate(GRIDCODE = as.numeric(rownames(dis_future_1))) %>% dplyr::relocate(GRIDCODE) %>% tibble()

china_grid <- terra::vect("./Gap_analysis/grid20kmChina_with_PAs/grid20kmChina_with_PAs.shp")
china_grid2 <- as.data.frame(china_grid) %>% dplyr::select("ID", "GRIDCODE", "AREA", "LON", "LAT", "PAcovarage")
names(china_grid2)
dis_future_1 <- left_join(dis_future_1, china_grid2) %>% dplyr::relocate(colnames(china_grid2))

sp_pas <- tibble(sp = spp, total_area = 0, within_pa = 0)
cell_area <- dis_future_1$AREA
PAcovarage <- dis_future_1$PAcovarage
for (i in 1:length(spp)) {
  sp_name <- spp[i]
  dist <- dis_future_1 %>% pull(sp_name)
  sp_pas[sp_pas$sp == sp_name, "total_area"] <- sum(dist * cell_area)
  sp_pas[sp_pas$sp == sp_name, "within_pa"] <- sum(dist * cell_area * PAcovarage)
}

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

sp_pas$within_pa[is.na(sp_pas$within_pa)] <- 0
sp_pas <- sp_pas %>% mutate(proportion = within_pa / total_area)

target <- na.omit(sp_pas$total_area) %>% unique() %>% sort()
target <- tibble(total_area = c(1, 1000, target, 250000) %>% sort(), target = NA)
target$target[target$total_area <= 1000] <- 1
target$target[target$total_area >= 250000] <- 0.10

da_model <- target %>% dplyr::filter(total_area %in% c(1000,250000))
inter <-approx(da_model$total_area, da_model$target, target$total_area[is.na(target$target)])
target$target[is.na(target$target)] <- inter$y

sp_pas <- left_join(sp_pas, target)
sp_pas$target_achieved <- (sp_pas$proportion) / (sp_pas$target)
sp_pas$category <-cut(sp_pas$target_achieved,breaks = c(-Inf, 0, 0.2, .9, Inf),labels = c("No protected", "Gap", "Partial gap", "Protected"))

#### Environmental novelty based on AOA ####

current <- readr::read_csv("./Current_condition.csv")
future <- readr::read_csv("./future_rcp8.5.csv")

env_var <- c("bio4", "bio10", "bio11", "bio15", "bio16", "bio17") 
current2 <- current %>% dplyr::select(env_var)
future2 <- future %>% dplyr::select(env_var)

cl <- makeCluster(30)
registerDoParallel(cl)

aoa_metric <- CAST::aoa(newdata = future2, train=current2, variables="all", cl = cl)
stopCluster(cl)

aoa_metric$parameters # this thrshold was calculated based on 
# thres <- grDevices::boxplot.stats(TrainDI)$stats[5]
aoa_metric$DI #dissimilarity metric (environmental novelty degree)
aoa_metric$AOA #Area of Applicability 

novelity_future <- tibble(current %>% dplyr::select(Gridcode, LON, LAT), bind_rows(aoa_metric[c("DI", "AOA")]))










