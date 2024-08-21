
slurmId = Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(slurmId )


# Loading libraries and modules
setwd("/rds/projects/2018/pughtam-acil-wind/Colleagues/Cornelius/Neighbourhood/Test_Neighbourhood/")
library(OSMscale)
library(raster)
library(tidyverse)
library(foreign)
library(dplyr)
library(RcppThread)
Rcpp::sourceCpp('ngb_rcpp.cpp')


`%!in%` = Negate(`%in%`)

##################
### Settings
##################

years <- 1:18
radius <- 5000 # in meters (10 000m, 5000m, 1000m)


##################
## Source raster
##################

tiles = '/rds/projects/p/pughtam-environ/Acil/Merged/TIFsGrp/Lowest/'
setwd(tiles)

f = list.files(tiles, full.names = F, pattern="*.tif$", all.files=F)
print(f[[i]])

tile = substr(as.character(f[[i]]),1,3)
dist_map = raster(f[[i]])
crs(dist_map)


##################
## Centroids
##################
l = list()
for (y in years){
  #sn = paste0('W:/Acil/Merged/TIFsGrp/Lowest/Annual/20',sprintf("%02d",y),'/PIDsGrp/SHP/',tile,'_UTM_30m_grp_20',sprintf("%02d",y),'_pid_lsm.dbf')
  print(y)
  sn = paste0('/rds/projects/p/pughtam-environ/Acil/Merged/TIFsGrp/Lowest/Annual/20',sprintf("%02d",y),'/PIDsGrp/SHP/',tile,'_UTM_30m_grp_20',sprintf("%02d",y),'_pid_lsm.dbf')
  if (!file.exists(sn)) next
  tb = read.dbf(sn)
  
  tb$Tile = tile#substr(as.character(x),1,3)
  tb$year=y
  #tb$year= sprintf("%02d", y)
  tb$PID = paste0(tb$Tile,'_',sprintf("%02d", y),'_',sprintf("%08d",tb$gridcode))
  tb$year = as.integer(tb$year)
  
  tb=tb[,c('PID', 'year','INSIDE_X', 'INSIDE_Y','X','Y','gridcode')]
  colnames(tb)[3]='X_inside_wgs84'
  colnames(tb)[4]='Y_inside_wgs84'
  colnames(tb)[5]='X_UTM'
  colnames(tb)[6]='Y_UTM'
  pts_utm = projectPoints(tb$Y_inside_wgs84, tb$X_inside_wgs84, to = crs(dist_map), spout = FALSE, dfout = TRUE)
  tb = cbind(tb, pts_utm)
  
  
  #rn = paste0('/rds/projects/p/pughtam-environ/Acil/Merged/TIFsGrp/Lowest/Annual/20',sprintf("%02d",y),'/PIDsGrp/',tile,'_UTM_30m_grp_20',sprintf("%02d",y),'_pid_lsm.tif')
  #rtb = raster(rn) 
  
  #index_values <- expand.grid(col = 1:ncol(rtb), row = 1:nrow(rtb))
  #index_values <- cbind(index_values, patch = getValues(rtb))
  
  te <- cbind(tb, rowColFromCell(dist_map, cellFromXY(dist_map, pts_utm)))
  #te = merge(tb,index_values,by.x='gridcode',by.y='patch', all.x=T)
  l[[y]] = te
  
}

patches_centroids = do.call(rbind, l)
rm(l,te,tb, sn)
str(patches_centroids)

##################
## Kernel
##################

maxd <- round(radius / 30)

dmat <- matrix(0, nrow= 2*maxd+1, ncol=2*maxd+1)

for (i in 1:(2 * maxd + 1))
  for (j in 1:(2 * maxd + 1))
    dmat[j, i] <- sqrt((maxd - i) * (maxd - i) + (maxd - j) * (maxd - j) ) * 30

dtab <- data.frame()

for (i in 1:(2 * maxd + 1))
  for (j in 1:(2 * maxd + 1))
    if (dmat[j, i] < radius)
      dtab <- rbind(dtab, data.frame(ix = j - maxd, iy = i - maxd))



##################
## Metric table
##################


dist_map_matrix <- as.matrix(dist_map)

ngb_out <- vector("list", length(years))

for (y in years) {
  
  print(y)
  
  #if (y %!in% values(dist_map)) next
  patches_centroids_tmp = patches_centroids[patches_centroids$year==y,]
  rn = paste0('/rds/projects/p/pughtam-environ/Acil/Merged/TIFsGrp/Lowest/Annual/20',sprintf("%02d",y),'/PIDsGrp/',tile,'_UTM_30m_grp_20',sprintf("%02d",y),'_pid_lsm.tif')
  if (!file.exists(rn)) next
  patches <- raster(rn)

  
  # -> This is the actual calculation, which is quite fast!
  
  ngb <- countPxCentroidMT(dist_map_matrix, 
                           patches_centroids_tmp$col, 
                           patches_centroids_tmp$row,
                           patches_centroids_tmp$year,
                           dtab$ix, 
                           dtab$iy,
                           getValues(patches))
  
  
  
  colnames(ngb) <- c("n_t0", # Disturbance area in year 0
                     "n_total", # Total disturbance area over the whole period
                     "n_tminus1", # Disturbance area in year -1
                     "n_tplus1", # Disturbance area in year +1
                     "n_patches") # Number of patches in year 0
  ngb <- cbind(patches_centroids_tmp, ngb)
  ngb_out[[y]] <- as.data.frame(ngb)
}

### Write into final dataframe and save

ngb_metrics <- ngb_out %>%
  bind_rows() %>%
  mutate(perc_t0 = n_t0 / n_total, # Percent of total disturbance area falling in year 0
         perc_tminus1 = n_tminus1 / n_total, # Percent of total disturbance area falling in year -1
         perc_tplus1 = n_tplus1 / n_total) %>% # Percent of total disturbance area falling in year +1
  as_tibble()

head(ngb_metrics)


write_csv(ngb_metrics, paste0("/rds/projects/p/pughtam-environ/Acil/Merged/TIFsGrp/Lowest/Annual/Ancillary_Tables/NGB/",as.character(radius),"/",tile, "_nbg_5000m.csv"))


