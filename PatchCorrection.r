#########################################
############# PATCH CORRECTION ##########
#########################################

# Authors:
#----------
# Cornelius Senf, Technical University of Munich
# Nezha Acil, University of Birmingham & University of Leicester


rm(list=ls())

### Job array ###
#================

slurmId = Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(slurmId )


#### Libraries ####
#==================
library(tidyverse)
library(raster)
library(landscapemetrics)
library(data.table)
library(dplyr)



#### raster conversion to data.table ####
#========================================

as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = canProcessInMemory(x, 2), ...) {
  stopifnot(require("data.table"))
  if(inmem) {
    v <- as.data.table(as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
  } else {
    tr <- blockSize(x, n=2)
    l <- lapply(1:tr$n, function(w) 
      as.data.table(as.data.frame(getValues(x, 
                                            row=tr$row[w], 
                                            nrows=tr$nrows[w]), 
                                  row.names=row.names, optional=optional, xy=xy, ...)))
    v <- rbindlist(l)
  }
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  setnames(v, coln)
  v
}

#### Negating %in% ####
#======================
`%!in%` = Negate(`%in%`)


#### Loading disturbance years ####
#==================================

folder = "path to the UTM folder"

tifs = paste0(folder, "TIFs/")
setwd(tifs)

f = list.files(tifs, full.names = F, pattern="*.tif$", all.files=F)
print(f[[i]])

years = raster(f[[i]])



#### Looping through years and grouping adjacent patches of consecutive years ####
#====================================================================================

year_range <- 1:18 

for (y in year_range) {
    
print(y)
  

  # Checking if the raster contains pixels in 2 consecutive years (y and y+1)
  #-----------------------------------------------------------------------
  
if (y %!in% values(years) || (y+1) %!in% values(years)) next 
  print('There are pixels for both years y and y+1') 
  
  
  # Reclassifying consecutive years into binary raster
  #------------------------------------------------
  reclass_matrix <- matrix(c(year_range, year_range %in% c(y, y + 1)), ncol = 2) 
  selection <- reclassify(years, reclass_matrix, filename=paste0('selection_', substr(as.character(f[[i]]),1,3),y,'.tif'), overwrite=TRUE) 

  
  # Delineating patches of consecutive years
  #-----------------------------------------
  patches <- get_patches(selection, class = 1, directions = 8) 
  years_dt <- as.data.table.raster(years) 
  names(years_dt) <- "year"
  patches_dt <- as.data.table.raster(patches$`1`) 
  names(patches_dt) <- "patch"

  dt <- cbind(years_dt, patches_dt) 
  dt[patch>0 , n := length(unique(na.omit(year))), by = patch]
 
  # Checking patch adjacency
  #-------------------------   
  if (2 %!in% dt$n) next 
  print('Patches are adjacent')
  dt$year = as.integer(dt$year)
  
  
  # Patch metrics conditions (too slow)
  #-----------------------------------------
  
  #ci = spatialize_lsm(selection, what = "lsm_p_circle",directions = 8) circularity
  #fr = spatialize_lsm(selection, what = "lsm_p_frac",directions = 8) frac
  
  #c_dt <- as.data.table.raster(ci[[1]][[1]]) # add the corresponding CIRC index
  #names(c_dt) <- "circ"
  #f_dt <- as.data.table.raster(fr[[1]][[1]]) # add the corresponding FRAC index
  #names(f_dt) <- "frac"
  
  #dt <- cbind(years_dt, patches_dt, c_dt, f_dt) 
  
  
  # Re-assigning year to combined patches
  #--------------------------------------
  dt[n == 2, year := as.integer(modal(year, na.rm=T, ties='lowest')), by = patch] # 'lowest' means when there are ties, assign to earliest year
  #dt[n == 2 & circ > 0.6 & frac < 1.2, year := as.integer(modal(year)), by = patch] 
  
  
  # Overwriting disturbance years raster 
  #--------------------------------------
  values(years) <- dt$year 
  print('Merged')

}

# Masking no loss pixels and saving
#-----------------------------------
years=calc(years, fun = function(x) {x[x == 0] = NA; return(x)}, 
           filename=paste0("output folder", gsub(".tif", "_grp.tif", as.character(f[[i]]))), datatype="INT1U", overwrite=TRUE)

print(f[[i]])

