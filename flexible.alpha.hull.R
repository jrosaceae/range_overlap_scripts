# this function calculates alpha hulls and associated areas for a set of points
# importantly, there is flexibility for cases with less than three points which
# is often the case for narrowly endemic species byt which the functions borrowed
# from other packages here do nor account for. 
# data should be of UNIQUE records in the form of a dataframe of all points, with the first column species,
# second column longitude, and third column latitude as is required by the used alpha hull
# function. IMPORTANT. coordinate columns should be "decimalLongitude" and "decimalLatitude"

options(warn = -1)

`flexible.alpha.hull` <-
function(data, singleton_buffer = 2414.02, max_distance = 3000)
{

# first remove all species with 2 unique records
	coo_2pts <- setDT(data)[, if (.N == 2) .SD, by = species]
	coo_2pts <- droplevels(coo_2pts)
	
# get distances between points
	dist2<-list()
	for (i in levels(coo_2pts)){
  		pts2 <- ruellia_2pts[which(ruellia_2pts$species == i),]
  		d2 <- distm(pts[1,2:3], pts[2,2:3], fun = distHaversine)
  		dist2[[i]] <- d2
  		}
	
	dist2<-do.call(rbind.data.frame, dist2)
	dist2$species <- rownames(dist2)
	coo_2pts <- merge(coo_2pts,dist2,by="species",all=T)
	
# remove species with points more than 3km away from each other, this is assumed to be bad taxonomy or sparse records
# select only odd rows for species close by (a random point). We are going to assume this is one locality and merge them with the next step

	coo_2pts <- coo_2pts[ which(coo_2pts$V1 <= max_distance), ]
	coo_2pts <- coo_2pts %>% dplyr::filter(row_number() %% 2 == 1)
	coo_2pts <- coo_2pts[,1:3]
	

# now remove all species with 1 unique record and merge with the data frame we just created (coo_2pts)

	coo_1pts <- setDT(data)[, if (.N == 1) .SD, by = species]
	coo_1pts <- rbind.data.frame(coo_1pts,coo_2pts)
	coo_1pts <- droplevels(coo_1pts)
	
# create a buffer around these points to generate a polygon that we can use downstream
# the buffer we have chosen here is 1.5 mi but this can be changed
	hull_buff<-list()
	for (i in levels(coo_1pts$species)){
	 sp1 <- coo_1pts[which(coo_1pts$species == i),]
	 point_obj<-SpatialPoints(sp1[1,2:3],proj4string=CRS("+proj=longlat +datum=WGS84"))
	 point_obj<-buffer(point_obj, width = singleton_buffer)
 	 spdf_obj <- SpatialPolygonsDataFrame(point_obj, sp1[1,], match.ID = F) 
 	 hull_buff[[i]]<-spdf_obj
 	 }

# calculate alpha hulls for all species with >= 3 records
# getDynamicAlphaHull throws WAY to many (repetitive) errors or provides inconsistent results
# as a workaround, this function tries multiple times. Twice to generate a hull that encompasses 0.95
# of points and if both fail, a hull that encompasses 0.90 of points. If all fail, NA is returned.
# higher values of alpha are more complex shapes a=0 is a convex hull!!

coo <- setDT(data)[, if (.N > 2) .SD, by = species]
coo <- droplevels(coo)
hulls <- list()

# calculate hulls loop

for (i in levels(coo$species)){
  sp1 <- coo[which(coo$species == i),]
  sp1 <- sp1[sample(nrow(sp1)),]
  sp1 <- droplevels(sp1)
  calc_hull <- try(getDynamicAlphaHull(sp1, coordHeaders = c('decimalLongitude',
                  'decimalLatitude'),initialAlpha = 0,verbose=F,fraction = 1,
                  partCount = 5,alphaIncrement=1))
                  
  if (class(calc_hull) == "try-error") {
    sp1 <- sp1[sample(nrow(sp1)),]
    calc_hull <-  try(getDynamicAlphaHull(sp1, coordHeaders = c('decimalLongitude',
                      'decimalLatitude'),initialAlpha = 0,verbose=F,fraction = 0.95,
                      partCount = 3,alphaIncrement=1))
                      
    if (class(calc_hull) == "try-error") {
      sp1 <- sp1[sample(nrow(sp1)),]
      calc_hull <-  try(getDynamicAlphaHull(sp1, coordHeaders = c('decimalLongitude',
                        'decimalLatitude'),initialAlpha = 0,verbose=F,fraction = 0.9,
                         partCount = 3,alphaIncrement=1))
                         
      if (class(calc_hull) == "try-error") {
        calc_hull <-  NA
      }}}
      
  calc_hull<-calc_hull[[1]]
  hulls[[i]]<-calc_hull
}
  
# merge both hull objects so all polygons are in one list
  hulls_clean<-hulls[!is.na(hulls)] 
  hulls_clean<-c(hulls_clean,hull_buff)

# calculate areas loop
	areas<-list()
	for (i in 1:length(hulls_clean)){
 	 area_temp<-hulls_clean[[i]]
 	 if(class(area_temp) == "character"){
     area_temp <- NA
  	}
  else{
    area_temp <- geosphere::areaPolygon(hulls_clean[[i]])
  }
# convert m2 to km2
  areas[[i]]<-area_temp/1e6
}

names(areas) = names(hulls_clean)
areas_df <- do.call(rbind.data.frame, areas)
rownames(areas_df) = names(areas)
colnames(areas_df) = c("area (km2)")
areas_df <- as.data.frame(unlist(areas))
areas_df$species <- rownames(areas_df)

return(list(hulls_clean,areas_df))
}
