# this function calculates overlap for a series of hulls
# it requires a list of hulls and a datafame of areas for each hull
# it also calculates sypmatry as a function of the area of intersection divided by the area of the less wide ranging species
# it requires a list of hulls as a SpatialPolygons object

`hulls.sympatry` <-
 function(hulls, areas){
	combos <- combn(names(hulls),2, simplify=T)

# set up lists 	
	overlap_ith <- list(); overlap_jth <- list()
	areas_intersect <- list()
	centroid_lat <- list(); centroid_lon <- list()
	sp_i <- list(); sp_j <- list()

	for(k in seq_along(combos[1,])){
  		i <- combos[1,k]
  		j <- combos[2,k]
  		
  		area_i <- areas[which(areas$species == i),1]
  		area_j <- areas[which(areas$species == j),1]

# calculate polygon of overlap
  		inter <- gIntersection(hulls[[which(names(hulls) == i)]], 
                        	 	hulls[[which(names(hulls) == j)]],
                       		    byid = FALSE)
# calculate area of overlapping polygon in km2. If NA (polygons do not intersect) then return 0  
  	if(is.null(inter)){
    	areas2<-0
    	center_inter<-data.frame(NA,NA)
  										}
  	else{
    	areas2 <- areaPolygon(inter)/1e6
    	center_inter <- gCentroid(inter)
    	center_inter <- data.frame(center_inter@coords)
 		
 		areas2 <- areas2/pmin(area_i,area_j)		
 										}
  		areas_intersect[k] <- areas2
  		sp_i[k] <- i
  		sp_j[k] <- j
  		
  		centroid_lon[k]<-center_inter[1]
  		centroid_lat[k]<-center_inter[2]
											}
											
	symp_vals <- data.frame(matrix(unlist(areas_intersect)))
	sp_i <- data.frame(matrix(unlist(sp_i)))
	sp_j <- data.frame(matrix(unlist(sp_j)))
	symp_vals <- cbind.data.frame(sp_i,sp_j,symp_vals);colnames(symp_vals)=c("sp1","sp2","SYMP")	
	return(symp_vals) 
 }
  
  
  
  