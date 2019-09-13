# this function takes a distribution of elevations by species. It then takes the distribution of points,
# computes kernel density estimates using the stats density() function and using these densities, calculates 
# the percentage overlap in elevation
# input takes the form of a dataframe of unique records
# with the first column as taxon names and the second column as the mean/median elevation for that collection
# it also removes outlier points to account for errors in data collection, such as labels for an entire site 
# with a wide range in elevations, so the 25th and 75th quantiles are excluded, 

`elevation.overlap` <-
function(data)
{

# get all combinations of taxa
combos_elev <- combn(unique(data$species),2, simplify=T)

# get min and max elevations for the whole data set with a 100m buffer for the upper limit
lower <- min(data$elevation)
upper <- max(data$elevation) + 100

elev_overlap<-list(); elev_i<-list(); elev_j<-list()
for(k in seq_along(combos_elev[1,])){

# extract ith and jth taxon in each pair	
	
  i <- combos_elev[1,k]
  elev_i[k]<- i
  j <- combos_elev[2,k]
  elev_j[k]<- j
  
# take just ith species for now
  elev_int_1 <- data[data$species == i,]

# identify outliers
  outliers_elev1 <- boxplot(elev_int_1$elevation, plot=FALSE)$out

# if no outliers, proceed. Otherwise, exclude them 
  ifelse(length(outliers_elev1) > 0, elev_int_1 <- elev_int_1[-which(elev_int_1$elevation %in%  outliers_elev1),],elev_int_1)

# for singleton records, create 2 points in a +/-100m buffer
  ifelse(nrow(elev_int_1) == 1,elev_int_1 <- rbind.data.frame(elev_int_1, c(elev_int_1[1,1],elev_int_1[1,2]+100),c(elev_int_1[1,1],elev_int_1[1,2]-100)),elev_int_1)
  elev_int_1$elevation <- as.numeric(elev_int_1$elevation)
  
# do the same for the jth species  
  elev_int_2 <- data[data$species == j,]
  outliers_elev2 <- boxplot(elev_int_2$elevation, plot=FALSE)$out
  ifelse(length(outliers_elev2) > 0,elev_int_2 <- elev_int_2[-which(elev_int_2$elevation %in%  outliers_elev2),],elev_int_2)
  ifelse(nrow(elev_int_2) == 1,elev_int_2 <- rbind.data.frame(elev_int_2, c(elev_int_2[1,1],elev_int_2[1,2]+100),c(elev_int_2[1,1],elev_int_2[1,2]-100)),elev_int_2)
  elev_int_2$elevation <- as.numeric(elev_int_2$elevation)
 
# generate kernel density estimates for each species in the pair 
  
  da = density(elev_int_1$elevation,from = lower, to = upper, adjust = 1)
  db = density(elev_int_2$elevation, from = lower, to = upper, adjust = 1)
  
# calculate the true coefficient of overlapping between two distributions (fraction of  area under the lower of the two curves).   
  elev_frac_overlap <- overlapTrue(da$y,db$y)
  elev_overlap[k]<- elev_frac_overlap
  }
  
# get result into a good format

	diff_elev <- data.frame(matrix(unlist(elev_overlap)))
	diff_elev <- cbind.data.frame(unlist(elev_i),unlist(elev_j),diff_elev)
	colnames(diff_elev) <- c("sp1","sp2","elevation_overlap")

# make sure this can be merged with other datasets so reverse order of species and append
	diff_elev2 <- diff_elev[,c(2,1,3)]; colnames(diff_elev2)<-c("sp1","sp2","elevation_overlap")
	diff_elev<-rbind(diff_elev,diff_elev2)
	return(diff_elev)
}