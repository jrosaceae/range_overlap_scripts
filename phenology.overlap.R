# this function takes a distribution of collection dates by species
# ideally, specimens can be sorted by their phenology (i.e. flowering, fruiting).
# it then takes the distribution of collection dates for UNIQUE collections,
# input takes the from a 3-4 column dataframe of unique records in the form of "species", "month", "day", "year"
# year is optinal, as its specific value is irrelvant for the calculation.
# this converts clander date to Julian day, conversts this to radians, and then claculates percent overlap
# using a function designed for time series (overlapEst) in the package overlap. This function calcuilates a 
# kernal density estimate. Based on the reccomendations of the finction, we calculate the Dhat1 statistic, which is eqivalet 
# to Dhat2 for circular distributions (as these are), and best for small sample sizes
# note that here, we are assuming that all specimens (post outlier removal) are flowering

`phenology.overlap` <-
function(data)
{

# first, fill in missing data (in various possible values). If month is mising, this is probematic
# if is missing, assume it was collected mid-month, and assing it to the 15th day

	data$day[data$day == 99] <- 15
	data$day[data$day == 0] <- 15
	data$day[is.na(data$day)] <- 15

# as year is irrelvent, lets set all values to a common one
	data$year<-"2000"

# create a new column for a concatenated date and return as Julian day
# if month is missing, NA will be returned

	data$date <- do.call(paste, list(data$month, data$day, data$year))
	data$date <- as.Date(data$date, format=c("%m %d %Y"))
	data$julian.date <- as.numeric(format(data$date, "%j"))
	data <- na.omit(data)

# because this dataset lacks data about reproductive state, lets remove specimens collected at outlying months
# we assume that specimens with at least 1 flower are most commonly collected and that fruiting specimens are rarely collected.
	counts <- ddply(data, .(data$species, data$month), nrow)
	colnames(counts)=c("species","month","rec_month")
	counts2 <- count(data$species) ; colnames(counts2)=c("species","total_rec")
	coll_sp <- merge(counts,counts2,by="species",all=T)
# count percent collected per month/species
	coll_sp$percent_per_month<-coll_sp$rec_month/coll_sp$total_rec
	data2 <- merge(data,coll_sp,by=c("species","month"),all=T)
# remove months for which >= 10% of all records/species are from
	data2 <- subset(data2, percent_per_month >= 0.1)
	data2<-data2[,c(1,6)]; names(data2)=c("species","time")

# convert Julian day to radians	to use the function
	
	timeRad <- (data2$time)/366 * 2 * pi

# get all combinations of species and caluclate phenology overlap	
	combos_time <- combn(unique(data2$species),2, simplify=T)
	time_overlap<-list(); time_i<-list(); time_j<-list()

	for(k in seq_along(combos_time[1,])){
		i <- combos_time[1,k]; time_i[k]<- i
		j <- combos_time[2,k]; time_j[k]<- j
	
# extract radian values for ith species 		
		time_int_1 <- timeRad[data2 == i]
  
# if 1 or 2 data points, generate a +/-10 day (in radians) buffer around the mean radian   
		ifelse(length(time_int_1) <= 2,time_int_1<-c(time_int_1,mean(time_int_1) + 0.1716717,mean(time_int_1) - 0.1716717),time_int_1)
  
# same for jth  
  		time_int_2 <- timeRad[data2 == j]
  		ifelse(length(time_int_2) <= 2,time_int_2<-c(time_int_2,mean(time_int_2) + 0.1716717,mean(time_int_2) - 0.1716717),time_int_2)
  	
# calcualte overlap coefficent for i,j 		
  		time_est <- overlapEst(time_int_1, time_int_2, type="Dhat1",kmax = 6, adjust=3)
  		time_overlap[k] <- time_est
  		}
  		
pheno_o <- data.frame(matrix(unlist(time_overlap)))
time_i <- data.frame(matrix(unlist(time_i)))
time_j <- data.frame(matrix(unlist(time_j)))

pheno_o_df <- cbind.data.frame(time_i,time_j,pheno_o)
colnames(pheno_o_df)=c("sp1","sp2","time_overlap")
pheno_o_df2<-pheno_o_df[,c(2,1,3)];colnames(pheno_o_df2)=c("sp1","sp2","time_overlap")
pheno_o_df<-rbind.data.frame(pheno_o_df,pheno_o_df2)
return(pheno_o_df)
}