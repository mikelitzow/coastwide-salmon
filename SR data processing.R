setwd("/Users/MikeLitzow/Documents/R/NSF-GOA")
require(gplots)
require(nlme)
require(scales) 
require(marmap)
require(dplyr)
require(zoo)
library(maps)
library(mapdata)
library(ncdf4)
library(chron)
library(fields)

dat <- read.csv("data.update.Nov.2018.csv")
head(dat)
pink <- filter(dat, species=="Pink")

# this is a bit verbose, but selects keeper runs!
start <- tapply(pink$broodYr, pink$stock, min, na.rm=T)

end <- tapply(pink$broodYr, pink$stock, max, na.rm=T)
identical(names(start), names(end)) #T
pink.runs <- data.frame(run=as.character(names(end)), start=start, end=end)
pink.runs <- na.omit(pink.runs) # there isn't any improvement over what we used before

##########
# check chum data!
chum <- filter(dat, species=="Chum")

# this is a bit verbose, but selects keeper runs!
start <- tapply(chum$broodYr, chum$stock, min, na.rm=T)

end <- tapply(chum$broodYr, chum$stock, max, na.rm=T)
identical(names(start), names(end)) #T
chum.runs <- data.frame(run=as.character(names(end)), start=start, end=end)
chum.runs <- na.omit(chum.runs)
chum.runs # also no improvement here...
#########################
# check sockeye
sock <- filter(dat, species=="Sockeye")

# this is a bit verbose, but selects keeper runs!
start <- tapply(sock$broodYr, sock$stock, min, na.rm=T)

end <- tapply(sock$broodYr, sock$stock, max, na.rm=T)
identical(names(start), names(end)) #T
sock.runs <- data.frame(run=as.character(names(end)), start=start, end=end)
sock.runs <- na.omit(sock.runs)

# data must at least run between 1972 and 2007
sockeye.use <- sock.runs[sock.runs$start <= 1972 & sock.runs$end >= 2007,]
sockeye.use # this is much better than the last version wer had!

sock <- sock[sock$stock %in% sock.use$run,]

# and... drop the runs that Rich Brenner identified as hatchery subsidized,
# as well as the combined Chignik
unique(sock$stock)
drop <- c("Coghill", "Eshamy", "Copper")
sockeye.use <- sockeye.use[!sockeye.use$run %in% drop,]
sockeye <- sock[sock$stock %in% sockeye.use$run,]

# rename some columns
colnames(sockeye)[7:9] <- c("brood.yr", "spawners", "recruits")
# add entry yr and ln.rs
sockeye$entry.yr <- sockeye$brood.yr+2 # this is a nominal entry age!
# Mueter et al. 2002 say all southern stocks enter at age 2; AK at 2-3
sockeye$ln.rs <- log(sockeye$recruits/sockeye$spawners)

head(sockeye)

# need lat/long!
sockeye$lat <- sockeye$long <- NA

info <- read.csv("/Users/MikeLitzow 1/Documents/R/NSF-GOA/sockeye.lat.long.csv")

colnames(info)
# change a name
info$stock <- as.character(info$stock)

change <- grep("Black", info$stock)
info$stock[change] <- "Black Lake"
change <- grep("Chignik", info$stock)
info$stock[change] <- "Chignik Lake"
change <- grep("E.Upper Station", info$stock)
info$stock[change] <- "E. Upper Station"
change <- grep("L.Upper Station", info$stock)
info$stock[change] <- "L. Upper Station"
change <- grep("Lake Washington", info$stock)
info$stock[change] <- "Washington"

colnames(sockeye)

# and! limit to 1960-present following GOA analysis
sockeye <- sockeye[sockeye$brood.yr >= 1960,]

# and limit to 2010
sockeye <- sockeye[sockeye$brood.yr <= 2010,]

# I know there must be a better way to do this!

for(i in 1:nrow(sockeye)){
  # i <- 1
  index <- as.character(sockeye$stock[i]) 
  sockeye$lat[i] <- info$lat[as.character(info$stock)==index]
  sockeye$long[i] <- info$long[as.character(info$stock)==index]
}

# now set for all the Fraser runs
fill <- grep("Fraser", sockeye$sub.region)
sockeye$lat[fill] <- 49.12
sockeye$long[fill] <- 123.06

# check
check <- is.na(sockeye$lat)
sum(check) # looks good!

unique(sockeye$stock)

par(las=1)
plot(360-tapply(sockeye$long, sockeye$stock, mean), tapply(sockeye$lat, sockeye$stock, mean), 
     ylab="Latitude", xlab="Longitude", xaxt="n", pch=21, bg="red", ylim=c(46,68), 
     xlim=c(180,240))
axis(1, at=seq(180,240,10), labels=seq(180, 120,-10))
map('world2Hires',fill=F,xlim=c(132,250), ylim=c(20,68),add=T, lwd=1)
mtext("sockeye salmon locations")

#########################################
# now get sst data!
# make column for local sst at ocean entry
sockeye$loc.sst <- NA


# now load the file
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/NSF-GOA/sst.mnmean.v4.nc")

# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates

# save months and years for use later on
m <- months(d)
yrs <- years(d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

####
# begin with EBS
# EBS
# 54-62 deg. N, 186-202 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=9)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,14,1153), count=c(9,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

EBSx <- x
EBSy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out as needed
blank <- c("N54E186", "N54E188", "N54E190",  "N54E196", "N54E198", "N54E200", "N54E202",
           "N56E186", "N56E188", "N58E186",  "N56E200", "N56E202")
SST[,blank] <- NA

# plot mean temperature pattern to check
EBS.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(50,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!


plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

###########
# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# # define seasons for Alaska
# win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# # and define winter years
# win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# # now we need to assign Nov and Dec data to the year corresponding to Jan
# win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# # calculate winter mean SST
# SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
# use <- m %in% win # logical vector that's true for winter months
# SST.mean <- SST.mean[use] # select winter means only
# win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
# SST.win <- tapply(SST.mean, win.yrs, mean)
# f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
# win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
# SST.win <- SST.win[win.count]

# check with a plot
# plot(names(SST.win), SST.win, type="b")
# looks right!

# ok, now need to define larger region!
sockeye$our.region <- "GOA"

# identify ebs runs
stocks <- unique(sockeye$stock)

sockeye$our.region[sockeye$stock %in% stocks[1:9]] <- "EBS"

sockeye$our.region[sockeye$lat<52] <- "South"

# check
tapply(sockeye$lat, list(sockeye$stock, sockeye$our.region), mean)

# # now smooth with 2-yr running mean
# win.3 <- rollmean(SST.win, 3, fill=NA)
# smooth <- rollmean(SST.win, 2, align="left", fill=NA)
# 
# head(sockeye)
# 
# # now set smooth sst to the winter before and the winter after ocean entry
# for(i in 1:nrow(sockeye)){
#   #  i <- 1
#   if(sockeye$our.region[i] == "EBS") sockeye$w.sst.3[i] <- 
#       win.3[names(win.3)==sockeye$entry.yr[i]] else sockeye$w.sst.3[i] <- sockeye$w.sst.3[i]
#       
#       if(sockeye$our.region[i] == "EBS") sockeye$w.sst.2[i] <- 
#           smooth[names(smooth)==sockeye$entry.yr[i]] else sockeye$w.sst.2[i] <- sockeye$w.sst.2[i]
#           
#           if(sockeye$our.region[i] == "EBS") sockeye$w.sst.1[i] <- 
#               SST.win[names(SST.win)==sockeye$entry.yr[i]] else sockeye$w.sst.1[i] <- sockeye$w.sst.1[i]
# }

############################################################
# step back - look at entire coastal range for sst values!
# 46-68 deg. N, 186-232 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=26)
y <- ncvar_get(nc, "lat", start=11, count=12)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,11,1153), count=c(26,12,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,12:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now get anomalies
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std


# choose coastal stations

coast <- c("N68E196", "N68E194", "N66E194", "N66E192", "N64E200", "N64E198", "N64E196",
           "N64E194", "N62E194", "N60E194", "N60E196", "N60E198", "N58E202", "N58E200",
           "N58E198", "N56E198", "N56E196",
           "N52E228", "N52E230", "N52E232", "N54E196", 
           "N54E198", "N54E200","N54E226", "N54E228", "N54E230", "N56E202", "N56E204", 
           "N56E206", "N56E224", "N56E226", "N56E228", "N58E204", "N58E206", "N58E208", 
           "N58E222", "N58E224", "N58E226", "N60E208", "N60E210", "N60E212", "N60E214", 
           "N60E216", "N60E218", "N60E220", "N50E232", "N50E234", "N50E236" ,"N48E236",
           "N48E234", "N46E236")

coast.sst <- SST
use <- colnames(coast.sst) %in% coast
coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# let's get an object with the lat and long for each cell!


temp1 <- strsplit(coast,"E")
temp2 <- matrix(unlist(temp1), ncol=2, byrow=TRUE)
temp3 <- strsplit(temp2[,1], "N")
temp4 <- matrix(unlist(temp3), ncol=2, byrow=TRUE)

coast.cells <- data.frame(cell=coast, lat=as.numeric(temp4[,2]), long=as.numeric(temp2[,2]))
str(coast.cells)

coast.cells <- coast.cells[order(coast.cells$lat, decreasing = T),]

# add EBS/GOA factor!
coast.cells$region <- "GOA"
change <- coast.cells$lat >= 57 & coast.cells$long <=202
coast.cells$region[change] <- "EBS"

change <- coast.cells$lat == 56 & coast.cells$long <=200
coast.cells$region[change] <- "EBS"

# and check!
EBS.coast.cells <- coast.cells$cell[coast.cells$region=="EBS"]

EBS.coast.sst <- SST

use <- colnames(EBS.coast.sst) %in% EBS.coast.cells
EBS.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(EBS.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

GOA.coast.cells <- coast.cells$cell[coast.cells$region=="GOA"]

GOA.coast.sst <- SST

use <- colnames(GOA.coast.sst) %in% GOA.coast.cells
GOA.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(GOA.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now need an object with the lat longs of each sockeye run!
p.lat <-tapply(sockeye$lat, sockeye$stock, mean, na.rm=T)
keep <- !is.na(p.lat)
p.lat <- p.lat[keep]

p.long <-tapply(sockeye$long, sockeye$stock, mean, na.rm=T)
keep <- !is.na(p.long)
p.long <- p.long[keep]

stocks <- unique(sockeye$stock)
stocks <- stocks[order(stocks)]

sockeye.locations <- data.frame(stock=stocks, lat=p.lat, long=p.long, region=NA, coast.cells=NA)
str(sockeye.locations)
rownames(sockeye.locations) <- 1:nrow(sockeye.locations)

ebs <- c(1,10,14,17,21,22,29,30,33)

sockeye.locations$region[-ebs] <- "GOA"
sockeye.locations$region[ebs] <- "EBS"
# looks good...

# change coast.cells long back to degrees W!
coast.cells$Wlong <- 180-(coast.cells$long-180)

# get cells <= 500 km away for ebs runs...
ebs.x1 <- as.matrix(cbind(sockeye.locations$long[sockeye.locations$region=="EBS"], 
                          sockeye.locations$lat[sockeye.locations$region=="EBS"]))
rownames(ebs.x1) <- sockeye.locations$stock[sockeye.locations$region=="EBS"]
ebs.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="EBS"], 
                          coast.cells$lat[coast.cells$region=="EBS"]))


dist.ebs <- rdist.earth(ebs.x1, ebs.x2, miles=F)

colnames(dist.ebs) <- coast.cells$cell[coast.cells$region=="EBS"]


for(i in 1:nrow(dist.ebs)){
 #  i <- 1
  run <- rownames(dist.ebs)[i]
  use <- dist.ebs[i,] <= 500
  cells <- as.vector(colnames(dist.ebs)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  sockeye.locations$coast.cells[sockeye.locations$stock==run] <- paste(cells, collapse=",")
}

###
# get cells <= 500 km away for goa runs...
goa.x1 <- as.matrix(cbind(sockeye.locations$long[sockeye.locations$region=="GOA"], 
                          sockeye.locations$lat[sockeye.locations$region=="GOA"]))
rownames(goa.x1) <- sockeye.locations$stock[sockeye.locations$region=="GOA"]
goa.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="GOA"], 
                          coast.cells$lat[coast.cells$region=="GOA"]))


dist.goa <- rdist.earth(goa.x1, goa.x2, miles=F)

colnames(dist.goa) <- coast.cells$cell[coast.cells$region=="GOA"]

for(i in 1:nrow(dist.goa)){
  # i <- 1
  run <- rownames(dist.goa)[i]
  use <- dist.goa[i,] <= 500
  cells <- as.vector(colnames(dist.goa)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  sockeye.locations$coast.cells[sockeye.locations$stock==run] <- paste(cells, collapse=",")
}

# now, some random-ish checks!
pdf("sockeye runs and local sst cells.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

for(i in 1:nrow(sockeye.locations)){
 # i <- 3
  sub <- sockeye.locations[i,]
  sub.sst <- SST
  temp <- strsplit(sub$coast.cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  lim <- range(colMeans(SST), na.rm = T)
  
  sub.mean <- colMeans(sub.sst)
  z <- t(matrix(sub.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=tim.colors(64), zlim=lim, xlim=c(160,240), ylim=c(44,70), ylab="", xlab="", yaxt="n", xaxt="n")
  #contour(x, y, z, add=T)  
  map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
  points(360-sub$long, sub$lat, pch=21, bg="red", cex=1.5)
  mtext(sub$stock)
}
dev.off()


# looks good!
# actually....looks GREAT!

# now, we need to get seasonal means based on the periods identified by Mueter et al. 2005!
# these are years relative to brood year!!
AK.mo.l1 <- c("Jan", "Feb")
AK.mo.l2 <- c("Jan", "Feb", "Mar", "Apr", "May")
AK.mo.l3 <- c("Jan", "Feb", "Mar", "Apr", "May", "Aug", "Sep", "Oct", "Nov", "Dec")


S.mo.l1 <- c("Jan", "Mar", "Apr", "May")
S.mo.l2 <- c("Jan", "Mar", "Apr")

sockeye.locations$region2 <- "AK"
south <- c(3,5,7,9,11,12,13,19,20,23:28,31,32)
sockeye.locations$region2[south] <- "South"

# now we need to go through and select the proper months to make a local sst time series!!
stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
   # i <- 2
  sub <- sockeye.locations[i,]
  
  if(sub$region2=="AK"){
    # get correct cells
    sub.sst <- SST.anom # so these are scaled sst time series! (mean 0, unit variance)
    temp <- strsplit(sub$coast.cells,",")
    use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
    
    keep <- colnames(sub.sst) %in% use
    sub.sst[,!keep] <- NA
    
    # year 1 (wrt b.y.)
    keep <- months(rownames(sub.sst)) %in% AK.mo.l1
    sub.sst1 <- sub.sst
    sub.sst1[!keep,] <- NA 
    
    keep <- months(rownames(sub.sst)) %in% AK.mo.l2
    sub.sst2 <- sub.sst
    sub.sst2[!keep,] <- NA
    
    keep <- months(rownames(sub.sst)) %in% AK.mo.l3
    sub.sst3 <- sub.sst
    sub.sst3[!keep,] <- NA
    
    # now get means!
    years <- years(row.names(sub.sst))
    yr <- data.frame(yr=1960:2012)
    yr$mean.loc <- NA
    for(y in 1960:2012){
      # y <- 1960
      temp <- c(as.vector(sub.sst2[years==y,]), 
                as.vector(sub.sst1[years==(y-1),]), 
                as.vector(sub.sst3[years==(y+1),]), na.rm=T)
      temp <- na.omit(temp)
      yr$mean.loc[yr$yr==y] <- mean(temp)
    }
    
    # now plug into the sockeye data frame!
    sub.sub <- sockeye[sockeye$stock==sub$stock,]
    for(j in 1:nrow(sub.sub)){
      #j <- 1
      sub.sub$loc.sst[j] <- yr$mean.loc[yr$yr==sub.sub$entry.yr[j]] 
      
    }
    
    # and now put the local sst time series into the sockeye data frame
    
    sockeye$loc.sst[sockeye$stock==sub$stock] <- sub.sub$loc.sst
    # temp <- rowMeans(sub.sst, na.rm=T) # monthly mean anomaly for this region!
    # ts <- tapply(temp, yr, mean)
    
  }
  
  # now the southern stocks!!!
  if(sub$region2=="South"){
    # get correct cells
    sub.sst <- SST.anom # so these are scaled sst time series! (mean 0, unit variance)
    temp <- strsplit(sub$coast.cells,",")
    use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
    
    keep <- colnames(sub.sst) %in% use
    sub.sst[,!keep] <- NA
    
    # year 1 (wrt b.y.)
    keep <- months(rownames(sub.sst)) %in% S.mo.l1
    sub.sst1 <- sub.sst
    sub.sst1[!keep,] <- NA 
    
    keep <- months(rownames(sub.sst)) %in% S.mo.l2
    sub.sst2 <- sub.sst
    sub.sst2[!keep,] <- NA
    
    # now get means!
    years <- years(row.names(sub.sst))
    yr <- data.frame(yr=1960:2012)
    yr$mean.loc <- NA
    for(y in 1960:2012){
      # y <- 1960
      temp <- c(as.vector(sub.sst2[years==y,]), 
                as.vector(sub.sst1[years==(y-1),]), na.rm=T)
      temp <- na.omit(temp)
      yr$mean.loc[yr$yr==y] <- mean(temp)
    }
    
    # now plug into the sockeye data frame!
    sub.sub <- sockeye[sockeye$stock==sub$stock,]
    for(j in 1:nrow(sub.sub)){
      #j <- 1
      sub.sub$loc.sst[j] <- yr$mean.loc[yr$yr==sub.sub$entry.yr[j]] 
      
    }
    
    # and now put the local sst time series into the sockeye data frame
    
    sockeye$loc.sst[sockeye$stock==sub$stock] <- sub.sub$loc.sst
    # temp <- rowMeans(sub.sst, na.rm=T) # monthly mean anomaly for this region!
    # ts <- tapply(temp, yr, mean)
    
  }
}

# check time series!
pdf("sockeye local sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
#   i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$loc.sst, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()
# these look right to me!

# save sockeye data
write.csv(sockeye,"updated.sockeye.data.csv")

######################################

# # now plug in winter goa temps!!
# 
# # Extract GOA SST, 54-61 deg. N, 200-226 deg. E:
# # now load the file
# nc <- nc_open("sst.mnmean.v4.nc")
# # get info
# nc
# 
# # view dates (middle of month):
# ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)
# 
# # assign dates
# d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))
# 
# # we will start with Jan 1950 and end with most recent month
# d[c(1153,length(d))] # check - those are correct
# d <- d[1153:length(d)] # restrict to the desired dates
# 
# 
# ncvar_get(nc, "lon")     # view longitudes (degrees East)
# ncvar_get(nc, "lat")     # view latitudes
# 
# x <- ncvar_get(nc, "lon", start=99, count=16)
# y <- ncvar_get(nc, "lat", start=14, count=5)
# x; y # check
# 
# # get required sst data
# SST <- ncvar_get(nc, "sst", start=c(99,14,1141), count=c(16,5,length(d)), verbose = T)
# 
# # check
# dim(SST) # 14 longitudes, 5 latitudes, 798 months
# 
# # need to change to matrix for easier use
# SST <- aperm(SST, 3:1) # transpose array
# SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
# y <- rev(y)  # Also reverse corresponding vector of lattidues
# SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
# dim(SST)  # Matrix with column for each grid point, rows for monthly means
# 
# GOAx <- x
# GOAy <- y
# 
# # Keep track of corresponding latitudes and longitudes of each column:
# lat <- rep(y, length(x))   # Vector of latitudes
# lon <- rep(x, each = length(y))   # Vector of longitudes
# 
# # set names
# dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# 
# #blank out Bristol Bay
# blank <- c("N56E200", "N58E200", "N58E202", "N56E196", "N56E198", "N58E196",
#            "N58E198", "N60E196", "N60E198")
# SST[,blank] <- NA
# 
# # plot mean temperature pattern to check
# GOA.mean <- SST.mean <- colMeans(SST)
# z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(x,y,z, col=tim.colors(64), xlim=c(190,230), ylim=c(48,64))
# contour(x, y, z, add=T)  
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# # looks good!
# 
# 
# 
# # and another check!
# 
# plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")
# 
# yrs <- years(d)
# ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)
# 
# plot(1950:2016, ann.mean[1:67], type="l")
# 
# # correcto!
# 
# # now remove seasonal signal and scale 
# mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
# mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 
# 
# #now need to account for trailing months (i.e., fractions of a year that are available with updated data)
# 
# add <- length(d)-12*floor(length(d)/12)
# 
# # add in additional required months
# mu <- rbind(mu, mu[1:add,])
# # check
# identical(nrow(mu), nrow(SST)) #true!
# 
# std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# # and stack as for mu
# std <- std[rep(1:12, floor(length(d)/12)),] 
# # add in additional required months
# std <- rbind(std, std[1:add,])
# 
# # now calculate anomalies...
# SST.anom <- (SST - mu)/std
# 
# # define seasons for Alaska
# win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# # and define winter years
# win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# # now we need to assign Nov and Dec data to the year corresponding to Jan
# win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1
# 
# # calculate winter mean SST
# SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
# use <- m %in% win # logical vector that's true for winter months
# SST.mean <- SST.mean[use] # select winter means only
# win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
# SST.win <- tapply(SST.mean, win.yrs, mean)
# f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
# win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
# SST.win <- SST.win[win.count]
# 
# # check with a plot
# plot(names(SST.win), SST.win, type="b")
# # looks right!
# 
# 
# # now smooth with 2-yr  and 3-yr running mean
# smooth <- rollmean(SST.win, 2, align="left", fill=NA)
# win.3 <- rollmean(SST.win, 3, fill=NA)
# 
# head(sockeye)
# 
# # now set smooth sst to the winter before and the winter after ocean entry
# for(i in 1:nrow(sockeye)){
#   #  i <- 1
#   
#   if(sockeye$our.region[i] == "GOA") sockeye$w.sst.3[i] <- 
#       win.3[names(win.3)==sockeye$entry.yr[i]] else sockeye$w.sst.3[i] <- sockeye$w.sst.3[i]
#       
#       if(sockeye$our.region[i] == "GOA") sockeye$w.sst.2[i] <- 
#           smooth[names(smooth)==sockeye$entry.yr[i]] else sockeye$w.sst.2[i] <- sockeye$w.sst.2[i]
#           
#           if(sockeye$our.region[i] == "GOA") sockeye$w.sst.1[i] <- 
#               SST.win[names(SST.win)==sockeye$entry.yr[i]] else sockeye$w.sst.1[i] <- sockeye$w.sst.1[i]
# }
# 
# ###############################################
# # next, get winter sst for southern stocks!
# 
# ncvar_get(nc, "lon")     # view longitudes (degrees East)
# ncvar_get(nc, "lat")     # view latitudes
# 
# # 46-68 deg. N, 186-232 deg. E
# x <- ncvar_get(nc, "lon", start=115, count=7)
# y <- ncvar_get(nc, "lat", start=19, count=4)
# x; y # check
# 
# # get required sst data
# SST <- ncvar_get(nc, "sst", start=c(115,19,1141), count=c(7,4,length(d)), verbose = T)
# 
# # check
# dim(SST) # 14 longitudes, 5 latitudes, 798 months
# 
# # need to change to matrix for easier use
# SST <- aperm(SST, 3:1) # transpose array
# SST <- SST[,4:1,]  # Reverse order of latitudes to be increasing for plotting
# y <- rev(y)  # Also reverse corresponding vector of lattidues
# SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
# dim(SST)  # Matrix with column for each grid point, rows for monthly means
# 
# Southx <- x
# Southy <- y
# 
# # Keep track of corresponding latitudes and longitudes of each column:
# lat <- rep(y, length(x))   # Vector of latitudes
# lon <- rep(x, each = length(y))   # Vector of longitudes
# 
# # set names
# dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# 
# blank <- c("N50E228", "N48E228", "N46E228", "N48E230", "N46E230")
# 
# SST[,blank] <- NA
# 
# # plot mean temperature pattern to check
# South.mean <- SST.mean <- colMeans(SST)
# z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
# image(x,y,z, col=tim.colors(64), xlim=c(190,238), ylim=c(42,64))
# 
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# # looks good!
# 
# # and another check!
# 
# plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")
# 
# yrs <- years(d)
# ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)
# 
# plot(1950:2016, ann.mean[1:67], type="l")
# 
# # correcto!
# 
# # now remove seasonal signal and scale 
# mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
# mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 
# 
# #now need to account for trailing months (i.e., fractions of a year that are available with updated data)
# 
# add <- length(d)-12*floor(length(d)/12)
# 
# # add in additional required months
# mu <- rbind(mu, mu[1:add,])
# # check
# identical(nrow(mu), nrow(SST)) #true!
# 
# std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# # and stack as for mu
# std <- std[rep(1:12, floor(length(d)/12)),] 
# # add in additional required months
# std <- rbind(std, std[1:add,])
# 
# # now calculate anomalies...
# SST.anom <- (SST - mu)/std
# 
# # define seasons
# win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# # and define winter years
# win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# # now we need to assign Nov and Dec data to the year corresponding to Jan
# win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1
# 
# # calculate winter mean SST
# SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
# use <- m %in% win # logical vector that's true for winter months
# SST.mean <- SST.mean[use] # select winter means only
# win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
# SST.win <- tapply(SST.mean, win.yrs, mean)
# f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
# win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
# SST.win <- SST.win[win.count]
# 
# # check with a plot
# plot(names(SST.win), SST.win, type="b")
# # looks right!
# 
# 
# # now smooth with 2-yr and 3-yr running mean
# smooth <- rollmean(SST.win, 2, align="left", fill=NA)
# win.3 <- rollmean(SST.win, 3, fill=NA)
# head(sockeye)
# 
# # now set smooth sst to the winter before and the winter after ocean entry
# for(i in 1:nrow(sockeye)){
#   #  i <- 1
#   
#   if(sockeye$our.region[i] == "South") sockeye$w.sst.3[i] <- 
#       win.3[names(win.3)==sockeye$entry.yr[i]] else sockeye$w.sst.3[i] <- sockeye$w.sst.3[i]
#       
#       if(sockeye$our.region[i] == "South") sockeye$w.sst.2[i] <- 
#           smooth[names(smooth)==sockeye$entry.yr[i]] else sockeye$w.sst.2[i] <- sockeye$w.sst.2[i]
#           
#           if(sockeye$our.region[i] == "South") sockeye$w.sst.1[i] <- 
#               SST.win[names(SST.win)==sockeye$entry.yr[i]] else sockeye$w.sst.1[i] <- sockeye$w.sst.1[i]
# }
# 
# # check!
# 
# pdf("sockeye 3-yr winter sst time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 3
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$w.sst.3, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# 
# pdf("sockeye 2-yr winter sst time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 3
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$w.sst.2, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# pdf("sockeye 1-yr winter sst time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 10
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$w.sst.1, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# # all looks good!
# 
# # finally, add PDO and NPGO
# 
# pdo <- read.csv("pdo.8.25.15.csv")
# head(pdo)
# pdo$sm <- rollapply(pdo$NDJFM, 2, mean, fill=NA, align="left")
# pdo$pdo3 <- rollapply(pdo$NDJFM, 3, mean, fill=NA)
# 
# npgo <- read.csv("npgo8.25.15.csv")
# head(npgo)
# npgo$sm <- rollapply(npgo$NDJFM, 2, mean, fill=NA, align="left")
# npgo$npgo3 <- rollapply(npgo$NDJFM, 3, mean, fill=NA)
# 
# sockeye$pdo <- sockeye$pdo2 <- sockeye$pdo3 <- sockeye$npgo <- sockeye$npgo2 <- sockeye$npgo3 <- NA
# 
# for(i in 1:nrow(sockeye)){
#   # i <- 1
#   sockeye$pdo[i] <- pdo$NDJFM[pdo$YEAR==sockeye$entry.yr[i]]
#   sockeye$pdo2[i] <- pdo$sm[pdo$YEAR==sockeye$entry.yr[i]]
#   sockeye$pdo3[i] <- pdo$pdo3[pdo$YEAR==sockeye$entry.yr[i]]
#   
#   sockeye$npgo[i] <- npgo$NDJFM[npgo$YEAR==sockeye$entry.yr[i]]
#   sockeye$npgo2[i] <- npgo$sm[npgo$YEAR==sockeye$entry.yr[i]]
#   sockeye$npgo3[i] <- npgo$npgo3[npgo$YEAR==sockeye$entry.yr[i]]
# }
# 
# # final checks!
# 
# 
# pdf("sockeye 3-yr winter pdo time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 3
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$pdo3, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# 
# pdf("sockeye 2-yr winter pdo time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 3
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$pdo2, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# pdf("sockeye 1-yr winter pdo time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 10
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$pdo, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# # final checks!
# pdf("sockeye 3-yr winter npgo time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))
# 
# stocks <- unique(sockeye$stock)
# 
# for(i in 1:length(stocks)){
#   # i <- 3
#   sub <- sockeye[sockeye$stock==stocks[i],]
#   plot(sub$entry.yr, sub$npgo3, type="b", ylab="", xlab="")
#   
#   mtext(stocks[i])
# }
# dev.off()
# 
# pdf("sockeye 2-yr winter npgo time series.pdf", 10,8)
# par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("sockeye 1-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# export!
write.csv(sockeye, "coastwide sockeye data UPDATED 6-20-17.csv")

# plot all 3 winter areas!
pdf("winter sst areas.pdf", 5, 10)
par(mfrow=c(3,1), mar=c(2,2,2,0.5), oma=c(0,0,2,0))

lim=range(EBS.mean, GOA.mean, South.mean, na.rm=T)
# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)

z <- t(matrix(EBS.mean,length(EBSy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(EBSx,EBSy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("EBS")

sub <- sockeye.locations[sockeye.locations$region=="EBS",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

#
z <- t(matrix(GOA.mean,length(GOAy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(GOAx,GOAy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("GOA")

sub <- sockeye.locations[sockeye.locations$region=="GOA" & sockeye.locations$region2=="AK",]
unique(sub$stock)
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")
##NB - ERRORS HERE WITH REGIONAL DEFINITIONS, NEEDS TO BE FIXED! (6/20/17)
#
z <- t(matrix(South.mean,length(Southy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(Southx,Southy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("South")

sub <- sockeye.locations[sockeye.locations$region2=="South",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

mtext("Winter SST areas", outer=T)
dev.off()

# need lat/long!
pink$lat <- pink$long <- NA

setwd("/Users/MikeLitzow/Documents/R/NSF-GOA/MichaelMalick-pink-chum-database-af9da57")
info <- read.csv("pink_info.csv")

setwd("/Users/MikeLitzow/Documents/R/NSF-GOA")

colnames(info)
colnames(pink)
lats <- tapply(info$lat, info$stock.id, mean)
lons <- tapply(info$lon, info$stock.id, mean)

# I know there must be a better way to do this!

for(i in 1:nrow(pink)){
  #  i <- 1
  index <- pink$stock.id[i] 
  pink$lat[i] <- lats[names(lats)==index]
  pink$long[i] <- lons[names(lons)==index]
}

unique(pink$stock)

par(las=1)
plot(360-tapply(pink$long, pink$full.stock, mean), tapply(pink$lat, pink$full.stock, mean), 
     ylab="Latitude", xlab="Longitude", xaxt="n", pch=21, bg="red", ylim=c(46,68), 
     xlim=c(180,240))
axis(1, at=seq(180,240,10), labels=seq(180, 120,-10))
map('world2Hires',fill=F,xlim=c(132,250), ylim=c(20,68),add=T, lwd=1)
mtext("Pink salmon locations")

#########################################
# now get sst data!
# make columns for 2-yr winter sst, 1-yr winter sst, and local sst for ocean entry
pink$w.sst.2 <- pink$w.sst.1 <- pink$loc.sst <- NA


# now load the file
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/sst.mnmean.v4.nc")

# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates

# save months and years for use later on
m <- months(d)
yrs <- years(d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

####
# begin with EBS
# EBS
# 54-68 deg. N, 186-202 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=9)
y <- ncvar_get(nc, "lat", start=11, count=8)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,11,1153), count=c(9,8,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,8:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

EBSx <- x
EBSy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out as needed
blank <- c("N54E186", "N54E188", "N54E190",  "N54E196", "N54E198", "N54E200", "N54E202",
           "N56E186", "N56E188", "N58E186",  "N56E200", "N56E202", "N68E186", "N68E188", 
           "N68E190")
SST[,blank] <- NA

# plot mean temperature pattern to check
EBS.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(50,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!


plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

###########
# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!

# ok, now need to define larger region!
pink$our.region <- "GOA"

# identify ebs runs
stocks <- unique(pink$full.stock)

pink$our.region[pink$full.stock %in% stocks[c(32:37)]] <- "EBS"

pink$our.region[pink$lat<52] <- "South"

# now smooth with 2-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(pink)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(pink)){
  #  i <- 1
  if(pink$our.region[i] == "EBS") pink$w.sst.2[i] <- 
      smooth[names(smooth)==pink$entry.yr[i]] else pink$w.sst.2[i] <- pink$w.sst.2[i]
      
      if(pink$our.region[i] == "EBS") pink$w.sst.1[i] <- 
          SST.win[names(SST.win)==pink$entry.yr[i]] else pink$w.sst.1[i] <- pink$w.sst.1[i]
}

############################################################
# step back - look at entire coastal range for sst values!
# 46-68 deg. N, 186-232 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=26)
y <- ncvar_get(nc, "lat", start=11, count=12)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,11,1153), count=c(26,12,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,12:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now get anomalies
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std



# choose coastal stations

coast <- c("N68E196", "N68E194", "N66E194", "N66E192", "N64E200", "N64E198", "N64E196",
           "N64E194", "N62E194", "N60E194", "N60E196", "N60E198", "N58E202", "N58E200",
           "N58E198", "N56E198", "N56E196",
           "N52E228", "N52E230", "N52E232", "N54E196", 
           "N54E198", "N54E200","N54E226", "N54E228", "N54E230", "N56E202", "N56E204", 
           "N56E206", "N56E224", "N56E226", "N56E228", "N58E204", "N58E206", "N58E208", 
           "N58E222", "N58E224", "N58E226", "N60E208", "N60E210", "N60E212", "N60E214", 
           "N60E216", "N60E218", "N60E220", "N50E232", "N50E234", "N50E236" ,"N48E236",
           "N48E234", "N46E236")

coast.sst <- SST
use <- colnames(coast.sst) %in% coast
coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# let's get an object with the lat and long for each cell!


temp1 <- strsplit(coast,"E")
temp2 <- matrix(unlist(temp1), ncol=2, byrow=TRUE)
temp3 <- strsplit(temp2[,1], "N")
temp4 <- matrix(unlist(temp3), ncol=2, byrow=TRUE)

coast.cells <- data.frame(cell=coast, lat=as.numeric(temp4[,2]), long=as.numeric(temp2[,2]))
str(coast.cells)

coast.cells <- coast.cells[order(coast.cells$lat, decreasing = T),]

# add EBS/GOA factor!
coast.cells$region <- "GOA"
change <- coast.cells$lat >= 57 & coast.cells$long <=202
coast.cells$region[change] <- "EBS"

change <- coast.cells$lat == 56 & coast.cells$long <=200
coast.cells$region[change] <- "EBS"

# and check!
EBS.coast.cells <- coast.cells$cell[coast.cells$region=="EBS"]

EBS.coast.sst <- SST

use <- colnames(EBS.coast.sst) %in% EBS.coast.cells
EBS.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(EBS.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

GOA.coast.cells <- coast.cells$cell[coast.cells$region=="GOA"]

GOA.coast.sst <- SST

use <- colnames(GOA.coast.sst) %in% GOA.coast.cells
GOA.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(GOA.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

###################
# aside - get restricted GOA cells for NSF CNH proposal!
coast.cells$region2 <- coast.cells$region
coast.cells$region2[coast.cells$lat<54] <- NA
coast.cells$region2[coast.cells$lat<56 & coast.cells$long>220] <- NA


proposal.coast.cells <- coast.cells$cell[coast.cells$region2=="GOA"]

proposal.coast.sst <- SST

use <- colnames(proposal.coast.sst) %in% proposal.coast.cells
proposal.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(proposal.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# remove the NA columns to compress!
drop <- apply(proposal.coast.sst, 2, function(x) sum(is.na(x)))

keep <- drop == 0
proposal.coast.sst <- proposal.coast.sst[,keep]
head(proposal.coast.sst)

# export!
write.csv(proposal.coast.sst, "proposal.coast.sst.csv")


####################

# now need an object with the lat longs of each pink run!
p.lat <-tapply(pink$lat, pink$full.stock, mean, na.rm=T)
keep <- !is.na(p.lat)
p.lat <- p.lat[keep]

p.long <-tapply(pink$long, pink$full.stock, mean, na.rm=T)
keep <- !is.na(p.long)
p.long <- p.long[keep]

stocks <- unique(pink$full.stock)
stocks <- stocks[order(stocks)]

pink.locations <- data.frame(stock=stocks, lat=p.lat, long=p.long, region=NA, coast.cells=NA)
str(pink.locations)
rownames(pink.locations) <- 1:nrow(pink.locations)

ebs <- c(11,12,17,18,20,21)

pink.locations$region[-ebs] <- "GOA"
pink.locations$region[ebs] <- "EBS"
# looks good...

# change coast.cells long back to degrees W!
coast.cells$Wlong <- 180-(coast.cells$long-180)

# get cells <= 500 km away for ebs runs...
ebs.x1 <- as.matrix(cbind(pink.locations$long[pink.locations$region=="EBS"], 
                          pink.locations$lat[pink.locations$region=="EBS"]))
ebs.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="EBS"], 
                          coast.cells$lat[coast.cells$region=="EBS"]))


dist.ebs <- rdist.earth(ebs.x1, ebs.x2, miles=F)

colnames(dist.ebs) <- coast.cells$cell[coast.cells$region=="EBS"]

for(i in 1:nrow(dist.ebs)){
  # i <- 1
  run <- rownames(dist.ebs)[i]
  use <- dist.ebs[i,] <= 500
  cells <- as.vector(colnames(dist.ebs)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  pink.locations$coast.cells[pink.locations$stock==run] <- paste(cells, collapse=",")
}

###
# get cells <= 500 km away for goa runs...
goa.x1 <- as.matrix(cbind(pink.locations$long[pink.locations$region=="GOA"], 
                          pink.locations$lat[pink.locations$region=="GOA"]))
goa.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="GOA"], 
                          coast.cells$lat[coast.cells$region=="GOA"]))


dist.goa <- rdist.earth(goa.x1, goa.x2, miles=F)

colnames(dist.goa) <- coast.cells$cell[coast.cells$region=="GOA"]

for(i in 1:nrow(dist.goa)){
  # i <- 1
  run <- rownames(dist.goa)[i]
  use <- dist.goa[i,] <= 500
  cells <- as.vector(colnames(dist.goa)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  pink.locations$coast.cells[pink.locations$stock==run] <- paste(cells, collapse=",")
}

# now, some random-ish checks!
pdf("pink runs and local sst cells.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

for(i in 1:nrow(pink.locations)){
  #i <- 30
  sub <- pink.locations[i,]
  sub.sst <- SST
  temp <- strsplit(sub$coast.cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  lim <- range(colMeans(SST), na.rm = T)
  
  sub.mean <- colMeans(sub.sst)
  z <- t(matrix(sub.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=tim.colors(64), zlim=lim, xlim=c(160,240), ylim=c(44,70), ylab="", xlab="", yaxt="n", xaxt="n")
  #contour(x, y, z, add=T)  
  map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
  points(360-sub$long, sub$lat, pch=21, bg="red", cex=1.5)
  mtext(sub$stock)
}
dev.off()


# looks good!
# actually....looks GREAT!

# now, we need to get seasonal means based on the periods identified by Mueter et al. 2005!
# these are for entry year, i.e., the year after brood year!!
AK.mo <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")
S.mo <- c("Mar", "Apr", "May", "Jun", "Jul") # these aren't significant in Mueter et al 2005! 
# but are the strongest correlations, consistent with the AK results (i.e., around ocean entry)

pink.locations$region2 <- "AK"
south <- c(1,2,7,10,19,24,31:33)
pink.locations$region2[south] <- "South"

# now we need to go through and select the proper months to make a local sst time series!!
stocks

for(i in 1:length(stocks)){
  # i <- 1
  sub <- pink.locations[i,]
  sub.sst <- SST.anom # so these are scaled sst time series! (mean 0, unit variance)
  temp <- strsplit(sub$coast.cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  if(sub$region2=="AK") m.use <- AK.mo else m.use <- S.mo
  keep <- months(rownames(sub.sst)) %in% m.use
  sub.sst <- sub.sst[keep,] 
  yr <- years(rownames(sub.sst))
  temp <- rowMeans(sub.sst, na.rm=T) # monthly mean anomaly for this region!
  ts <- tapply(temp, yr, mean)
  
  # now plug into the pink data frame!
  sub.sub <- pink[pink$full.stock==sub$stock,]
  for(j in 1:nrow(sub.sub)){
    #j <- 1
    sub.sub$loc.sst[j] <- ts[names(ts)==sub.sub$entry.yr[j]] 
    
  }
  
  # and now put the local sst time series into the pink data frame
  
  pink$loc.sst[pink$full.stock==sub$stock] <- sub.sub$loc.sst
}

# check time series!
pdf("pink local sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$loc.sst, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()
# these look right to me!

######################################

# now plug in winter goa temps!!

# Extract GOA SST, 54-61 deg. N, 200-226 deg. E:
# now load the file
nc <- nc_open("sst.mnmean.v4.nc")
# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates


ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=99, count=16)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(99,14,1141), count=c(16,5,length(d)), verbose = T)

# check
dim(SST) # 14 longitudes, 5 latitudes, 798 months

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

GOAx <- x
GOAy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out Bristol Bay
blank <- c("N56E200", "N58E200", "N58E202", "N56E196", "N56E198", "N58E196",
           "N58E198", "N60E196", "N60E198")
SST[,blank] <- NA

# plot mean temperature pattern to check
GOA.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,230), ylim=c(48,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!



# and another check!

plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")

yrs <- years(d)
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

# correcto!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!


# now smooth with 2-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(pink)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(pink)){
  #  i <- 1
  if(pink$our.region[i] == "GOA") pink$w.sst.2[i] <- 
      smooth[names(smooth)==pink$entry.yr[i]] else pink$w.sst.2[i] <- pink$w.sst.2[i]
      
      if(pink$our.region[i] == "GOA") pink$w.sst.1[i] <- 
          SST.win[names(SST.win)==pink$entry.yr[i]] else pink$w.sst.1[i] <- pink$w.sst.1[i]
}

###############################################
# next, get winter sst for southern stocks!

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# 46-68 deg. N, 186-232 deg. E
x <- ncvar_get(nc, "lon", start=115, count=7)
y <- ncvar_get(nc, "lat", start=19, count=4)
x; y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(115,19,1141), count=c(7,4,length(d)), verbose = T)

# check
dim(SST) # 14 longitudes, 5 latitudes, 798 months

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,4:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

Southx <- x
Southy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

blank <- c("N50E228", "N48E228", "N46E228", "N48E230", "N46E230")

SST[,blank] <- NA

# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,238), ylim=c(42,64))

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# and another check!

plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")

yrs <- years(d)
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

# correcto!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!


# now smooth with 2-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(pink)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(pink)){
  #  i <- 1
  if(pink$our.region[i] == "South") pink$w.sst.2[i] <- 
      smooth[names(smooth)==pink$entry.yr[i]] else pink$w.sst.2[i] <- pink$w.sst.2[i]
      
      if(pink$our.region[i] == "South") pink$w.sst.1[i] <- 
          SST.win[names(SST.win)==pink$entry.yr[i]] else pink$w.sst.1[i] <- pink$w.sst.1[i]
}

# check!

pdf("pink 2-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("pink 1-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.1, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# all looks good!

# finally, add PDO and NPGO

pdo <- read.csv("pdo.8.25.15.csv")
head(pdo)
pdo$sm <- rollapply(pdo$NDJFM, 2, mean, fill=NA, align="left")

npgo <- read.csv("npgo8.25.15.csv")
head(npgo)
npgo$sm <- rollapply(npgo$NDJFM, 2, mean, fill=NA, align="left")

pink$pdo <- pink$pdo2 <- pink$npgo <- pink$npgo2 <- NA

for(i in 1:nrow(pink)){
  # i <- 1
  pink$pdo[i] <- pdo$NDJFM[pdo$YEAR==pink$entry.yr[i]]
  pink$pdo2[i] <- pdo$sm[pdo$YEAR==pink$entry.yr[i]]
  
  pink$npgo[i] <- npgo$NDJFM[npgo$YEAR==pink$entry.yr[i]]
  pink$npgo2[i] <- npgo$sm[npgo$YEAR==pink$entry.yr[i]]
  
}

# final checks!


pdf("pink 2-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("pink 1-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# final checks!


pdf("pink 2-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("pink 1-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(pink$full.stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- pink[pink$full.stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# export!
write.csv(pink, "coastwide pink data.csv")

# plot all 3 winter areas!
pdf("winter sst areas.pdf", 5, 10)
par(mfrow=c(3,1), mar=c(2,2,2,0.5), oma=c(0,0,2,0))

lim=range(EBS.mean, GOA.mean, South.mean, na.rm=T)
# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)

z <- t(matrix(EBS.mean,length(EBSy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(EBSx,EBSy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("EBS")

sub <- pink.locations[pink.locations$region=="EBS",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

#
z <- t(matrix(GOA.mean,length(GOAy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(GOAx,GOAy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("GOA")

sub <- pink.locations[pink.locations$region=="GOA" & pink.locations$region2=="AK",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

#
z <- t(matrix(South.mean,length(Southy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(Southx,Southy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("South")

sub <- pink.locations[pink.locations$region2=="South",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

mtext("Winter SST areas", outer=T)
dev.off()
###################################################################################
###################################################################################

# chums!

data <- read.csv("chum_data clean.csv")
chum <- filter(data, use==1)

# this is a bit verbose, but selects keeper runs!
start <- tapply(chum$brood.yr, chum$stock, min, na.rm=T)

end <- tapply(chum$brood.yr, chum$stock, max, na.rm=T)
identical(names(start), names(end)) #T
chum.runs <- data.frame(run=as.character(names(end)), start=start, end=end)
chum.runs <- na.omit(chum.runs)
chum.runs
# data must at least run between 1972 and 2005
#### NOTE THAT THIS MIGHT BE INTERESTING TO LOOK AT W/SENSITIVITY TEST:
# DO 2005 RUNS BIAS RESULTS??

chum.use <- chum.runs[chum.runs$start <= 1972 & chum.runs$end >= 2005,]
chum.use

chum <- chum[chum$stock %in% chum.use$run,]

# add entry yr and ln.rs
chum$entry.yr <- chum$brood.yr+1
chum$ln.rs <- log(chum$recruits/chum$spawners)

head(chum)

# need lat/long!
chum$lat <- chum$long <- NA

setwd("/Users/MikeLitzow/Documents/R/NSF-GOA/MichaelMalick-pink-chum-database-af9da57")
info <- read.csv("chum_info.csv")

setwd("/Users/MikeLitzow/Documents/R/NSF-GOA")

colnames(info)
colnames(chum)
lats <- tapply(info$lat, info$stock.id, mean)
lons <- tapply(info$lon, info$stock.id, mean)

# I know there must be a better way to do this!

for(i in 1:nrow(chum)){
  #  i <- 1
  index <- chum$stock.id[i] 
  chum$lat[i] <- lats[names(lats)==index]
  chum$long[i] <- lons[names(lons)==index]
}

unique(chum$stock)

par(las=1)
plot(360-tapply(chum$long, chum$stock, mean), tapply(chum$lat, chum$stock, mean), 
     ylab="Latitude", xlab="Longitude", xaxt="n", pch=21, bg="red", ylim=c(46,68), 
     xlim=c(180,240))
axis(1, at=seq(180,240,10), labels=seq(180, 120,-10))
map('world2Hires',fill=F,xlim=c(132,250), ylim=c(20,68),add=T, lwd=1)
mtext("chum salmon locations")

#########################################
# now get sst data!
# make columns for 2-yr winter sst, 1-yr winter sst, and local sst for ocean entry
chum$w.sst.2 <- chum$w.sst.1 <- chum$loc.sst <- NA


# now load the file
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/sst.mnmean.v4.nc")

# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates

# save months and years for use later on
m <- months(d)
yrs <- years(d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

####
# begin with EBS
# EBS
# 54-62 deg. N, 186-202 deg. E
# (only 62 N as the chum runs are in the south!)
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=9)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,14,1153), count=c(9,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

EBSx <- x
EBSy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out as needed
blank <- c("N54E186", "N54E188", "N54E190",  "N54E196", "N54E198", "N54E200", "N54E202",
           "N56E186", "N56E188", "N58E186",  "N56E200", "N56E202")
SST[,blank] <- NA

# plot mean temperature pattern to check
EBS.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(50,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!


plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

###########
# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!

# ok, now need to define larger region!
chum$our.region <- "GOA"

# identify ebs runs
stocks <- unique(chum$stock)
stocks
chum$our.region[chum$stock %in% stocks[c(18,19,21)]] <- "EBS"

chum$our.region[chum$lat<52] <- "South"

# now smooth with 2-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(chum)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(chum)){
  #  i <- 1
  if(chum$our.region[i] == "EBS") chum$w.sst.2[i] <- 
      smooth[names(smooth)==chum$entry.yr[i]] else chum$w.sst.2[i] <- chum$w.sst.2[i]
      
      if(chum$our.region[i] == "EBS") chum$w.sst.1[i] <- 
          SST.win[names(SST.win)==chum$entry.yr[i]] else chum$w.sst.1[i] <- chum$w.sst.1[i]
}

############################################################
# step back - look at entire coastal range for sst values!
# 46-68 deg. N, 186-232 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=26)
y <- ncvar_get(nc, "lat", start=11, count=12)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,11,1153), count=c(26,12,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,12:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now get anomalies
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std



# choose coastal stations

coast <- c("N68E196", "N68E194", "N66E194", "N66E192", "N64E200", "N64E198", "N64E196",
           "N64E194", "N62E194", "N60E194", "N60E196", "N60E198", "N58E202", "N58E200",
           "N58E198", "N56E198", "N56E196",
           "N52E228", "N52E230", "N52E232", "N54E196", 
           "N54E198", "N54E200","N54E226", "N54E228", "N54E230", "N56E202", "N56E204", 
           "N56E206", "N56E224", "N56E226", "N56E228", "N58E204", "N58E206", "N58E208", 
           "N58E222", "N58E224", "N58E226", "N60E208", "N60E210", "N60E212", "N60E214", 
           "N60E216", "N60E218", "N60E220", "N50E232", "N50E234", "N50E236" ,"N48E236",
           "N48E234", "N46E236")

coast.sst <- SST
use <- colnames(coast.sst) %in% coast
coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# let's get an object with the lat and long for each cell!


temp1 <- strsplit(coast,"E")
temp2 <- matrix(unlist(temp1), ncol=2, byrow=TRUE)
temp3 <- strsplit(temp2[,1], "N")
temp4 <- matrix(unlist(temp3), ncol=2, byrow=TRUE)

coast.cells <- data.frame(cell=coast, lat=as.numeric(temp4[,2]), long=as.numeric(temp2[,2]))
str(coast.cells)

coast.cells <- coast.cells[order(coast.cells$lat, decreasing = T),]

# add EBS/GOA factor!
coast.cells$region <- "GOA"
change <- coast.cells$lat >= 57 & coast.cells$long <=202
coast.cells$region[change] <- "EBS"

change <- coast.cells$lat == 56 & coast.cells$long <=200
coast.cells$region[change] <- "EBS"

# and check!
EBS.coast.cells <- coast.cells$cell[coast.cells$region=="EBS"]

EBS.coast.sst <- SST

use <- colnames(EBS.coast.sst) %in% EBS.coast.cells
EBS.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(EBS.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

GOA.coast.cells <- coast.cells$cell[coast.cells$region=="GOA"]

GOA.coast.sst <- SST

use <- colnames(GOA.coast.sst) %in% GOA.coast.cells
GOA.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(GOA.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now need an object with the lat longs of each chum run!
p.lat <-tapply(chum$lat, chum$stock, mean, na.rm=T)
keep <- !is.na(p.lat)
p.lat <- p.lat[keep]

p.long <-tapply(chum$long, chum$stock, mean, na.rm=T)
keep <- !is.na(p.long)
p.long <- p.long[keep]

stocks <- unique(chum$stock)
stocks <- stocks[order(stocks)]

chum.locations <- data.frame(stock=stocks, lat=p.lat, long=p.long, region=NA, coast.cells=NA)
str(chum.locations)
rownames(chum.locations) <- 1:nrow(chum.locations)

ebs <- c(4,8,20)

chum.locations$region[-ebs] <- "GOA"
chum.locations$region[ebs] <- "EBS"
# looks good...

# change coast.cells long back to degrees W!
coast.cells$Wlong <- 180-(coast.cells$long-180)

# get cells <= 500 km away for ebs runs...
ebs.x1 <- as.matrix(cbind(chum.locations$long[chum.locations$region=="EBS"], 
                          chum.locations$lat[chum.locations$region=="EBS"]))
ebs.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="EBS"], 
                          coast.cells$lat[coast.cells$region=="EBS"]))


dist.ebs <- rdist.earth(ebs.x1, ebs.x2, miles=F)

colnames(dist.ebs) <- coast.cells$cell[coast.cells$region=="EBS"]

for(i in 1:nrow(dist.ebs)){
  # i <- 1
  run <- rownames(dist.ebs)[i]
  use <- dist.ebs[i,] <= 500
  cells <- as.vector(colnames(dist.ebs)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  chum.locations$coast.cells[chum.locations$stock==run] <- paste(cells, collapse=",")
}

###
# get cells <= 500 km away for goa runs...
goa.x1 <- as.matrix(cbind(chum.locations$long[chum.locations$region=="GOA"], 
                          chum.locations$lat[chum.locations$region=="GOA"]))
goa.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="GOA"], 
                          coast.cells$lat[coast.cells$region=="GOA"]))


dist.goa <- rdist.earth(goa.x1, goa.x2, miles=F)

colnames(dist.goa) <- coast.cells$cell[coast.cells$region=="GOA"]

for(i in 1:nrow(dist.goa)){
  # i <- 1
  run <- rownames(dist.goa)[i]
  use <- dist.goa[i,] <= 500
  cells <- as.vector(colnames(dist.goa)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  chum.locations$coast.cells[chum.locations$stock==run] <- paste(cells, collapse=",")
}

# now, some random-ish checks!
pdf("chum runs and local sst cells.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

for(i in 1:nrow(chum.locations)){
  #i <- 30
  sub <- chum.locations[i,]
  sub.sst <- SST
  temp <- strsplit(sub$coast.cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  lim <- range(colMeans(SST), na.rm = T)
  
  sub.mean <- colMeans(sub.sst)
  z <- t(matrix(sub.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=tim.colors(64), zlim=lim, xlim=c(160,240), ylim=c(44,70), ylab="", xlab="", yaxt="n", xaxt="n")
  #contour(x, y, z, add=T)  
  map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
  points(360-sub$long, sub$lat, pch=21, bg="red", cex=1.5)
  mtext(sub$stock)
}
dev.off()


# looks good!
# actually....looks GREAT!

# now, we need to get seasonal means based on the periods identified by Mueter et al. 2005!
# these are for entry year, i.e., the year after brood year!!
AK.mo <- c("Jul", "Aug", "Sep")
S.mo <- c("Jul", "Aug", "Sep") # these aren't significant in Mueter et al 2005! 
# but are the strongest correlations, consistent with the AK results (i.e., around ocean entry)

chum.locations$region2 <- "AK"
south <- c(1,5,6,11,13,14,16,21)
chum.locations$region2[south] <- "South"

# now we need to go through and select the proper months to make a local sst time series!!
stocks

for(i in 1:length(stocks)){
  # i <- 1
  sub <- chum.locations[i,]
  sub.sst <- SST.anom # so these are scaled sst time series! (mean 0, unit variance)
  temp <- strsplit(sub$coast.cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  if(sub$region2=="AK") m.use <- AK.mo else m.use <- S.mo
  keep <- months(rownames(sub.sst)) %in% m.use
  sub.sst <- sub.sst[keep,] 
  yr <- years(rownames(sub.sst))
  temp <- rowMeans(sub.sst, na.rm=T) # monthly mean anomaly for this region!
  ts <- tapply(temp, yr, mean)
  
  # now plug into the chum data frame!
  sub.sub <- chum[chum$stock==sub$stock,]
  for(j in 1:nrow(sub.sub)){
    #j <- 1
    sub.sub$loc.sst[j] <- ts[names(ts)==sub.sub$entry.yr[j]] 
    
  }
  
  # and now put the local sst time series into the chum data frame
  
  chum$loc.sst[chum$stock==sub$stock] <- sub.sub$loc.sst
}

# check time series!
pdf("chum local sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$loc.sst, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()
# these look right to me!

######################################

# now plug in winter goa temps!!

# Extract GOA SST, 54-61 deg. N, 200-226 deg. E:
# now load the file
nc <- nc_open("sst.mnmean.v4.nc")
# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates


ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=99, count=16)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(99,14,1141), count=c(16,5,length(d)), verbose = T)

# check
dim(SST) # 14 longitudes, 5 latitudes, 798 months

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

GOAx <- x
GOAy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out Bristol Bay
blank <- c("N56E200", "N58E200", "N58E202", "N56E196", "N56E198", "N58E196",
           "N58E198", "N60E196", "N60E198")
SST[,blank] <- NA

# plot mean temperature pattern to check
GOA.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,230), ylim=c(48,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!



# and another check!

plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")

yrs <- years(d)
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

# correcto!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!


# now smooth with 2-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(chum)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(chum)){
  #  i <- 1
  if(chum$our.region[i] == "GOA") chum$w.sst.2[i] <- 
      smooth[names(smooth)==chum$entry.yr[i]] else chum$w.sst.2[i] <- chum$w.sst.2[i]
      
      if(chum$our.region[i] == "GOA") chum$w.sst.1[i] <- 
          SST.win[names(SST.win)==chum$entry.yr[i]] else chum$w.sst.1[i] <- chum$w.sst.1[i]
}

###############################################
# next, get winter sst for southern stocks!

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# 46-68 deg. N, 186-232 deg. E
x <- ncvar_get(nc, "lon", start=115, count=7)
y <- ncvar_get(nc, "lat", start=19, count=4)
x; y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(115,19,1141), count=c(7,4,length(d)), verbose = T)

# check
dim(SST) # 14 longitudes, 5 latitudes, 798 months

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,4:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

Southx <- x
Southy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

blank <- c("N50E228", "N48E228", "N46E228", "N48E230", "N46E230")

SST[,blank] <- NA

# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,238), ylim=c(42,64))

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# and another check!

plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")

yrs <- years(d)
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

# correcto!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!


# now smooth with 2-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(chum)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(chum)){
  #  i <- 1
  if(chum$our.region[i] == "South") chum$w.sst.2[i] <- 
      smooth[names(smooth)==chum$entry.yr[i]] else chum$w.sst.2[i] <- chum$w.sst.2[i]
      
      if(chum$our.region[i] == "South") chum$w.sst.1[i] <- 
          SST.win[names(SST.win)==chum$entry.yr[i]] else chum$w.sst.1[i] <- chum$w.sst.1[i]
}

# check!

pdf("chum 2-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("chum 1-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.1, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# all looks good!

# finally, add PDO and NPGO

pdo <- read.csv("pdo.8.25.15.csv")
head(pdo)
pdo$sm <- rollapply(pdo$NDJFM, 2, mean, fill=NA, align="left")

npgo <- read.csv("npgo8.25.15.csv")
head(npgo)
npgo$sm <- rollapply(npgo$NDJFM, 2, mean, fill=NA, align="left")

chum$pdo <- chum$pdo2 <- chum$npgo <- chum$npgo2 <- NA

for(i in 1:nrow(chum)){
  # i <- 1
  chum$pdo[i] <- pdo$NDJFM[pdo$YEAR==chum$entry.yr[i]]
  chum$pdo2[i] <- pdo$sm[pdo$YEAR==chum$entry.yr[i]]
  
  chum$npgo[i] <- npgo$NDJFM[npgo$YEAR==chum$entry.yr[i]]
  chum$npgo2[i] <- npgo$sm[npgo$YEAR==chum$entry.yr[i]]
  
}

# final checks!


pdf("chum 2-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("chum 1-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# final checks!


pdf("chum 2-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("chum 1-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(chum$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- chum[chum$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# export!
write.csv(chum, "coastwide chum data.csv")

# plot all 3 winter areas!
pdf("winter sst areas.pdf", 5, 10)
par(mfrow=c(3,1), mar=c(2,2,2,0.5), oma=c(0,0,2,0))

lim=range(EBS.mean, GOA.mean, South.mean, na.rm=T)
# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)

z <- t(matrix(EBS.mean,length(EBSy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(EBSx,EBSy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("EBS")

sub <- chum.locations[chum.locations$region=="EBS",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

#
z <- t(matrix(GOA.mean,length(GOAy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(GOAx,GOAy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("GOA")

sub <- chum.locations[chum.locations$region=="GOA" & chum.locations$region2=="AK",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

#
z <- t(matrix(South.mean,length(Southy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(Southx,Southy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("South")

sub <- chum.locations[chum.locations$region2=="South",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

mtext("Winter SST areas", outer=T)
dev.off()

#############################################################################
# now sockeyes!
#############################################################################
### NB! Alagnak seems to somehow have been excluded in the data that I sent to Bethany in early June!

# I have (quickly) updated with most recent available data

#data <- read.csv("sockeye_data clean.csv")

data <- read.csv("sockeye_data_extra_years.csv")

head(data)

# limit to the good stuff
data <- data[data$use==1,]



ender <- tapply(data$brood.yr, data$stock, max)
starter <- tapply(data$brood.yr, data$stock, min)
overview <- cbind(starter, ender)


sockeye <- filter(data, use==1)
# this is a bit verbose, but selects keeper runs!
start <- tapply(sockeye$brood.yr, sockeye$stock, min, na.rm=T)

end <- tapply(sockeye$brood.yr, sockeye$stock, max, na.rm=T)
identical(names(start), names(end)) #T
sockeye.runs <- data.frame(run=as.character(names(end)), start=start, end=end)
sockeye.runs <- na.omit(sockeye.runs)

# data must at least run between 1972 and 2005 
# (more liberal at the end, as sockeye take so long to go to sea)
sockeye.use <- sockeye.runs[sockeye.runs$start <= 1972 & sockeye.runs$end >= 2005,]
sockeye.use

# and... drop the runs that Rich Brenner identified as hatchery subsidized,
# as well as the combined Chignik

drop <- c("Coghill", "Eshamy", "Copper", "Chignik.combined")
sockeye.use <- sockeye.use[!sockeye.use$run %in% drop,]
sockeye <- sockeye[sockeye$stock %in% sockeye.use$run,]

# add entry yr and ln.rs
sockeye$entry.yr <- sockeye$brood.yr+2 # this is a nominal entry age!
# Mueter et al. 2002 say all southern stocks enter at age 2; AK at 2-3
sockeye$ln.rs <- log(sockeye$recruits/sockeye$spawners)

head(sockeye)

# need lat/long!
sockeye$lat <- sockeye$long <- NA

info <- read.csv("sockeye.lat.long.csv")



colnames(info)
colnames(sockeye)

# and! limit to 1960-present following GOA analysis
sockeye <- sockeye[sockeye$brood.yr >= 1960,]

# I know there must be a better way to do this!

for(i in 1:nrow(sockeye)){
  # i <- 1
  index <- as.character(sockeye$stock[i]) 
  sockeye$lat[i] <- info$lat[as.character(info$stock)==index]
  sockeye$long[i] <- info$long[as.character(info$stock)==index]
}

unique(sockeye$stock)

par(las=1)
plot(360-tapply(sockeye$long, sockeye$stock, mean), tapply(sockeye$lat, sockeye$stock, mean), 
     ylab="Latitude", xlab="Longitude", xaxt="n", pch=21, bg="red", ylim=c(46,68), 
     xlim=c(180,240))
axis(1, at=seq(180,240,10), labels=seq(180, 120,-10))
map('world2Hires',fill=F,xlim=c(132,250), ylim=c(20,68),add=T, lwd=1)
mtext("sockeye salmon locations")

#########################################
# now get sst data!
# make columns for 2-yr winter sst, 1-yr winter sst, and local sst for ocean entry (and adding 3-yr as well!)
sockeye$w.sst.3 <- sockeye$w.sst.2 <- sockeye$w.sst.1 <- sockeye$loc.sst <- NA


# now load the file
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/sst.mnmean.v4.nc")

# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates

# save months and years for use later on
m <- months(d)
yrs <- years(d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

####
# begin with EBS
# EBS
# 54-62 deg. N, 186-202 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=9)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,14,1153), count=c(9,5,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

EBSx <- x
EBSy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out as needed
blank <- c("N54E186", "N54E188", "N54E190",  "N54E196", "N54E198", "N54E200", "N54E202",
           "N56E186", "N56E188", "N58E186",  "N56E200", "N56E202")
SST[,blank] <- NA

# plot mean temperature pattern to check
EBS.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(50,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!


plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

###########
# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!

# ok, now need to define larger region!
sockeye$our.region <- "GOA"

# identify ebs runs
stocks <- unique(sockeye$stock)

sockeye$our.region[sockeye$stock %in% stocks[c(1,9,13,16,20,21,25,26,28)]] <- "EBS"

sockeye$our.region[sockeye$lat<52] <- "South"

# now smooth with 2-yr running mean
win.3 <- rollmean(SST.win, 3, fill=NA)
smooth <- rollmean(SST.win, 2, align="left", fill=NA)

head(sockeye)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(sockeye)){
  #  i <- 1
  if(sockeye$our.region[i] == "EBS") sockeye$w.sst.3[i] <- 
      win.3[names(win.3)==sockeye$entry.yr[i]] else sockeye$w.sst.3[i] <- sockeye$w.sst.3[i]
      
      if(sockeye$our.region[i] == "EBS") sockeye$w.sst.2[i] <- 
          smooth[names(smooth)==sockeye$entry.yr[i]] else sockeye$w.sst.2[i] <- sockeye$w.sst.2[i]
          
          if(sockeye$our.region[i] == "EBS") sockeye$w.sst.1[i] <- 
              SST.win[names(SST.win)==sockeye$entry.yr[i]] else sockeye$w.sst.1[i] <- sockeye$w.sst.1[i]
}

############################################################
# step back - look at entire coastal range for sst values!
# 46-68 deg. N, 186-232 deg. E
ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=94, count=26)
y <- ncvar_get(nc, "lat", start=11, count=12)
x; y

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(94,11,1153), count=c(26,12,length(d)), verbose = T)

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,12:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now get anomalies
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std



# choose coastal stations

coast <- c("N68E196", "N68E194", "N66E194", "N66E192", "N64E200", "N64E198", "N64E196",
           "N64E194", "N62E194", "N60E194", "N60E196", "N60E198", "N58E202", "N58E200",
           "N58E198", "N56E198", "N56E196",
           "N52E228", "N52E230", "N52E232", "N54E196", 
           "N54E198", "N54E200","N54E226", "N54E228", "N54E230", "N56E202", "N56E204", 
           "N56E206", "N56E224", "N56E226", "N56E228", "N58E204", "N58E206", "N58E208", 
           "N58E222", "N58E224", "N58E226", "N60E208", "N60E210", "N60E212", "N60E214", 
           "N60E216", "N60E218", "N60E220", "N50E232", "N50E234", "N50E236" ,"N48E236",
           "N48E234", "N46E236")

coast.sst <- SST
use <- colnames(coast.sst) %in% coast
coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# let's get an object with the lat and long for each cell!


temp1 <- strsplit(coast,"E")
temp2 <- matrix(unlist(temp1), ncol=2, byrow=TRUE)
temp3 <- strsplit(temp2[,1], "N")
temp4 <- matrix(unlist(temp3), ncol=2, byrow=TRUE)

coast.cells <- data.frame(cell=coast, lat=as.numeric(temp4[,2]), long=as.numeric(temp2[,2]))
str(coast.cells)

coast.cells <- coast.cells[order(coast.cells$lat, decreasing = T),]

# add EBS/GOA factor!
coast.cells$region <- "GOA"
change <- coast.cells$lat >= 57 & coast.cells$long <=202
coast.cells$region[change] <- "EBS"

change <- coast.cells$lat == 56 & coast.cells$long <=200
coast.cells$region[change] <- "EBS"

# and check!
EBS.coast.cells <- coast.cells$cell[coast.cells$region=="EBS"]

EBS.coast.sst <- SST

use <- colnames(EBS.coast.sst) %in% EBS.coast.cells
EBS.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(EBS.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

GOA.coast.cells <- coast.cells$cell[coast.cells$region=="GOA"]

GOA.coast.sst <- SST

use <- colnames(GOA.coast.sst) %in% GOA.coast.cells
GOA.coast.sst[,!use] <- NA
# plot mean temperature pattern to check
coast.mean <- colMeans(GOA.coast.sst)
z <- t(matrix(coast.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(44,70))
#contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
# looks good!

# now need an object with the lat longs of each sockeye run!
p.lat <-tapply(sockeye$lat, sockeye$stock, mean, na.rm=T)
keep <- !is.na(p.lat)
p.lat <- p.lat[keep]

p.long <-tapply(sockeye$long, sockeye$stock, mean, na.rm=T)
keep <- !is.na(p.long)
p.long <- p.long[keep]

stocks <- unique(sockeye$stock)
stocks <- stocks[order(stocks)]

sockeye.locations <- data.frame(stock=stocks, lat=p.lat, long=p.long, region=NA, coast.cells=NA)
str(sockeye.locations)
rownames(sockeye.locations) <- 1:nrow(sockeye.locations)

ebs <- c(1,9,13,16,20,21,25,26,28)

sockeye.locations$region[-ebs] <- "GOA"
sockeye.locations$region[ebs] <- "EBS"
# looks good...

# change coast.cells long back to degrees W!
coast.cells$Wlong <- 180-(coast.cells$long-180)

# get cells <= 500 km away for ebs runs...
ebs.x1 <- as.matrix(cbind(sockeye.locations$long[sockeye.locations$region=="EBS"], 
                          sockeye.locations$lat[sockeye.locations$region=="EBS"]))
ebs.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="EBS"], 
                          coast.cells$lat[coast.cells$region=="EBS"]))


dist.ebs <- rdist.earth(ebs.x1, ebs.x2, miles=F)

colnames(dist.ebs) <- coast.cells$cell[coast.cells$region=="EBS"]

for(i in 1:nrow(dist.ebs)){
  # i <- 1
  run <- rownames(dist.ebs)[i]
  use <- dist.ebs[i,] <= 500
  cells <- as.vector(colnames(dist.ebs)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  sockeye.locations$coast.cells[sockeye.locations$stock==run] <- paste(cells, collapse=",")
}

###
# get cells <= 500 km away for goa runs...
goa.x1 <- as.matrix(cbind(sockeye.locations$long[sockeye.locations$region=="GOA"], 
                          sockeye.locations$lat[sockeye.locations$region=="GOA"]))
goa.x2 <- as.matrix(cbind(coast.cells$Wlong[coast.cells$region=="GOA"], 
                          coast.cells$lat[coast.cells$region=="GOA"]))


dist.goa <- rdist.earth(goa.x1, goa.x2, miles=F)

colnames(dist.goa) <- coast.cells$cell[coast.cells$region=="GOA"]

for(i in 1:nrow(dist.goa)){
  # i <- 1
  run <- rownames(dist.goa)[i]
  use <- dist.goa[i,] <= 500
  cells <- as.vector(colnames(dist.goa)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  sockeye.locations$coast.cells[sockeye.locations$stock==run] <- paste(cells, collapse=",")
}

# now, some random-ish checks!
pdf("sockeye runs and local sst cells.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

for(i in 1:nrow(sockeye.locations)){
  #i <- 30
  sub <- sockeye.locations[i,]
  sub.sst <- SST
  temp <- strsplit(sub$coast.cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  lim <- range(colMeans(SST), na.rm = T)
  
  sub.mean <- colMeans(sub.sst)
  z <- t(matrix(sub.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=tim.colors(64), zlim=lim, xlim=c(160,240), ylim=c(44,70), ylab="", xlab="", yaxt="n", xaxt="n")
  #contour(x, y, z, add=T)  
  map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,72),add=T, lwd=2)
  points(360-sub$long, sub$lat, pch=21, bg="red", cex=1.5)
  mtext(sub$stock)
}
dev.off()


# looks good!
# actually....looks GREAT!

# now, we need to get seasonal means based on the periods identified by Mueter et al. 2005!
# these are years relative to brood year!!
AK.mo.l1 <- c("Jan", "Feb")
AK.mo.l2 <- c("Jan", "Feb", "Mar", "Apr", "May")
AK.mo.l3 <- c("Jan", "Feb", "Mar", "Apr", "May", "Aug", "Sep", "Oct", "Nov", "Dec")


S.mo.l1 <- c("Jan", "Mar", "Apr", "May")
S.mo.l2 <- c("Jan", "Mar", "Apr")

sockeye.locations$region2 <- "AK"
south <- c(3,5,6,9,11,16,17,21,22,26)
sockeye.locations$region2[south] <- "South"

# now we need to go through and select the proper months to make a local sst time series!!
stocks

for(i in 1:length(stocks)){
  #  i <- 3
  sub <- sockeye.locations[i,]
  
  if(sub$region2=="AK"){
    # get correct cells
    sub.sst <- SST.anom # so these are scaled sst time series! (mean 0, unit variance)
    temp <- strsplit(sub$coast.cells,",")
    use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
    
    keep <- colnames(sub.sst) %in% use
    sub.sst[,!keep] <- NA
    
    # year 1 (wrt b.y.)
    keep <- months(rownames(sub.sst)) %in% AK.mo.l1
    sub.sst1 <- sub.sst
    sub.sst1[!keep,] <- NA 
    
    keep <- months(rownames(sub.sst)) %in% AK.mo.l2
    sub.sst2 <- sub.sst
    sub.sst2[!keep,] <- NA
    
    keep <- months(rownames(sub.sst)) %in% AK.mo.l3
    sub.sst3 <- sub.sst
    sub.sst3[!keep,] <- NA
    
    # now get means!
    years <- years(row.names(sub.sst))
    yr <- data.frame(yr=1960:2011)
    yr$mean.loc <- NA
    for(y in 1960:2011){
      # y <- 1960
      temp <- c(as.vector(sub.sst2[years==y,]), 
                as.vector(sub.sst1[years==(y-1),]), 
                as.vector(sub.sst3[years==(y+1),]), na.rm=T)
      temp <- na.omit(temp)
      yr$mean.loc[yr$yr==y] <- mean(temp)
    }
    
    # now plug into the sockeye data frame!
    sub.sub <- sockeye[sockeye$stock==sub$stock,]
    for(j in 1:nrow(sub.sub)){
      #j <- 1
      sub.sub$loc.sst[j] <- yr$mean.loc[yr$yr==sub.sub$entry.yr[j]] 
      
    }
    
    # and now put the local sst time series into the sockeye data frame
    
    sockeye$loc.sst[sockeye$stock==sub$stock] <- sub.sub$loc.sst
    # temp <- rowMeans(sub.sst, na.rm=T) # monthly mean anomaly for this region!
    # ts <- tapply(temp, yr, mean)
    
  }
  
  # now the southern stocks!!!
  if(sub$region2=="South"){
    # get correct cells
    sub.sst <- SST.anom # so these are scaled sst time series! (mean 0, unit variance)
    temp <- strsplit(sub$coast.cells,",")
    use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
    
    keep <- colnames(sub.sst) %in% use
    sub.sst[,!keep] <- NA
    
    # year 1 (wrt b.y.)
    keep <- months(rownames(sub.sst)) %in% S.mo.l1
    sub.sst1 <- sub.sst
    sub.sst1[!keep,] <- NA 
    
    keep <- months(rownames(sub.sst)) %in% S.mo.l2
    sub.sst2 <- sub.sst
    sub.sst2[!keep,] <- NA
    
    # now get means!
    years <- years(row.names(sub.sst))
    yr <- data.frame(yr=1960:2011)
    yr$mean.loc <- NA
    for(y in 1960:2011){
      # y <- 1960
      temp <- c(as.vector(sub.sst2[years==y,]), 
                as.vector(sub.sst1[years==(y-1),]), na.rm=T)
      temp <- na.omit(temp)
      yr$mean.loc[yr$yr==y] <- mean(temp)
    }
    
    # now plug into the sockeye data frame!
    sub.sub <- sockeye[sockeye$stock==sub$stock,]
    for(j in 1:nrow(sub.sub)){
      #j <- 1
      sub.sub$loc.sst[j] <- yr$mean.loc[yr$yr==sub.sub$entry.yr[j]] 
      
    }
    
    # and now put the local sst time series into the sockeye data frame
    
    sockeye$loc.sst[sockeye$stock==sub$stock] <- sub.sub$loc.sst
    # temp <- rowMeans(sub.sst, na.rm=T) # monthly mean anomaly for this region!
    # ts <- tapply(temp, yr, mean)
    
  }
  
}  
# check time series!
pdf("sockeye local sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$loc.sst, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()
# these look right to me!

######################################

# now plug in winter goa temps!!

# Extract GOA SST, 54-61 deg. N, 200-226 deg. E:
# now load the file
nc <- nc_open("sst.mnmean.v4.nc")
# get info
nc

# view dates (middle of month):
ncvar_get(nc, "time")  # Days since January 1, 1800 (see CDC documentation)

# assign dates
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# we will start with Jan 1950 and end with most recent month
d[c(1153,length(d))] # check - those are correct
d <- d[1153:length(d)] # restrict to the desired dates


ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

x <- ncvar_get(nc, "lon", start=99, count=16)
y <- ncvar_get(nc, "lat", start=14, count=5)
x; y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(99,14,1141), count=c(16,5,length(d)), verbose = T)

# check
dim(SST) # 14 longitudes, 5 latitudes, 798 months

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,5:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

GOAx <- x
GOAy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

#blank out Bristol Bay
blank <- c("N56E200", "N58E200", "N58E202", "N56E196", "N56E198", "N58E196",
           "N58E198", "N60E196", "N60E198")
SST[,blank] <- NA

# plot mean temperature pattern to check
GOA.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,230), ylim=c(48,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!



# and another check!

plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")

yrs <- years(d)
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

# correcto!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!


# now smooth with 2-yr  and 3-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)
win.3 <- rollmean(SST.win, 3, fill=NA)

head(sockeye)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(sockeye)){
  #  i <- 1
  
  if(sockeye$our.region[i] == "GOA") sockeye$w.sst.3[i] <- 
      win.3[names(win.3)==sockeye$entry.yr[i]] else sockeye$w.sst.3[i] <- sockeye$w.sst.3[i]
      
      if(sockeye$our.region[i] == "GOA") sockeye$w.sst.2[i] <- 
          smooth[names(smooth)==sockeye$entry.yr[i]] else sockeye$w.sst.2[i] <- sockeye$w.sst.2[i]
          
          if(sockeye$our.region[i] == "GOA") sockeye$w.sst.1[i] <- 
              SST.win[names(SST.win)==sockeye$entry.yr[i]] else sockeye$w.sst.1[i] <- sockeye$w.sst.1[i]
}

###############################################
# next, get winter sst for southern stocks!

ncvar_get(nc, "lon")     # view longitudes (degrees East)
ncvar_get(nc, "lat")     # view latitudes

# 46-68 deg. N, 186-232 deg. E
x <- ncvar_get(nc, "lon", start=115, count=7)
y <- ncvar_get(nc, "lat", start=19, count=4)
x; y # check

# get required sst data
SST <- ncvar_get(nc, "sst", start=c(115,19,1141), count=c(7,4,length(d)), verbose = T)

# check
dim(SST) # 14 longitudes, 5 latitudes, 798 months

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array
SST <- SST[,4:1,]  # Reverse order of latitudes to be increasing for plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

Southx <- x
Southy <- y

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

blank <- c("N50E228", "N48E228", "N46E228", "N48E230", "N46E230")

SST[,blank] <- NA

# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,238), ylim=c(42,64))

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# and another check!

plot(1:nrow(SST), rowMeans(SST, na.rm=T), type="l")

yrs <- years(d)
ann.mean <- tapply(rowMeans(SST, na.rm=T), yrs, mean)

plot(1950:2016, ann.mean[1:67], type="l")

# correcto!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)

add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# define seasons
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
SST.win <- tapply(SST.mean, win.yrs, mean)
f3 <- function(x) sum(!is.na(x)) # function to count number of cases present
win.count <- tapply(SST.mean, win.yrs, f3)==5 # logical vector to identify complete winters!
SST.win <- SST.win[win.count]

# check with a plot
plot(names(SST.win), SST.win, type="b")
# looks right!


# now smooth with 2-yr and 3-yr running mean
smooth <- rollmean(SST.win, 2, align="left", fill=NA)
win.3 <- rollmean(SST.win, 3, fill=NA)
head(sockeye)

# now set smooth sst to the winter before and the winter after ocean entry
for(i in 1:nrow(sockeye)){
  #  i <- 1
  
  if(sockeye$our.region[i] == "South") sockeye$w.sst.3[i] <- 
      win.3[names(win.3)==sockeye$entry.yr[i]] else sockeye$w.sst.3[i] <- sockeye$w.sst.3[i]
      
      if(sockeye$our.region[i] == "South") sockeye$w.sst.2[i] <- 
          smooth[names(smooth)==sockeye$entry.yr[i]] else sockeye$w.sst.2[i] <- sockeye$w.sst.2[i]
          
          if(sockeye$our.region[i] == "South") sockeye$w.sst.1[i] <- 
              SST.win[names(SST.win)==sockeye$entry.yr[i]] else sockeye$w.sst.1[i] <- sockeye$w.sst.1[i]
}

# check!

pdf("sockeye 3-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.3, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()


pdf("sockeye 2-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("sockeye 1-yr winter sst time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$w.sst.1, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# all looks good!

# finally, add PDO and NPGO

pdo <- read.csv("pdo.8.25.15.csv")
head(pdo)
pdo$sm <- rollapply(pdo$NDJFM, 2, mean, fill=NA, align="left")
pdo$pdo3 <- rollapply(pdo$NDJFM, 3, mean, fill=NA)

npgo <- read.csv("npgo8.25.15.csv")
head(npgo)
npgo$sm <- rollapply(npgo$NDJFM, 2, mean, fill=NA, align="left")
npgo$npgo3 <- rollapply(npgo$NDJFM, 3, mean, fill=NA)

sockeye$pdo <- sockeye$pdo2 <- sockeye$pdo3 <- sockeye$npgo <- sockeye$npgo2 <- sockeye$npgo3 <- NA

for(i in 1:nrow(sockeye)){
  # i <- 1
  sockeye$pdo[i] <- pdo$NDJFM[pdo$YEAR==sockeye$entry.yr[i]]
  sockeye$pdo2[i] <- pdo$sm[pdo$YEAR==sockeye$entry.yr[i]]
  sockeye$pdo3[i] <- pdo$pdo3[pdo$YEAR==sockeye$entry.yr[i]]
  
  sockeye$npgo[i] <- npgo$NDJFM[npgo$YEAR==sockeye$entry.yr[i]]
  sockeye$npgo2[i] <- npgo$sm[npgo$YEAR==sockeye$entry.yr[i]]
  sockeye$npgo3[i] <- npgo$npgo3[npgo$YEAR==sockeye$entry.yr[i]]
}

# final checks!


pdf("sockeye 3-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo3, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()


pdf("sockeye 2-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("sockeye 1-yr winter pdo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$pdo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# final checks!
pdf("sockeye 3-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo3, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("sockeye 2-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 3
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo2, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

pdf("sockeye 1-yr winter npgo time series.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

stocks <- unique(sockeye$stock)

for(i in 1:length(stocks)){
  # i <- 10
  sub <- sockeye[sockeye$stock==stocks[i],]
  plot(sub$entry.yr, sub$npgo, type="b", ylab="", xlab="")
  
  mtext(stocks[i])
}
dev.off()

# export!
write.csv(sockeye, "coastwide sockeye data UPDATED 6-20-17.csv")

# plot all 3 winter areas!
pdf("winter sst areas.pdf", 5, 10)
par(mfrow=c(3,1), mar=c(2,2,2,0.5), oma=c(0,0,2,0))

lim=range(EBS.mean, GOA.mean, South.mean, na.rm=T)
# plot mean temperature pattern to check
South.mean <- SST.mean <- colMeans(SST)

z <- t(matrix(EBS.mean,length(EBSy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(EBSx,EBSy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("EBS")

sub <- sockeye.locations[sockeye.locations$region=="EBS",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

#
z <- t(matrix(GOA.mean,length(GOAy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(GOAx,GOAy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("GOA")

sub <- sockeye.locations[sockeye.locations$region=="GOA" & sockeye.locations$region2=="AK",]
unique(sub$stock)
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")
##NB - ERRORS HERE WITH REGIONAL DEFINITIONS, NEEDS TO BE FIXED! (6/20/17)
#
z <- t(matrix(South.mean,length(Southy)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(Southx,Southy,z, col=tim.colors(64), xlim=c(180,238), ylim=c(42,70), zlim=lim, ylab="", xlab="",
      yaxt="n", xaxt="n")

map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,70),add=T, lwd=2)
mtext("South")

sub <- sockeye.locations[sockeye.locations$region2=="South",]
points(360-sub$long, sub$lat, pch=21, cex=1.5, bg="red")

mtext("Winter SST areas", outer=T)
dev.off()







#####################
# slightly different tack from above! 
# just use the local definitions used by Mueter et al. 2005!


# now identify the cells used for each region by Mueter et al. 2005!

BB <- c("N58E202", "N58E200", "N58E198", "N56E198", "N56E196")
SEAK <- c("N58E224", "N58E226","N56E226", "N56E228")
cGOA <- c("N60E210", "N60E212", "N60E214")
CI <- "N60E208"
Kod <- c("N58E206", "N58E208", "N56E204", "N56E206")
AkPen <- c("N54E196", "N54E198")

reg.cells <- matrix(nrow=5, ncol=4) # set up matrix of regional cells
reg.cells[1,1:4] <- SEAK
reg.cells[2,1:3] <- cGOA
reg.cells[3,1] <- CI
reg.cells[4,1:4] <- Kod
reg.cells[5,1:2] <- AkPen

row.names(reg.cells) <- c("SEAK", "cGOA", "CI", "Kod", "AkPen")

# and check each!!
yy <- 51:63
xx <- seq(190, 232, length.out=length(yy))
plot(xx,yy, type="n")
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
pl.cl <- cells[cells$cell %in% AkPen,]
points(pl.cl[,3], pl.cl[,2], pch=19, col="red", cex=2)

# they all look good!

# now get monthly anomalies for SST at each cell to enable valid averages across months
sst.m <- months(d)  # Extracts months from the date vector
sst.m  # Result is a factor
f <- function(x) tapply(x, sst.m, mean)  # function to compute monthly means for a single time series
sst.mu <- apply(SST, 2, f)	# Compute monthly means for each time series (location)
sst.mu
sst.mu <- sst.mu[rep(1:12, length(d)/12),]  # Replicate means matrix for each year at each location

SST.anom <- SST - sst.mu  # Compute matrix of monthly anomalies! 
# so, for pink, the months that Mueter et al. 2005 found as significant correlations were Jan-Sep

keep.m <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")
d
m <- months(d)
keep <- m %in% keep.m

SST.pink <- SST.anom[keep,] # so these are anomalies only the months that we want for all pink runs

y <- years(rownames(SST.pink))

# now need to go through and add the average for SST.pink for the correct region, for entry year, for each run!

# add entry year to pink data
pink$entry.yr <- pink$brood.yr+1

# and, we need  df of pink runs by region!!
run.reg <- data.frame(stock = unique(pink$stock), region = NA)
run.reg[,2] <- c(rep("SEAK", 6),"cGOA", "cGOA", rep("AkPen", 14))

# check
run.reg
# looks good!


pink$loc.sst <- NA
for(i in 1:length(unique(pink$stock))){ # loop through each run!
  # i <- 1
  region <- run.reg[run.reg$stock %in% unique(pink$stock)[i],2]
  # choose the correct cells!
  reg <- reg.cells[rownames(reg.cells) %in% region]
  
  temp.SST <- rowMeans(SST.pink[,colnames(SST.pink) %in% reg])
  loc.sst <- tapply(temp.SST, y, mean)
  
  # now pull out the run to fit local SST data!
  temp <- pink[pink$stock %in% unique(pink$stock)[i],]
  
  for(j in 1:nrow(temp)){ # loop through each year of data
    # j <- 1
    temp$loc.sst[j] <- loc.sst[names(loc.sst) %in% temp$entry.yr[j]]	
  }
  
  # now add loc.sst back into pink data frame!
  pink$loc.sst[pink$stock %in% unique(pink$stock)[i]] <- temp$loc.sst
}

# check
head(temp, n=20)
tail(pink, n=20)
# looks good!!

# add a few more columns for mixed-effects models!
pink$ln.rs <- log(pink$recruits/pink$spawners)
pink$era <- NA
pink$era[pink$entry.yr <= 1988] <- "early"
pink$era[pink$entry.yr > 1988] <- "late"



