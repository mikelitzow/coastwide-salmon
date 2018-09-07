library(geodist)
library(dplyr)
library(maps)
library(mapdata)

# load the pink salmon data!

data <- read.csv("coastwide pink data CHECKED.csv") 
# these are the data that I provided to Bethany/Lorenzo/Patri

# restrict to GOA, and use=1 to make sure we only include good data, limit to data from 1960 onwards...
pink <- filter(data, use==1)

# take a look at the included runs
unique(pink$full.stock)
tapply(pink$entry.yr, pink$full.stock, min)

tapply(pink$entry.yr, pink$full.stock, min)
tapply(pink$entry.yr, pink$full.stock, max)


# drop PWS
drop <- grep("PWS", pink$full.stock)
pink <- pink[-drop,]
# set "era" factor
pink$era <- "early"
pink$era[pink$entry.yr > 1988] <- "late"

# and log recruits
pink$ln.recr <- log(pink$recruits)

# now get alongshore distance, S -> N

pink.runs <- pink %>%
  group_by(stock) %>%
  summarise(stock.id = mean(stock.id), lat=mean(lat), long=mean(360-long))


# plot locations
plot(pink.runs$long, pink.runs$lat, type="n")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
map('world2Hires', 'Canada',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
text(pink.runs$long, pink.runs$lat, labels = pink.runs$stock.id, col="red", cex=0.8)

# and southern stocks only
plot(pink.runs$long, pink.runs$lat, type="n", xlim=c(230,240), ylim=c(45,52))
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
map('world2Hires', 'Canada',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
text(pink.runs$long, pink.runs$lat, labels = pink.runs$stock.id, col="red", cex=0.8)

# now define alongshore distance from stock 201

pink.runs$dist <- NA

pink.runs$dist[pink.runs$stock.id==201] <- 0

temp <- filter(pink.runs, stock.id==201)
p1 <- c(360-temp$long, temp$lat)

xx <- c(202:207,209, 221, 222, 223)

for(ii in 1:length(xx)){
  i <- xx[ii]
  temp2 <- filter(pink.runs, stock.id==i)
  p2 <- c(360-temp2$long, temp2$lat) 
  dd <- rbind(p1,p2)
  pink.runs$dist[pink.runs$stock.id==i] <- geodist(p1,p2)/1000
  
}


# next set
temp <- filter(pink.runs, stock.id==222)
p1 <- c(360-temp$long, temp$lat)

xx <- 237

for(ii in 1:length(xx)){
  i <- xx[ii]
  temp2 <- filter(pink.runs, stock.id==i)
  p2 <- c(360-temp2$long, temp2$lat) 
  dd <- rbind(p1,p2)
  pink.runs$dist[pink.runs$stock.id==i] <- geodist(p1,p2)/1000 + pink.runs$dist[pink.runs$stock.id==222]
  
}

# and next
temp <- filter(pink.runs, stock.id==237)
p1 <- c(360-temp$long, temp$lat)

xx <- c(236,235,238:241)

for(ii in 1:length(xx)){
  i <- xx[ii]
  temp2 <- filter(pink.runs, stock.id==i)
  p2 <- c(360-temp2$long, temp2$lat) 
  dd <- rbind(p1,p2)
  pink.runs$dist[pink.runs$stock.id==i] <- geodist(p1,p2)/1000 + pink.runs$dist[pink.runs$stock.id==237]
  
}

# and the last 3
temp <- filter(pink.runs, stock.id==241)
p1 <- c(360-temp$long, temp$lat)

xx <- 242

for(ii in 1:length(xx)){
  i <- xx[ii]
  temp2 <- filter(pink.runs, stock.id==i)
  p2 <- c(360-temp2$long, temp2$lat) 
  dd <- rbind(p1,p2)
  pink.runs$dist[pink.runs$stock.id==i] <- geodist(p1,p2)/1000 + pink.runs$dist[pink.runs$stock.id==241]
}


temp <- filter(pink.runs, stock.id==242)
p1 <- c(360-temp$long, temp$lat)

xx <- 245

for(ii in 1:length(xx)){
  i <- xx[ii]
  temp2 <- filter(pink.runs, stock.id==i)
  p2 <- c(360-temp2$long, temp2$lat) 
  dd <- rbind(p1,p2)
  pink.runs$dist[pink.runs$stock.id==i] <- geodist(p1,p2)/1000 + pink.runs$dist[pink.runs$stock.id==242]
}


temp <- filter(pink.runs, stock.id==245)
p1 <- c(360-temp$long, temp$lat)

xx <- 246

for(ii in 1:length(xx)){
  i <- xx[ii]
  temp2 <- filter(pink.runs, stock.id==i)
  p2 <- c(360-temp2$long, temp2$lat) 
  dd <- rbind(p1,p2)
  pink.runs$dist[pink.runs$stock.id==i] <- geodist(p1,p2)/1000 + pink.runs$dist[pink.runs$stock.id==245]
}

# and add back to pink data!

pink$dist <- NA

for(i in 1:nrow(pink)){
  # i <- 1
  pink$dist[i] <- pink.runs$dist[pink.runs$stock.id == pink$stock.id[i]]
}

# and save
write.csv(pink, "pink.data.csv")
