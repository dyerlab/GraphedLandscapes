rm(list=ls())

require(raster)
require(rgdal)
require(raster)
require(GeNetIt) 
require(spdep) 
require(maptools) 
require(RANN) 
require(sp) 
require(spatialEco) 
require(GeNetIt)
#data for exercise 
data(dps) 
data(ralu.site)  
data(rasters)
# Coerce #SpatialPixelsDataFrame to #raster stack 
xvars <- stack(rasters[-6]) 
land.cov <- stack(rasters[6]) 
#plot sample locations 

plot(ralu.site)
plot( xvars )


ralu.site@data <- data.frame(ralu.site@data, extract(xvars[[c(1,2)]], ralu.site))
head(ralu.site@data)

dist.graph <- knn.graph(ralu.site, row.names = ralu.site@data[,"SiteName"])
dist.graph@data$from.to <- paste(dist.graph$i, dist.graph$j, sep=".") 

dps$from.to <- paste(dps$FROM_SITE, dps$TO_SITE, sep=".") 

dist.graph <- merge(dist.graph, dps, by = "from.to")

dim( dist.graph@data )
dist.graph@data <- dist.graph@data[,-c(7,8)]
dim( dist.graph@data )


# Some other crap

stats <- graph.statistics(dist.graph, 
                          r = xvars, 
                          d=30, 
                          stats = c("min", "mean", "max", "var"),
                          sp = FALSE)
dist.graph@data <- data.frame(dist.graph@data, stats)

wet.pct <- function(x) {
  x<-ifelse(x==11| x==12|x==90|x==95,1,0)
  prop.table(table(x))[2] 
  }
lc.stats <- graph.statistics(dist.graph, r = land.cov, d=30, stats = "wet.pct")
lc.stats[is.na(lc.stats)] <- 0
dist.graph@data <- data.frame(dist.graph@data, lc.stats)
str(dist.graph@data)




# set up the gravity model stuff
from <- ralu.site@data[,c(1,6,8,18,19)]
names(from)[2:ncol(from)] <- paste("f", names(from)[2:ncol(from)], sep=".")
to <- ralu.site@data[,c(1,6,8,18,19)]
names(to)[2:ncol(to)] <- paste("t", names(to)[2:ncol(to)], sep=".")

site <- data.frame(from,to) 
site <- site[,-(dim(to)[2]+1)]
cdata <- merge(dist.graph, site, by.x="from_ID", by.y="SiteName") 
cdata$Dps <- 1 - cdata$Dps
cdata <- cdata@data