##########################################################################
# Gravity model lab exercise- Columbia spotted frog
##########################################################################
# CONTACT:
#     Melanie Murphy
#     University of Wyoming 
#     melanie.murphy@uwyo.edu
#     (208_874-3749
##########################################################################

#install new packages
install.packages("rgdal")
install.packages("spatialEco")
install.packages("GeNetIt")

#require packages (may also need to install some of the below #packages)
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
( xvars <- stack(rasters[-6]) )
( land.cov <- stack(rasters[6]) )

#plot sample locations
plot(ralu.site)

#plot raster data
plot(xvars)


###################################################################  
# Build at site covariates by extracting raster point values 
###################################################################  

#### Build at site covariates 
ralu.site@data <- data.frame(ralu.site@data, extract(xvars[[c(1,2)]], ralu.site))
#at each sample site, extracting values from the first and second rasters

#Take a look at the resulting data
head(ralu.site@data)  
# 1) Create kNN graph from the site data (ralu.site)
dist.graph <- knn.graph(ralu.site, row.names = ralu.site@data[,"SiteName"])

#you can limit the graph by maximum distance #(max.dist) if desired
#dist.graph <- knn.graph(ralu.site, row.names = ralu.site@data[,"SiteName"], max.dist = 5000)

# 2) Add "from.to" unique ID's and merge with genetic distance matrix 
  dist.graph@data$from.to <- paste(dist.graph$i, dist.graph$j, sep=".")
  dps$from.to <- paste(dps$FROM_SITE, dps$TO_SITE, sep=".") 
  dist.graph <- merge(dist.graph, dps, by = "from.to") 

# 3) Merge graph with at site nodes 
  dist.graph@data <- dist.graph@data[,-c(7,8)]

# Remove NA values
  na.index <- unique(as.data.frame(which(is.na(dist.graph@data), arr.ind = TRUE))[,1])
  dist.graph <- dist.graph[-na.index, ]

# Display columns and plot
str(dist.graph@data)  
plot(xvars[[2]])
 plot(dist.graph, add=T)
  points(ralu.site, pch=20, col="red")


###################################################################  
# Build between site covariates by extracting raster values 
#   and calculating statistics
###################################################################  
#can calculate any statistical moment to describe the values between sites
#This examples has min, mean, max and variance
#Example is all for floating point data

stats <- graph.statistics(dist.graph, r = xvars, d=30, 
            stats = c("min", "mean", "max", "var"),
            sp = FALSE) 
dist.graph@data <- data.frame(dist.graph@data, stats)

##Calculating statistical moments of categorical data does not make sense
#### Example of categorical raster data
####   function for percent wetland landcover (NLCD).

wet.pct <- function(x) { 
  x <- ifelse( x == 11 |  x == 12 | x == 90 | x == 95, 1, 0)
#have multiple numeric categories that are wetlands

   prop.table(table(x))[2] 
}   
lc.stats <- graph.statistics(dist.graph, r = land.cov, d=30, 
            stats = "wet.pct")
    lc.stats[is.na(lc.stats)] <- 0
dist.graph@data <- data.frame(dist.graph@data, lc.stats)

                
str(dist.graph@data)


###################################################################  
# Build data for gravity model 
###################################################################

from <- ralu.site@data[,c(1,6,8,18,19)]
  names(from)[2:ncol(from)] <- paste("f", names(from)[2:ncol(from)], sep=".") 
to <- ralu.site@data[,c(1,6,8,18,19)]
  names(to)[2:ncol(to)] <- paste("t", names(to)[2:ncol(to)], sep=".") 
site <- data.frame(from,to)  
  site <- site[,-(dim(to)[2]+1)]

cdata <- merge(dist.graph, site, by.x="from_ID", by.y="SiteName") 
cdata$Dps <- 1 - cdata$Dps
cdata <- cdata@data

###################################################################  
# Specify and fit gravity model 
###################################################################
# Specify parameters
#cti ? compound topographic index (wetness index)
#err27 ? elevation relief ratio (measure of topographic complexity, 27X27 window)
#gsp ? growing season precipitation
#f.AREA_m2 ? wetland area of the ?from? site in meters squared
#f.Depth_m ? wetland depth of the ?from? site, this is highly correlated with predatory fish
#f.err27 ? elevation relief ratio at the wetland
#length is distance between sites in meters

x = c("length", "mean.cti", "mean.err27", "mean.ffp", "mean.gsp",  
      "f.AREA_m2", "f.Depth_m", "f.err27")

# Gravity model
#group variable ? this sets the constraint.  This is a from node constraint (singly constrained #production model)

( gm <- gravity(y = "Dps", x = x, d = "length", group = "from_ID", data = cdata) )
# Plot gravity results
par(mfrow=c(2,3))
   for (i in 1:6) { plot(gm, type=i) } 

# Try a second model.  
x = c("length", "mean.cti", "mean.err27", "mean.ffp", ?wet.pct.nlcd?,   "f.Depth_m", "f.cti")
( gm.2 <- gravity(y = "Dps", x = x, d = "length", group = "from_ID", data = cdata) )
par(mfrow=c(2,3))
   for (i in 1:6) { plot(gm, type=i) } 



###################################################################
