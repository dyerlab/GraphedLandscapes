rm(list=ls())
require(raster)
require(rgeos)

# load the data file
load("cornus.offs.rda")


#####################################################
#                                                   #
#       Find some coordinates for the mothers       #
#                                                   #
#####################################################
#r is not recognizing the $ operator in the gstudio population object, so here is a workaround
Mom <- cornus.off[,2]
X <- cornus.off[,7]
Y <- cornus.off[,8]

tmp <- unique( cbind( Mom, X, Y ) )
coords <- data.frame(ID=as.character(tmp[,1]), X=tmp[,2], Y=tmp[,3] )
coords.sp <- SpatialPointsDataFrame( data=coords, coords=coords[,2:3])

# make a data frame
K <- dim(coords)[1]
N <- K*(K-1)/2

# load the raster into memory
open.raster <- raster( "open.asc" )
cornus.raster <- raster( "cornus.asc" )

roads.raster <- raster( "roads.asc" )
pines.raster <- raster( "pines.asc" )


raster.data <- data.frame( ID=character(N), 
                           open.mn=numeric(N),    open.var=numeric(N),
                           cornus.mn=numeric(N),  cornus.var=numeric(N),
                           roads.mn=numeric(N),   roads.var=numeric(N),
                           pines.mn=numeric(N),   pines.var=numeric(N),
                           decid.mn=numeric(N),   decid.var=numeric(N),
                           stringsAsFactors=FALSE )

# Go through the coordinates and find pairwise values
ctr <- 1
for(i in 1:K){
  for(j in i:K){
    if(j>i){
      line <- Line( rbind(coords[i,2:3],coords[j,2:3]) )
      spLines <- SpatialLines( list(Lines( line, ID="Lines")) )
      
      #generate bandwith around spLine. width = bandwith size
      bufLines <- gBuffer(spLines,width=1)
      
      vals.cornus <- extract(cornus.raster, bufLines )[[1]]
      
      vals.roads <- extract(roads.raster, bufLines )[[1]]
      
      vals.open <- extract(open.raster,bufLines )[[1]]
      
      vals.pines <- extract(pines.raster, bufLines )[[1]]
      
      vals.decid <- vals.open + vals.pines
      vals.decid[ vals.decid > 1.0 ] <- 1.0
      vals.decid <- 1.0 - vals.decid
      
      raster.data$ID[ctr] <- paste( coords$ID[i], coords$ID[j], sep="-")
      
      raster.data$open.mn[ctr] <- mean( vals.open )
      raster.data$open.var[ctr] <- var( vals.open )
      
      raster.data$cornus.mn[ctr] <- mean( vals.cornus )
      raster.data$cornus.var[ctr] <- var( vals.cornus )
	  raster.data$cornus.mn[ctr] <- mean( vals.cornus )
      raster.data$cornus.var[ctr] <- var( vals.cornus )
      
      raster.data$roads.mn[ctr] <- mean( vals.roads )
      raster.data$roads.var[ctr] <- var( vals.roads )
      
      raster.data$pines.mn[ctr] <- mean( vals.pines )
      raster.data$pines.var[ctr] <- var( vals.pines )
      
      raster.data$decid.mn[ctr] <- mean( vals.decid )
      raster.data$decid.var[ctr] <- var( vals.decid )
      
      #print(raster.data[ctr,1:11])
      
      ctr <- ctr + 1

    }
  }
}

