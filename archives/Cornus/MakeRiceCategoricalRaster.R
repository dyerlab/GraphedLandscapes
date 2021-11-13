rm(list=ls())
require(raster)

ropen <- raster("open.asc")
rpine <- raster("pines.asc")
rcornus <- raster("cornus.asc")


o <- as.matrix( ropen )
o[ o == 0 ] <- NA
o[ !is.na(o)] <- 1

p <- as.matrix( rpine )
p[ p==0] <- NA
p[ !is.na(p) ] <- 2

c <- as.matrix( rcornus )
c[ c==0 ] <- NA
c[ !is.na(c) ] <- 3


rice <- matrix( 0 , nrow = nrow(o), ncol=ncol(o) )
rice[ !is.na(p) ] <- 1
rice <- raster( rice )


plot(rice)
r <- list()
r[[1]] <- rice
r[[2]] <- aggregate(rice)
r[[3]] <- aggregate( r2, fac=4 )
r[[4]] <- aggregate( r3, fac=4 )
r[[5]] <- aggregate( r4, fac=4 )
r[[6]] <- aggregate( r5, fac=4 )
r[[7]] <- aggregate( r6 )

par(mfrow=c(2,3))
for( i in 7:2 ){
  fname <- paste("rice",i,"png", sep=".")
  plot( r[[i]],legend=FALSE, xlab="X-Coordinate", ylab="Y-Coordinate" )
}



rice <- aggregate(rice, fac=4)
rice[ rice > 0 ] <- 1

rice <- ratify( rice )
rat <- levels(rice)[[1]]
rat
rat$landcover <- c("Canopy","Open")
levels(rice) <- rat
rice

plot(rice)





