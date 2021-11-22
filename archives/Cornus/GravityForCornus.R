rm(list=ls())
require(gstudio)
require(nlme)




######################################################################
#
#  				I just grabbed this from GRAV_FUNC_REV
#
######################################################################
gravity <- function(y, x, group, data, ln=TRUE, constrained=TRUE) {
  if (!require(nlme)) stop("nlme PACKAGE MISSING")
  if (!require(lattice)) stop("lattice PACKAGE MISSING")
  gdata <- data[,c(group,y,x)]
  
  if( ln==TRUE ) { 
    gdata[,x] <- log(abs(gdata[,x]))
    gdata[,y] <- log(abs(gdata[,y]))
    gdata[gdata == -Inf] <- 0 
    gdata[gdata == Inf] <- 0  		
  }
  
  fmla <- as.formula( paste( paste( y, "~", sep="") , paste( x, collapse= "+" ) ) )       	
  
  if (constrained == FALSE) { 
    print("Running unconstrained model, defaulting to OLS. Please check assumptions")
    gvlmm <- lm(fmla, data=gdata)
    gvaic <- AIC(gvlmm)
    return( list(	formula=fmla, 
                  GravityModel=gvlmm, 
                  Summary=summary(gvlmm), 
                  AIC=gvaic, 
                  x=gdata[,x], 
                  y=gdata[,y] ) )  	  	  
  } 
  else {
    gvlmm <- lme(fmla, random = ~1 | group, data=gdata)
    gvaic <- AIC(gvlmm)
    return( list(	formula=fmla, 
                  GravityModel=gvlmm, 
                  Summary=summary(gvlmm), 
                  AIC=gvaic, 
                  x=gdata[,x], 
                  y=gdata[,y], 
                  groups=gdata[,group] ) )  	  
  }
}	




##########################################################################
#
#	Population-level Dps (code from Murphy for population-level)
#
#	This is a shortcut for Dps among populations producing a matrix
#		
##########################################################################
braycurtis <- function( freq1, freq2){
  Dps <- 0
  K <- length( names( freq1 ) )
  
  for(i in 1:K){
    alleles <- unique( c(names(freq1[[i]]), names(freq2[[i]])) )
    a1 <- get.frequencies( freq1[[i]], alleles )
    a2 <- get.frequencies( freq2[[i]], alleles )
    Dps <- Dps + sum(apply( cbind(a1,a2), 1, min ) )
  }
  return( -1*log(Dps/K) )
}



find.r2 <- function( lme.object ){
  yhat <- predict( lme.object )
  obs <- lme.object$data[,1]
  f <- lm( obs ~ yhat )
  summary(f)$adj
}


lme.sum <- function(mod){
  res <- data.frame(
    "Beta.CI" = paste(round(summary(mod)$coefficients$fixed, 3), 
                      " (",
                      round(summary(mod)$coefficients$fixed-1.96*sqrt(diag(mod$varFix)), 2), 
                      ",", 
                      round(summary(mod)$coefficients$fixed+1.96*sqrt(diag(mod$varFix)), 2),
                      ")", 
                      sep=""),
    "P.value" = round(2 * pt(-abs(summary(mod)$coefficients$fixed/sqrt(diag(mod$varFix))), 
                             summary(mod)$fixDF[[1]]), 3)
  )
  return(res)
}


load("cornus.offs.rda")

co <- attr(cornus.off, "values")

df <- data.frame( ID = co$ID,
                  Mom = co$Mom, 
                  X = co$X,
                  Y = co$Y,
                  PctSky = co$PctSky,
                  Clumping = co$Clumping, 
                  DBH = co$DBH, 
                  FloralOutput = co$FloralOutput ) 

# Remake the loci in the new geneticStudio framework

for (locus in c("Cf.G8", "Cf.H18", "Cf.N10", "Cf.05") ) { 
  x <- c()
  loc <- co[[locus]]
  for( item in loc ) { 
    vals <- attr(item, "alleles")
    x <- c( x, locus( vals ) )
  }
  class( x ) <- "locus"
  df[[locus]] <- x 
}


summary( df )

#  Saturated Graph
df$Mom <- as.factor(df$Mom )
freqs <- lapply( partition( df, stratum="Mom"), frequencies )
coords <- strata_coordinates( df, stratum="Mom", longitude="X", latitude="Y")
sdist <- strata_distance( coords )
K <- length( levels( df$Mom ) )
N <- K*(K-1)
GDist <- dist_bray(df, stratum="Mom")
PDist <- sdist
Graph <- rep(FALSE,N)
DBH <- rep(NA,N)
Flowers <- rep(NA,N)
Sky <- rep(NA,N)
Clumping <- rep(NA,N)
momFrom <- rep(NA,N)
momTO <- rep(NA,N)
ID <- rep(NA,N)
momNums <- as.numeric(levels( df$Mom ))
offspring.graph <- population_graph(df,stratum="Mom")
graph.adj <- igraph::as_adjacency_matrix( offspring.graph )

ctr <- 1
for(i in 1:K){
  freq1 <- freqs[[i]]
  for(j in 1:K){
    if( i != j ) {
      freq2 <- freqs[[j]]
      if( graph.adj[i,j])
        Graph[ctr] <- TRUE
      momFrom[ctr] <- momNums[j]
      momTO[ctr] <- momNums[i]
      if( momNums[i] < momNums[j])
        ID[ctr] <-  paste(momNums[i],momNums[j], sep="-")
      else
        ID[ctr] <-  paste(momNums[j],momNums[i], sep="-")
      
      ctr <- ctr+1
    }
  }
}



# Fill out saturated graph
data <- data.frame( GDist=GDist[lower.tri( GDist)]/2, 
                    PDist=PDist[lower.tri( PDist)], 
                    Grouping=as.factor( rep( levels(df$Mom ), each=16) ),
                    FROM=as.factor(momFrom), TO=as.factor(momTO), ID, InGraph=Graph )
data$DBH <- log( rep( df$DBH[ df$ID==1 ], each=16 ))
data$PctSky <- log( rep( df$PctSky[ df$ID==1 ], each=16 ))
data$Clumping <- log( rep( df$Clumping[ df$ID==1 ], each=16 ))
data$FloralOutput <- log( rep( df$FloralOutput[ df$ID==1 ], each=16 ))
data$FloralOutput[ is.na(data$FloralOutput) ] <- 4.40564

# single terms saturated
fit.dist <- lme( GDist ~ PDist, random = ~ 1 | Grouping, data=data )
fit.dbh <- lme( GDist ~ DBH, random = ~ 1 | Grouping, data=data)
fit.sky <- lme( GDist ~ PctSky, random = ~ 1 | Grouping, data=data )
fit.clumping <- lme( GDist ~ Clumping, random = ~ 1 | Grouping, data=data)
fit.flowers <- lme( GDist ~ FloralOutput, random = ~ 1 | Grouping, data=data )
AIC.vals <- c(AIC(fit.dist), AIC(fit.dbh), AIC(fit.sky), AIC(fit.clumping), AIC(fit.flowers) )
R2.vals <- c( find.r2(fit.dist), find.r2(fit.dbh), find.r2(fit.sky), 
              find.r2(fit.clumping), find.r2(fit.flowers) )
names(AIC.vals) <- c("PDist","DBH","Sky","Clumping","Flowers")



# two terms saturated
fit.clump.dist <- lme( GDist ~ Clumping+PDist, random =~1|Grouping, data=data)
fit.clump.dbh <- lme( GDist ~ Clumping+DBH, random=~1|Grouping, data=data)
fit.clump.sky <- lme( GDist ~ Clumping+PctSky, random=~1|Grouping, data=data)
fit.clump.flowers <- lme( GDist ~ Clumping+FloralOutput, random=~1|Grouping, data=data)
AIC.t <- c(AIC(fit.clump.dist), AIC(fit.clump.dbh), AIC(fit.clump.sky), AIC(fit.clump.flowers) )
names(AIC.t) <- c("Clumping+PDist","Clumping+DBH","Clumping+Sky","Clumping+Flowers")
AIC.vals <- c(AIC.vals, AIC.t)
R2.vals <- c(R2.vals, find.r2(fit.clump.dist), find.r2(fit.clump.dbh), find.r2(fit.clump.sky), 
             find.r2(fit.clump.flowers))


gdata <- data[ Graph,]
# single terms saturated
fit.dist.gr <- lme( GDist ~ PDist, random = ~ 1 | Grouping, data=gdata)
fit.dbh.gr <- lme( GDist ~ DBH, random = ~ 1 | Grouping, data=gdata)
fit.sky.gr <- lme( GDist ~ PctSky, random = ~ 1 | Grouping, data=gdata )
fit.clumping.gr <- lme( GDist ~ Clumping, random = ~ 1 | Grouping, data=gdata)
fit.flowers.gr <- lme( GDist ~ FloralOutput, random = ~ 1 | Grouping, data=gdata )
AIC.t <- c(AIC(fit.dist.gr), AIC(fit.dbh.gr), AIC(fit.sky.gr), AIC(fit.clumping.gr), 
           AIC(fit.flowers.gr) )
names(AIC.t) <- c("graph.PDist","graph.DBH","graph.Sky","graph.Clumping","graph.Flowers")
AIC.vals <- c(AIC.vals, AIC.t)
R2.vals <- c(R2.vals,find.r2(fit.dist.gr),find.r2(fit.dbh.gr),find.r2(fit.sky.gr),
             find.r2(fit.clumping.gr),find.r2(fit.flowers.gr))

# two terms saturated
fit.clump.dist.gr <- lme( GDist ~ Clumping+PDist, random =~1|Grouping, data=gdata)
fit.clump.dbh.gr <- lme( GDist ~ Clumping+DBH, random=~1|Grouping, data=gdata)
fit.clump.sky.gr <- lme( GDist ~ Clumping+PctSky, random=~1|Grouping, data=gdata)
fit.clump.flowers.gr <- lme( GDist ~ Clumping+FloralOutput, random=~1|Grouping, data=gdata)
AIC.t <- c(AIC(fit.clump.dist.gr), AIC(fit.clump.dbh.gr), AIC(fit.clump.sky.gr), 
           AIC(fit.clump.flowers.gr) )
names(AIC.t) <- c("graph.Clumping+PDist","graph.Clumping+DBH","graph.Clumping+Sky",
                  "graph.Clumping+Flowers")
AIC.vals <- c(AIC.vals, AIC.t)
R2.vals <- c(R2.vals, find.r2(fit.clump.dist.gr), find.r2(fit.clump.dbh.gr),find.r2(fit.clump.sky.gr), 
             find.r2(fit.clump.flowers.gr))

# Summarize
dAIC.vals <- AIC.vals
dAIC.vals[1:9] <- AIC.vals[1:9] - min(AIC.vals[1:9])
dAIC.vals[10:18] <- dAIC.vals[10:18] - min(AIC.vals[10:18])

cat("AIC\n")
print(AIC.vals)

cat("dAIC\n")
print(dAIC.vals)

cat("R2\n")
names(R2.vals) <- names(AIC.vals)
print(R2.vals)


