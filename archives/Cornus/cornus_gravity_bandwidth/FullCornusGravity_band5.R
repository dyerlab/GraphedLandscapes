##########################################################################
#
#  This is the whole enchelada for the Cornus stuff
#
##########################################################################



##########################################################################
# PROGRAM: GravityModel
# USE: NETWORK GRAVITY MODEL
# REQUIRES: DATAFRAME GROUPED BY NODE ID FACTOR
#
#  gravity(y, x, group, ln=TRUE)
#
# ARGUMENTS:  
#    y               DEPENDENT VARIABLE  
#    x               LIST OF INDEPENDENT VARIABLES
#    group           NAME OF GROUPING COLUMN
#    ln              NATURAL LOG TRANSFORM DATA (DEFAULT TRUE)
#    constrained     RUN CONSTRAINED MODEL (DEFAULT TRUE). IF FALSE A LINEAR MODEL IS RUN
#
# RETURNS: 
#  LIST OBJECT CONTANING: 
#    formula         GARVITY MODEL FORMULA  
#    GravityModel    RESULTING MODEL
#    Summary         SUMMARY OF MODEL
#    AIC             AIC VALUE FOR MODEL
#    x               INDEPENDENT VARIABLE(S)
#    y               DEPENDENT VARIABLE
#    groups          GROUPING VARIABLE
#
# EXAMPLE:   
#  cdata <- read.csv("C:/ANALYSIS/GenNet/mdata.csv")
#    cdata$group <- factor(cdata$FROM_SITE)  
#    cdata <- groupedData( DPS ~ DISTANCE | group, data=cdata)
#  y="DPS" 
#  x=c("DISTANCE", "DEPTH_F", "HLI_F", "CTI_F", "cti", "ffp")  
#  gm <- gravity(y=y, x=x, group=group, ln=TRUE)
#   
# REFERENCES:
#    Murphy M.A., R. Dezzani, D.S. Pilliod & A.S. Storfer (2010) Landscape genetics of  
#       high mountain frog metapopulations. Molecular Ecology 19(17):3634-3649 
#
#     Murphy, M. A., J.S. Evans. (in prep) GenNetIT: gravity analysis in R for landscape 
#       genetics. Molecular Ecology Resources
#
# CONTACT: 
#     Jeffrey S. Evans 
#     Senior Landscape Ecologist 
#     The Nature Conservancy - Central Science
#     Adjunct Faculty
#     University of Wyoming
#     Laramie, WY
#     (970)672-6766
#     jeffrey_evans@tnc.org
##########################################################################
gravity <- function(y, x, group, data, ln=TRUE, constrained=TRUE) {
  if (!require(nlme)) stop("nlme PACKAGE MISSING")
  if (!require(lattice)) stop("lattice PACKAGE MISSING")
  gdata <- data[,c(group,y,x)] 

  if(ln==TRUE) { 
    gdata[,x] <- log(abs(gdata[,x]))
    gdata[,y] <- log(abs(gdata[,y]))
    gdata[gdata == -Inf] <- 0 
    gdata[gdata == Inf] <- 0  		
  } 		   		
  fmla <- as.formula(paste(paste(y, "~", sep="") , paste(x, collapse= "+")))       
  if (constrained == FALSE) { 
    print("Running unconstrained model, defaulting to OLS. Please check assumptions")
    gvlmm <- lm(fmla, data=gdata)
    gvaic <- AIC(gvlmm)
    return( list(formula=fmla, GravityModel=gvlmm, Summary=summary(gvlmm), 
                 AIC=gvaic, x=gdata[,x], y=gdata[,y]) )  	  	  
  } else {
    gvlmm <- lme(fmla, random = ~1 | group, data=gdata)
    gvaic <- AIC(gvlmm)
    return( list(formula=fmla, GravityModel=gvlmm, Summary=summary(gvlmm), 
                 AIC=gvaic, x=gdata[,x], y=gdata[,y], groups=gdata[,group]) )  	  
  }
}	

##########################################################################
# PROGRAM: plot.gm
# USE: FUNCTION TO PLOT GRAVITY MODEL OBJECTS
# REQUIRES: OBJECT OF CLASS gm
#
#  plot.gm(x, type=1)
#
# ARGUMENTS:  
#    x         OBJECT OF CLASS gm 
#    type      TYPE OF PLOT
#           1 - MODEL STRUCTURE I
#           2 - MODEL STRUCTURE II
#           3 - Q-Q NORMAL - ORIGIN RANDOM EFFECTS
#           4 - Q-Q NORMAL - RESIDUALS
#           5 - FITTED VALUES	
#           6 - DISTRIBUTION OF OBSERVED VS. PRED		
#
# RETURNS: 
#    SPECIFIED PLOT
#
# EXAMPLE:
#  # CREATES PDF OF ALL MODEL PLOTS
#   pdf("GravityModel_Plots.pdf")
#     for (i in 1:6) { plot.gm(gm, type=i) }
#   dev.off() 
#   
# REFERENCES:
#    Murphy M.A., R. Dezzani, D.S. Pilliod & A.S. Storfer (2010) Landscape genetics of  
#       high mountain frog metapopulations. Molecular Ecology 19(17):3634-3649
#
#     Murphy, M. A., J.S. Evans. (in prep). "GenNetIT: gravity analysis in R for landscape 
#       genetics" Molecular Ecology Resources
#
# CONTACT: 
#     Jeffrey S. Evans 
#     Senior Landscape Ecologist 
#     The Nature Conservancy - Central Science
#     Adjunct Faculty
#     University of Wyoming
#     Laramie, WY
#     (970)672-6766
#     jeffrey_evans@tnc.org
##########################################################################
plot.gm <- function(x, type=1) {
  if(type == 1) {
    # MODEL STRUCTURE I  
    plot(fitted(x$GravityModel, level=0), x$y, xlab = "Fitted Values (DPS)",
         ylab="Observed Values", main="Model Structure (I)", pch=16)
    abline(0, 1, col = "blue")
  }
  if(type == 2) {
    # MODEL STRUCTURE II
    scatter.smooth(fitted(x$GravityModel), residuals(x$GravityModel, type="pearson"), 
                   ylab="Innermost Residuals", main="Model Structure (II)",
                   xlab="Fitted Values", pch=16)
    abline(h = 0, col = "red")
  }
  if(type == 3) {			  
    # Q-Q NORMAL - ORIGIN RANDOM EFFECTS		   
    qqnorm(ranef(x$GravityModel)[[1]], main="Q-Q Normal - Origin Random Effects", pch=16)
    qqline(ranef(x$GravityModel)[[1]], col="red")
  }
  if(type == 4) {
    # Q-Q NORMAL - RESIDUALS		   
    qqnorm(residuals(x$GravityModel, type="pearson"), main="Q-Q Normal - Residuals", pch=16)
    qqline(residuals(x$GravityModel, type="pearson"), col="red")
  }  
  if(type == 5) {
    # FITTED VALUES	
    options(warn=-1)	   
    boxplot(residuals(x$GravityModel, type="pearson", level=1) ~ x$groups,
            ylab="Innermost Residuals", xlab="Origin",
            notch=T, varwidth = T, at=rank(ranef(x$GravityModel)[[1]]))
    axis(3, labels=format(ranef(x$GravityModel)[[1]], dig=2), cex.axis=0.8,
         at=rank(ranef(x$GravityModel)[[1]]))
    abline(h=0, col="darkgreen")  
  }  
  if(type == 6) {
    # DISTRIBUTION OF OBSERVED VS. PRED		
    oden <- density(x$y)
    pden <- density(predict(x$GravityModel)) 
    plot(oden, type="n", main="", xlim=c(min(x$y), max(x$y)), 
         ylim=c(min(oden$y,pden$y), max(oden$y,pden$y)))     
    polygon(oden, col=rgb(1,0,0,0.5))
    polygon(pden, col=rgb(0,0,1,0.5))
    legend("topright", legend=c("Obs","Pred"), fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))  )
  }
}



##########################################################################
#
#       A little function to print out a nice lme object
#
##########################################################################
lme.sum <- function(mod){
  res <- data.frame(
    "Beta" = summary(mod)$coefficients$fixed,
    "Lower.CI" = summary(mod)$coefficients$fixed - 1.96*sqrt(diag(mod$varFix)),
    "Upper.CI" = summary(mod)$coefficients$fixed + 1.96*sqrt(diag(mod$varFix)),
    "P.value" = 2 * pt(-abs(summary(mod)$coefficients$fixed/sqrt(diag(mod$varFix))), 
                             summary(mod)$fixDF[[1]])
  )
  return(res)
}

find.r2 <- function( lme.object ){
  yhat <- predict( lme.object )
  obs <- lme.object$data$GDist
  f <- lm( obs ~ yhat )
  summary(f)$adj
}



require(nlme)
load("Cornus.data.band5.rda")

graph.data <- cornus.data.band5[ cornus.data.band5$InGraph, ]
graph.data$group <- as.factor(graph.data$Grouping)



# Backward Selection
gravity.1 <- gravity(y="GDist", x=c("DBH","PctSky","Clumping","FloralOutput","PDist",
                                    "open.var","cornus.mn","roads.mn","pines.var","decid.mn"),
                     group="group", data=graph.data)
print(gravity.1$formula)
print(lme.sum(gravity.1$GravityModel))

gravity.2 <- gravity(y="GDist", x=c("DBH","PctSky","Clumping","FloralOutput","PDist",
                                    "open.var","roads.mn","pines.var","decid.mn"),
                     group="group", data=graph.data)
print(gravity.2$formula)
print(lme.sum(gravity.2$GravityModel))

gravity.3 <- gravity(y="GDist", x=c("DBH","PctSky","Clumping","FloralOutput","PDist",
                                    "open.var","roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.3$formula)
print(lme.sum(gravity.3$GravityModel))

gravity.4 <- gravity(y="GDist", x=c("DBH","Clumping","FloralOutput","PDist",
                                    "open.var","roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.4$formula)
print(lme.sum(gravity.4$GravityModel))

gravity.5 <- gravity(y="GDist", x=c("DBH","FloralOutput","PDist",
                                    "open.var","roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.5$formula)
print(lme.sum(gravity.5$GravityModel))

gravity.6 <- gravity(y="GDist", x=c("DBH","PDist",
                                    "open.var","roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.6$formula)
print(lme.sum(gravity.6$GravityModel))

gravity.7 <- gravity(y="GDist", x=c("DBH",
                                    "open.var","roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.7$formula)
print(lme.sum(gravity.7$GravityModel))

gravity.8 <- gravity(y="GDist", x=c("DBH",
                                    "roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.8$formula)
print(lme.sum(gravity.8$GravityModel))

gravity.9 <- gravity(y="GDist", x=c(
                                    "roads.mn","decid.mn"),
                     group="group", data=graph.data)
print(gravity.9$formula)
print(lme.sum(gravity.9$GravityModel))

gravity.10 <- gravity(y="GDist", x=c("roads.mn"),
                     group="group", data=graph.data)
print(gravity.10$formula)
print(lme.sum(gravity.10$GravityModel))

#extract AIC
gravity.1$AIC
gravity.2$AIC
gravity.3$AIC
gravity.4$AIC
gravity.5$AIC
gravity.6$AIC
gravity.7$AIC
gravity.8$AIC
gravity.9$AIC
gravity.10$AIC