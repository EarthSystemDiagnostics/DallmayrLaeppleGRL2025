## Aim: Fit the statistical model between isotopes and topography and make the corresponding figures
## Status/belongs to: Topography hypothesis manuscript
#d18O vs altitude
#Full dataset vs. dataset restricted to locations where the slope is well defined (e.g. slope.sd<0.4)
#d18O vs. slope; 
#Deming regression with the uncertainties given by the replicability in d18O (d18O.sde) and the slope uncertainty (slope.sd)
#return the models  (slope, intercept) and the standard error of the slope

## Result: 
## Next step: 
basedrive="/Users/tlaepple/data/" 
library(dplyr)
library(ggplot2)
library(sf)
library(stars)
library(deming)
setwd(paste(basedrive,"topohypothesis",sep=""))
load(file="./save/data.dat",verbose=TRUE)



summary(lm(d18O~altitude,data=data)) #-3.39 +- 0.934 permille per meter, r2 = 0.243


ylab.d18O <- expression(delta^{18}*O~"["*"\u2030"*"]")
SC=2

pdf(file="./plots/Fig2c.pdf",width=3*SC,height=3*SC)
par(mai=c(1,1.1,0.8,0.4))
# simply Plot the d18O vs altitude,
plot(data$altitude, data$d18O, col = "black", pch = 19, cex = 1, 
     xlab = "Elevation [m]", ylab = ylab.d18O,xlim=range(data$altitude)+c(0,0),ylim=range(data$d18O)+c(-0.7,0.7))
model<-lm(d18O~altitude,data=data)
abline(model,lwd=2,col="grey")

arrows(x0 = data$altitude, y0 = data$d18O - 2*data$d18O.sde, 
       x1 = data$altitude, y1 = data$d18O + 2*data$d18O.sde, 
       length = 0.05, angle = 90, code = 3, col = "grey")

points(data$altitude, data$d18O, col = "black", pch = 19, cex = 1)

dev.off()



data.filtered <- data %>% filter(slope.sd < 0.4)

model.deming <- deming::deming( d18O~slope, data=data,xstd=c(slope.sd),ystd=d18O.sde) #0.38 +- 0.28; #0.47 +- 0.15.  0.18-0.75
model.deming.filtered <- deming::deming( d18O~slope, data=data.filtered,xstd=c(slope.sd),ystd=d18O.sde)# 0.42 +- 0.09 #0.42 +- 0.08. 0.27-0.60


cor.test(data.filtered$slope,data.filtered$d18O) #p=1e-5, df=28, N=30, R^2=0.504
cor.test(data$slope,data$d18O) #p=0.006, df=41, N=43 R^2=0.16
cor.test(data$slope,data$altitude)

summary(lm(data.filtered$d18O~data.filtered$slope)) #0.36 


###

pdf(file="./plots/Fig2d.pdf",width=3*SC,height=3*SC)
par(mai=c(1,1,1,1))
par(mai=c(1,1.1,0.8,0.4))
par(mfcol=c(1,1))


plot(data$slope, data$d18O, ylim=c(-2,2)-44.3, xlim=c(-2,5), 
     xlab="Slope [m/km]", ylab = ylab.d18O, 
     main="", pch=19, col="grey", cex=1)

# Add error bars on both sides
arrows(x0 = c(data$slope) - 2*data$slope.sd, y0 = data$d18O, 
       x1 = c(data$slope) + 2*data$slope.sd, y1 = data$d18O, 
       length = 0.05, angle = 90, code = 3, col = "grey")

arrows(x0 = data$slope, y0 = data$d18O - 2*data$d18O.sde, 
       x1 = data$slope, y1 = data$d18O + 2*data$d18O.sde, 
       length = 0.05, angle = 90, code = 3, col = "grey")


  abline(model.deming$coefficient[1], model.deming$coefficient[2], col="black", lwd=2)
  abline(model.deming.filtered$coefficient[1], model.deming.filtered$coefficient[2], col="black", lwd=2, lty=2)

# Add the filtered data points with error bars in red
points(data.filtered$slope, data.filtered$d18O, col="black", cex=1,pch=19)

dev.off()
