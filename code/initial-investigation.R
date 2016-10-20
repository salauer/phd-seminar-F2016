## Flu data investigation by Stephen Lauer and Sara Nunez
## See: https://predict.phiresearchlab.org/legacy/flu/targets.html and https://github.com/cdcepi/flu-tools-2015-16

library(cdcfluview)
library(ggplot2)
library(readr)
library(astsa)

#usflu <- get_flu_data("national", "ilinet", years=1997:2010)
#write.csv(usflu, "data/usflu.csv")
usflu <- read_csv("data/usflu.csv", na="X")

qplot(data=usflu, x=YEAR+(WEEK-1)/53,y=X.UNWEIGHTED.ILI, geom="line")

usflu2 <- usflu[,c("X1","X.UNWEIGHTED.ILI")]
usflu3 <- usflu2[usflu$YEAR>=2003,]

## ORIGINAL
par(mfrow=c(2,1))
acf(usflu3$X.UNWEIGHTED.ILI,lag.max = 53)
acf(usflu3$X.UNWEIGHTED.ILI,lag.max = 53,type="partial")

## FIRST DIFFERENCE
dflu <- diff(usflu3$X.UNWEIGHTED.ILI)
acf(dflu,lag.max= 53)
acf(dflu,lag.max= 53,type=c("partial"))

## FIRST DIFFERENCE QPLOT
qplot(x=usflu3$X1[-c(1)],y=dflu, geom="line")

## SEASONAL DIFFERENCING
ddflu <- diff(diff(usflu3$X.UNWEIGHTED.ILI),52)
acf(ddflu,lag.max = 53)
acf(ddflu,lag.max = 53,type="partial")

## ALL PLOTS
par(mfrow=c(2,2))
acf(dflu,lag.max= 53)
acf(dflu,lag.max= 53,type=c("partial"))
acf(ddflu,lag.max = 53)
acf(ddflu,lag.max = 53,type="partial")

par(mfrow=c(2,1))
dlogflu <- diff(log(usflu3$X.UNWEIGHTED.ILI))
acf(dlogflu,lag.max= 53)
acf(dlogflu,lag.max= 53,type=c("partial"))

ddlogflu <- diff(diff(log(usflu3$X.UNWEIGHTED.ILI)),52)
acf(ddlogflu,lag.max= 53)
acf(ddlogflu,lag.max= 53,type=c("partial"))
