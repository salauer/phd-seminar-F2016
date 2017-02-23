## Flu data investigation by Stephen Lauer and Sara Nunez
## See: https://predict.phiresearchlab.org/legacy/flu/targets.html and https://github.com/cdcepi/flu-tools-2015-16

library(ggplot2)
library(readr)
library(astsa)
library(dplyr)


train_series <- read_csv("data/train_data.csv", na="X") %>% 
    filter(YEAR>=2003) %>% 
    select(time=X1, ILI=X.UNWEIGHTED.ILI)

qplot(data=train_series, x=time, y=ILI, geom="line")
# qplot(data=train_series, x=ILI)
# qplot(data=train_series, x=log(ILI))

usflu <- read_csv("data/train_data.csv", na="X")

qplot(data=usflu, x=YEAR+(WEEK-1)/53,y=X.UNWEIGHTED.ILI, geom="line")

usflu2 <- usflu[,c("X1","X.UNWEIGHTED.ILI")]
usflu3 <- usflu2[usflu$YEAR>=2003,]

## ORIGINAL
par(mfrow=c(2,1))
acf2(usflu3$X.UNWEIGHTED.ILI,max.lag = 53)
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
qplot(x=usflu3$X1[-1], y=dlogflu, geom="line")
acf(dlogflu,lag.max= 53)
acf(dlogflu,lag.max= 53,type=c("partial"))

ddlogflu <- diff(diff(log(usflu3$X.UNWEIGHTED.ILI)),52)
acf(ddlogflu,lag.max= 53)
acf(ddlogflu,lag.max= 53,type=c("partial"))

##### with some logging stuff
## ORIGINAL
acf2(train_series$ILI, 200)
acf2(train_series$log_ILI, 200)
acf2(diff(train_series$ILI), 200)
acf2(diff(train_series$log_ILI), 200)

## FIRST DIFFERENCE QPLOT
qplot(x=train_series$time[-1],y=diff(train_series$ILI), geom="line")
qplot(x=train_series$time[-1],y=diff(log(train_series$ILI)), geom="line")

## SEASONAL DIFFERENCING
acf2(diff(diff(train_series$ILI),52), 200)
acf2(diff(diff(log(train_series$ILI)),52), 200)
qplot(y=diff(diff(train_series$ILI),52)) + geom_line()
