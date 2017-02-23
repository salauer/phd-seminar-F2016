## Flu data investigation by Stephen Lauer and Sara Nunez
## See: https://predict.phiresearchlab.org/legacy/flu/targets.html and https://github.com/cdcepi/flu-tools-2015-16

library(ggplot2)
library(readr)
library(astsa)
library(fGarch)
library(forecast)
library(dplyr)

train_series <- read_csv("data/train_data.csv", na="X") %>% 
    filter(YEAR>=2003) %>% 
    select(time=X1, ILI=X.UNWEIGHTED.ILI)

test_series <- read_csv("data/test_data.csv") %>% 
    select(time=X1, ILI=X.UNWEIGHTED.ILI)

## Fit and predict from a SARIMA model
ggplot() +
    geom_line(aes(x=train_series$time, y=train_series$ILI))
acf2(train_series$ILI,200)
ggplot() +
    geom_line(aes(x=train_series$time[-1], y=diff(train_series$ILI)))
acf2(diff(train_series$ILI),200)
ggplot() +
    geom_line(aes(x=train_series$time[-1], y=diff(log(train_series$ILI))))
acf2(diff(log(train_series$ILI)),200)
ggplot() +
    geom_line(aes(x=1:456, y=diff(diff(log(train_series$ILI),52))))
acf2(diff(diff(log(train_series$ILI),52)),200)
acf2(diff(diff(train_series$ILI),52),200)

dILI <- as.ts(diff(log(train_series$ILI)))

ILI_auto_arima <- auto.arima(dILI, stepwise = FALSE, approximation = FALSE,seasonal = TRUE, trace=TRUE)
ILI_auto_arima$arma
sarima(dILI, 1,1,0,0,1,0,0) ## don't know why the period is 0 here
sarima(dILI, 1,1,0,0,1,0,52) ## setting period to 52 messes everything up
ILI_sarima <- arima(dILI,
                    order=c(1,1,0),
                    seasonal=list(order=c(0,1,0),
                                  period=0))

testILI <- as.ts(diff(log(test_series$ILI)))
ILI_preds <- Arima(testILI, model=ILI_sarima)

ggplot() +
    geom_line(aes(x=test_series$time[-1], y=testILI)) +
    geom_line(aes(x=test_series$time[-1], y=testILI+ILI_preds$residuals), color="red", linetype="dashed")

mean(abs(ILI_preds$residuals))

## Build a GARCH Model, based on fit from SARIMA; try for different garch(m,r) values
ILI_garch10 <- garchFit(~arma(2,0)+garch(1,0), dILI)
summary(ILI_garch10)
ILI_garch20 <- garchFit(~arma(2,0)+garch(2,0), dILI)
summary(ILI_garch20)
ILI_garch11 <- garchFit(~arma(2,0)+garch(1,1), dILI)
summary(ILI_garch11)
ILI_garch21 <- garchFit(~arma(2,0)+garch(2,1), dILI)
summary(ILI_garch21)
ILI_garch10@fit$ics
ILI_garch20@fit$ics
ILI_garch11@fit$ics
ILI_garch21@fit$ics
ggplot() +
    geom_line(aes(x=train_series$time[-1],y=dILI)) +
    geom_line(aes(x=train_series$time[-1],y=dILI+ILI_garch11@residuals), color="red", linetype="dashed")

whole_series <- bind_rows(train_series, test_series)
garch_preds <- c()
for(i in 0:(nrow(test_series)-2)){
    new_train <- whole_series[1:(nrow(train_series)+i),]
    new_dILI <- as.ts(diff(log(new_train$ILI)))
    new_fit <- garchFit(~arma(2,0)+garch(1,1), new_dILI, trace=F)
    garch_preds[i+1] <- predict(new_fit, n.ahead=1)$meanForecast
}
ggplot() +
    geom_line(aes(x=1:208, y=testILI)) +
    geom_line(aes(x=1:208, y=garch_preds), color="red", linetype="dashed")

mean(abs(testILI-garch_preds))
