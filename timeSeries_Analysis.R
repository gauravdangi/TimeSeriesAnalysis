library(plyr);library(dplyr)
library(forecast)
library(astsa)
library(tseries)
data("AirPassengers")
AP<-AirPassengers;rm(AirPassengers)
class(AP)  # "ts"  <- time series object
AP
#'      Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
#' 1949 112 118 132 129 121 135 148 148 136 119 104 118
#' 1950 115 126 141 135 125 149 170 170 158 133 114 140
#' 1951 145 150 178 163 172 178 199 199 184 162 146 166
#' 1952 171 180 193 181 183 218 230 242 209 191 172 194

summary(AP)
start(AP)  # [1] 1949    1
end(AP)    # [1] 1960   12
frequency(AP)  # 12
cycle(AP)
plot(AP);abline(reg = lm(AP~time(AP)))
aggregate(AP)
aggregate(AP,FUN = mean)
plot(aggregate(AP))
plot(aggregate(AP,FUN = mean))
boxplot(AP)
boxplot(AP~cycle(AP))

# Exponential smoothing state space model
fit<-ets(AP)  #  Returns ets model applied to y.
fit
plot(forecast(fit))
plot(fit)

acf(AP)
acf2(AP)   # astsa package

plot(diff(AP));abline(reg=lm(diff(AP)~time(diff(AP))))
plot(diff(diff(AP)));abline(reg=lm(diff(diff(AP))~time(diff(diff(AP)))))
plot(diff(log(AP)));abline(reg=lm(diff(log(AP))~time(diff(log(AP)))))
acf2(diff(AP))
acf2(diff(diff(AP)))


#HoltWinters computes Holt-Winters Filtering of a given time series. 
#Unknown parameters are determined by minimizing the squared prediction error.

plot(HoltWinters(AP))
#it fits holt-winter therefore we can also predict/forcast using holt-winter 
# As AP is good seasonal data, therefore we can predict using holt-winter way

AP.pred <- predict(HoltWinters(AP),n.ahead=10*12)
# n.ahead=10*12 means predicting for next 10 years

AP.pred
plot(AP.pred)
ts.plot(AP,AP.pred,lty=1:2)
ts.plot(AP,AP.pred,lty=1:2,type = "o",col=c(1,2))

mean(AP);mean(log(AP));mean(diff(AP));mean(diff(diff(AP)));mean(diff(log(AP)))
#[1] 280.2986
#[1] 5.542176
#[1] 2.237762
#[1] 0.2535211
#[1] 0.009440047

ts.plot(AP,diff(AP),diff(diff(AP)),col=c("blue","black","red"))
legend('topleft',c("AP","diff(AP)","diff(diff(AP))"),fill = c("blue","black","red"),bty='n')
   
decompose(AP)
plot(decompose(AP))
plot(decompose(diff(AP)))
                                                             
# Dicky fuller test (tseries is required, for adf.test)
adf.test(AP,alternative = "stationary",k=12)
adf.test(diff(AP),alternative = "stationary",k=12)
adf.test(diff(diff(AP)),alternative = "stationary",k=12)

library(fUnitRoots)
schwert.param <- trunc(12 * (length(na.omit(AP)) / 100) ^ (1 / 4))
adfTest(x = na.omit(AP), lags = schwert.param, type = "nc", title = NULL, description = NULL)


# ARIMA 

ARIMAfit <- auto.arima(AP,d=1,approximation=FALSE,trace=FALSE)

plot(AP,col="blue")
lines(fitted(ARIMAfit),col="red")
AP.pred <- predict(ARIMAfit,n.ahead = 10*12)
AP.pred<-as.ts(AP.pred)
plot(AP.pred$pred)
ts.plot(AP,AP.pred$pred,lty=1:2)

ddap<-auto.arima(diff(diff(AP)))
plot(diff(diff(AP)),col="blue")
lines(fitted(ddap),col="red")
legend('topleft',c("diff(diff(AP))","fitted"),fill = c("blue","red"),bty = "n")


#----------------------------------------------
#EXPONENTIAL SMOOTHING
#Simple Exponential Smoothing

##1.using ses() for forecast
library(fpp)
data(oil)

oildata <- window(oil,start=1996,end=2007)
plot(oildata, ylab="Oil (millions of tonnes)",xlab="Year")

fit1 <- ses(oildata, alpha=0.2, initial="simple", h=3)

plot(fit1, plot.conf=FALSE, ylab="Oil (millions of tonnes)", xlab="Year", main="", fcol="white", type="o")

lines(fitted(fit1), col="blue", type="o")
lines(fit1$mean, col="blue", type="o")

fit2 <- ses(oildata, alpha=0.6, initial="simple", h=3)
lines(fitted(fit2), col="red", type="o")
lines(fit2$mean, col="red", type="o")
forecast(fit2,h=3)
fitted(fit2)

fitted(fit2)[length(fitted(fit2))]

0.2*oildata[length(oildata)]+0.8*fitted(fit1)[length(fitted(fit1))] #=fit1$mean=484.8025
0.6*oildata[length(oildata)]+0.4*fitted(fit2)[length(fitted(fit2))] #=fit2$mean=501.8375

?ses

plot(oildata)
acf2(oildata)
plot(diff(oildata))

adf.test(oildata)  # p-value = 0.5566
adf.test(diff(oildata))  # p-value = 0.5566
adf.test(log(oildata))  # p-value = 0.5654

fit<-auto.arima(oildata,d=1)
ts.plot(oildata,fitted(fit),lty=1:2)
forecast(fit,h=3)

