library(Mcomp)
library(ggplot2)
library(fpp)
library(forecast)
library(smooth)
library(stats)
library(tibble)
library(dplyr)
library(nlme)
library(forecTheta)

M3[[1357]]$x


# Part 1
# My student ID is 209572823

y_raw <- M3[[1357]]

autoplot(y_raw)
str(y_in)
summary(y_in)


# Create in-sample / out-sample subset
y_in <- y_raw$x
y_out <- y_raw$xx
y_out
autoplot(y_in,
    xlab = "Time", ylab = "Employment",
    main="Quarterly employment data of Finland from 1978 to 1996"
    )

summary(y_in)

# autoregression
gglagplot(y_in)

ggtsdisplay(y_in, main="Time Series Data, ACF, & PACF")

# Classical decomposition

# MA approach
length(y_in)
k4 <- 4
k8 <- 8
m4 <- k4 %/% 2
m8 <- k8 %/% 2

MA4 <- array(NA, length(y_in))
MA8 <- array(NA, length(y_in))

# MA4
for (i in (m4 + 1):(length(y_in) - m4)){

    MA4[i] = mean(c(y_in[(i-m4):(i+m4-1)], y_in[(i-m4+1):(i+m4)]))

}

# MA8
for (i in (m8 + 1):(length(y_in) - m8)){

    MA8[i] = mean(c(y_in[(i-m8):(i+m8-1)], y_in[(i-m8+1):(i+m8)]))

}

MA4 <- ts(MA4, start=c(1978,2), frequency=4)
MA8 <- ts(MA8, start=c(1978,2), frequency=4)

autoplot(y_in, series = "Data", main="Time Series vs MA(4) vs MA(8)") +
    autolayer(MA4, series="MA(4)") +
    autolayer(MA8, series="MA(8)") 

# Decomposition
dym <- decompose(y_in, type = "multiplicative")

autoplot(dym)

plot(dym)

checkresiduals(y_in)


nsdiffs(dy$seasonal)

adf.test(y_in)
kpss.test(y_in)





time(y_in)



# ETS 

# Compare ES models
# 1. Holt's Exponential Smoothing with linear trend
# 2. Damped Exponential Smoothing
# 3. Holt-Winters with additive seasonality
# 4. Holt-Winter's Exponential Smoothing with Damped Trend

no_ETS_models <- 4

models <- c("AAN", "AAN", "AAA", "MAM")
damping <- c(FALSE, TRUE, FALSE, TRUE)


ets(y_in, model="MAM", damped = TRUE)$par['phi']


MAPEs_ETS <- array(NA, no_ETS_models)
phi_ETS <- array(NA, no_ETS_models)
aic_ETS <- array(NA,no_ETS_models)

for (m in 1:no_ETS_models){
    
    fit <- ets(y_in, model=models[m], damped=damping[m])
    aic_ETS[m] <- fit$aic
    fcs <- forecast(fit, h=8)$mean
    MAPEs_ETS[m] <- accuracy(fcs, y_out)[5]
    phi_ETS[m] <- fit$par['phi']
    print(m)
}


ETS_Comparison <- tibble(model = models, damping = damping, 
    phi = phi_ETS, AIC =aic_ETS, MAPE = MAPEs_ETS)

ETS_Comparison

forecast(ets(y_in, model="AAN", damped=TRUE),h=8)
autoplot(forecast(ets(y_in, model="AAN", damped=TRUE),h=8))

# ARIMA


tsdisplay(y_in)
# No variation, no need for box-cox

Box.test(y_in, lag=1, type="Ljung")

tsdisplay(diff(y_in, 4))
tsdisplay(diff(diff(y_in, 4), 1))

diff_y_in <- diff(diff(y_in,4),1)

adf.test(diff_y_in)
tsdisplay(diff_y_in)

# OR ndiffs / nsdiffs
nsdiffs(y_in)
ndiffs(diff(y_in,4))

ndiffs(diff(diff(y_in,4),1))


ally <- cbind("Original data" = y_in,

              "Seasonal differences" = diff(y_in, 4),

              "First & seasonal diffs" = diff(diff(y_in, 4),1))

autoplot(ally, facets=TRUE) +

  xlab("Time") +

  ggtitle("Quarterly employment data of Finland from 1978 to 1996")


fit_arima <- Arima(y_in, order = c(0,1,0), seasonal = c(0,1,1))

accuracy(forecast::forecast(fit_arima,h=8), y_out)


checkresiduals(fit_arima)



plot(forecast(test_fit_ACF,h=8))


# Theta Method


fit_theta <- dotm(y_in,h=8,level=c(80,95))

fcs_theta <- accuracy(fit_theta$mean, y_out) 





plot(fit_theta)
accuracy(fcs_nn, y_out)

# Show Down of the three model

fcs_ETS_final <- forecast(ets(y_in, model="AAN", damped=TRUE),h=8)
fcs_ARIMA_final <- forecast::forecast(fit_arima, h=8)
fcs_Theta_final <- fcs_theta

final_model <- c("ETS","ARIMA","Theta")
final_MAPE <- c(accuracy(fcs_ETS_final, y_out)[2,5],
        accuracy(fcs_ARIMA_final, y_out)[2,5],
        fcs_Theta_final[5])

final_performance <- tibble(Model = final_model, MAPE= final_MAPE)

final_performance

plot(y_raw, main = "Original TS")
plot(fcs_ETS_final)
plot(fcs_ARIMA_final)
plot(fit_theta)

