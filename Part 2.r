library(Mcomp)
library(ggplot2)
library(fpp)
library(forecast)
library(base)
library(MAPA)
library(forecastHybrid)
library(Metrics)
library(stats)
library(dplyr)
library(tibble)
library(gtools)



origin <- 701:1400

# Abstract ID that ends with 3
array <- NA

for(i in origin){
    if(i%%10 == 3){
        array[which(i==origin)] <- i
    }
}

series_index <- data.frame(array)[!is.na(data.frame(array)$array),]


# Examine 70 series

M3[series_index]

# Create empty array including selected model and three error measures

final_array <- array(NA, c(length(series_index),4))

final_array


# Across all 70 time series
for (i in 1:length(series_index)){
    y_in <- M3[[series_index[i]]]$x
    y_out <- M3[[series_index[i]]]$xx
    
    fit_ets <- ets(y_in)
    fit_arima <- auto.arima(y_in)
    fit_mapa <- mapa(y_in, fh = 8, conf.lvl=c(0.8, 0.95))

    ets_in_MAE <- forecast::accuracy(fit_ets)[3]
    arima_in_MAE <- forecast::accuracy(fit_arima)[3]
    mapa_in_MAE <- fit_mapa$MAE
    
    temp_MAE <- c(ets_in_MAE, arima_in_MAE, mapa_in_MAE)

    if(temp_MAE[which.min(temp_MAE)]==ets_in_MAE){
        
        final_array[i, 1] <- "ETS"
        final_array[i, 2] <- forecast::accuracy(forecast::forecast(fit_ets, h = 8), y_out)[2,4]
        final_array[i, 3] <- forecast::accuracy(forecast::forecast(fit_ets, h = 8), y_out)[2,5]
        final_array[i, 4] <- forecast::accuracy(forecast::forecast(fit_ets, h = 8), y_out)[2,6]

    } else if(temp_MAE[which.min(temp_MAE)]==arima_in_MAE){
        
        final_array[i, 1] <- "ARIMA"
        final_array[i, 2] <- forecast::accuracy(forecast::forecast(fit_arima, h = 8), y_out)[2,4]
        final_array[i, 3] <- forecast::accuracy(forecast::forecast(fit_arima, h = 8), y_out)[2,5]
        final_array[i, 4] <- forecast::accuracy(forecast::forecast(fit_arima, h = 8), y_out)[2,6]

    } else if(temp_MAE[which.min(temp_MAE)]==mapa_in_MAE){
        
        final_array[i, 1] <- "MAPA"
        final_array[i, 2] <- forecast::accuracy(fit_mapa$outfor, y_out)[4]
        final_array[i, 3] <- forecast::accuracy(fit_mapa$outfor, y_out)[5]
        final_array[i, 4] <- mase(y_out, fit_mapa$outfor)

    }
    
}

# Checkpoint
final_array

# Create tibble form, rename column, and convert error measures 
# into numeric type, then calculate average error measures

temp_tibble <- final_array %>% as_tibble()
names(temp_tibble) <- c("Model", "MPE", "MAPE", "MASE")

temp_tibble$MPE <- as.numeric(temp_tibble$MPE)
temp_tibble$MAPE <- as.numeric(temp_tibble$MAPE)
temp_tibble$MASE <- as.numeric(temp_tibble$MASE)
temp_tibble

temp_group <- group_by(temp_tibble, Model)


error_table <- temp_group %>% summarise(
    MPE = mean(MPE),
    MAPE = mean(MAPE),
    MASE = mean(MASE)
)

# Keep count of the number of model chosen

cnt <- count(temp_tibble, Model)
error_table <- error_table %>%add_column(d=cnt['n'], .after ="Model")

names(error_table) <- c("Model", "Count", "MPE", "MAPE", "MASE")


error_table



# Combination Strategy and benchmark performance comparison

hybrid_errors <- array(NA, c(length(series_index),3))
naive_errors <- array(NA, c(length(series_index),3))
HES_errors <- array(NA, c(length(series_index),3))
Theta_errors <- array(NA, c(length(series_index),3))

for(i in 1:length(series_index)){

    y_in <- M3[[series_index[i]]]$x
    y_out <- M3[[series_index[i]]]$xx

    # Hybrid + Naive + Holts ES(Linear Trend)
    hybrid_fit <- hybridModel(y_in, models="ae", weights = "insample.errors")
    naive_fit <- naive(y_in, h=8, level=c(80,85,90,95,99))
    HES_fit <- ets(y_in, model="MAM", damped=TRUE)
    Theta_fit <- thief(y_in, h = 8, comb="struc", usemodel="theta")

    # Generate Forecast
    fcs_hybrid <- forecast::forecast(hybrid_fit, h=8)
    fcs_naive <- forecast::forecast(naive_fit, h=8)
    fcs_HES <- forecast::forecast(HES_fit, h=8)
    fcs_Theta <- forecast::forecast(Theta_fit, h=8)
     
    # Calculate out-sample model accuracy (MPE, MAPE, MASE)
    hybrid_errors[i, 1] <- forecast::accuracy(fcs_hybrid, y_out)[2,4]
    hybrid_errors[i, 2] <- forecast::accuracy(fcs_hybrid, y_out)[2,5] 
    hybrid_errors[i, 3] <- forecast::accuracy(fcs_hybrid, y_out)[2,6] 

    naive_errors[i, 1] <- forecast::accuracy(fcs_naive, y_out)[2,4]
    naive_errors[i, 2] <- forecast::accuracy(fcs_naive, y_out)[2,5]
    naive_errors[i, 3] <- forecast::accuracy(fcs_naive, y_out)[2,6]

    HES_errors[i, 1] <- forecast::accuracy(fcs_HES, y_out)[2,4]
    HES_errors[i, 2] <- forecast::accuracy(fcs_HES, y_out)[2,5]
    HES_errors[i, 3] <- forecast::accuracy(fcs_HES, y_out)[2,6]

    Theta_errors[i, 1] <- forecast::accuracy(fcs_Theta, y_out)[2,4]
    Theta_errors[i, 2] <- forecast::accuracy(fcs_Theta, y_out)[2,5]
    Theta_errors[i, 3] <- forecast::accuracy(fcs_Theta, y_out)[2,6]
    
}

# Change to tibble and rename column name
hybrid_errors <- hybrid_errors %>% as_tibble()
naive_errors <- naive_errors %>% as_tibble()
HES_errors <- HES_errors %>% as_tibble()
Theta_errors <- Theta_errors %>% as_tibble()

names(hybrid_errors) <- c("MPE", "MAPE", "MASE")
names(naive_errors) <- c("MPE", "MAPE", "MASE")
names(HES_errors) <- c("MPE", "MAPE", "MASE")
names(Theta_errors) <- c("MPE", "MAPE", "MASE")

# Calculate average error

HB_error_avg <- colMeans(hybrid_errors)
naive_error_avg <- colMeans(naive_errors)
HES_error_avg <- colMeans(HES_errors)
Theta_error_avg <- colMeans(Theta_errors)

Benchmark_Performance <- bind_rows(HB_error_avg, naive_error_avg, HES_error_avg, Theta_error_avg)

method <- c("Combination Strategy", "Naive", "Damped ES", "Theta")

Benchmark_Performance$Model <- method
Benchmark_Performance <- Benchmark_Performance[, c(4,1,2,3)]
Benchmark_Performance



# Decomposition

demographics <- subset(M3[series_index], 'DEMOGRAPHIC')
finance <- subset(M3[series_index], 'FINANCE')
industry <- subset(M3[series_index], 'INDUSTRY')
macro <- subset(M3[series_index], 'MACRO')
micro <- subset(M3[series_index], "MICRO")














ets(test)
forecast::forecast(ets(test))


