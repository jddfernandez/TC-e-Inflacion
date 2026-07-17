#install.packages("forecast")
#install.packages("tseries")

library(readxl)
library(tidyverse)
library(forecast)
library(tseries)
library(lubridate)
library(ggplot2)
library(dplyr)
library(zoo)
library(scales)

setwd("/Users/sergiosalazar/Documents/Ministerio de Economía/Sistema Financiero")
base <- read_excel("depos_cartera.xlsx")
base$Fecha <- as.Date(base$Fecha)
base <- base %>%
  filter(`Créditos` > 0)

base <- base %>%
  mutate(
    l_depositos = log(Depositos),
    l_creditos  = log(`Créditos`)
  )
ts_dep <- ts(base$l_depositos,
             start = c(2010,12),
             frequency = 12)

ts_cred <- ts(base$l_creditos,
              start = c(2010,12),
              frequency = 12)
autoplot(ts_dep)
autoplot(ts_cred)

# Raiz Unitaria 
adf.test(ts_dep)
adf.test(ts_cred)

adf.test(diff(ts_dep))
adf.test(diff(ts_cred))

autoplot(diff(ts_dep))
autoplot(diff(ts_cred))

ggAcf(diff(ts_dep))
ggPacf(diff(ts_dep))


fit1 <- Arima(ts_dep,
              order=c(0,1,1)) # ARIMA(0,1,1)
autoplot(fit1)
accuracy(fit1)
checkresiduals(fit1)

fit2 <- Arima(ts_dep,
              order=c(1,1,1)) # ARIMA(1,1,1)
autoplot(fit2)
accuracy(fit2)
checkresiduals(fit2)

fit3 <- Arima(ts_dep,
              order=c(0,1,1),
              seasonal=c(0,1,1)) # ARIMA(0,1,1)(0,1,1)
autoplot(fit3)
accuracy(fit3)
checkresiduals(fit3)

fit4 <- Arima(ts_dep,
              order=c(1,1,1),
              seasonal=c(1,1,1))
autoplot(fit4)
accuracy(fit4)
checkresiduals(fit4)

fit5 <- Arima(ts_cred,
              order=c(0,1,1)) # ARIMA(0,1,1)
autoplot(fit5)
accuracy(fit5)
checkresiduals(fit5)

fit6 <- Arima(ts_cred,
              order=c(1,1,1)) # ARIMA(1,1,1)
autoplot(fit6)
accuracy(fit6)
checkresiduals(fit6)

fit7 <- Arima(ts_cred,
              order=c(0,1,1),
              seasonal=c(0,1,1)) # ARIMA(0,1,1)(0,1,1)
autoplot(fit7)
accuracy(fit7)
checkresiduals(fit7)

fit8 <- Arima(ts_cred,
              order=c(1,1,1),
              seasonal=c(1,1,1))
autoplot(fit8)
accuracy(fit8)
checkresiduals(fit8)


serie<-ts_dep

results <- data.frame()

for(p in 0:4){
  for(q in 0:4){
    for(P in 0:0){
      for(Q in 0:2){
        nombre <- paste0(
          "sarima",
          p,q,P,Q
        )
        fit <- try(
          Arima(
            serie,
            order = c(p,1,q),
            seasonal = c(P,1,Q)
          ),
          silent = TRUE
        )
        if(!inherits(fit,"try-error")){
          assign(nombre, fit)
          fc <- forecast(fit)
          acc <- accuracy(fc)
          rmse <- acc[1,"RMSE"]
          mape <- acc[1,"MAPE"]
          lb <- Box.test(
            residuals(fit),
            lag = 12,
            type = "Ljung-Box"
          )
          lb_pvalue <- lb$p.value
          complexity <- p + q + P + Q
          results <- rbind(
            results,
            data.frame(
              modelo = nombre,
              p = p,
              q = q,
              P = P,
              Q = Q,
              RMSE = rmse,
              MAPE = mape,
              LjungBox = lb_pvalue,
              Complexity = complexity
            )
          )
        }
      }
    }
  }
}
results$RMSE_std <- as.numeric(scale(results$RMSE))
results$MAPE_std <- as.numeric(scale(results$MAPE))
results$Complexity_std <- as.numeric(scale(results$Complexity))
results$LB_penalty <- ifelse(
  results$LjungBox < 0.05, 1, 0)
results$Score <-
  results$RMSE_std +
  results$MAPE_std +
  0.3*results$Complexity_std +
  2*results$LB_penalty
results <- results %>%
  arrange(Score)
best_row <- results[1,]

best_model <- get(best_row$modelo)
best_fc <- forecast(
  best_model,
  h = 12,
  level = c(25,50)
)

fc_mean  <- exp(best_fc$mean)
fc_lower <- exp(best_fc$lower)
fc_upper <- exp(best_fc$upper)

hist_levels <- exp(serie)

checkresiduals(best_model)
accuracy(best_model)

autoplot(best_model)
autoplot(best_fc, 
         series = "Forecast") +
  autolayer(
    serie,
    series = c("Observed")
  ) +
  ggtitle(
    paste0(
      "Fan Chart Forecast - ",
      best_row$modelo
    )
  ) +
  xlab("") +
  ylab("Log Deposits") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

yoy_hist <- diff(serie, lag=12)*100
full_series <- ts(
  c(serie,best_fc$mean),
  start = start(serie),
  frequency = 12
)
yoy_full <- diff(
  full_series,
  lag = 12)*100
forecast_start <- length(yoy_hist)

yoy_fc <- window(
  yoy_full,
  start = time(yoy_full)[forecast_start]
)


autoplot(yoy_hist,
         series = "Observed") +
  autolayer(
    yoy_fc,
    series = "Forecast"
  ) +
  ggtitle(
    "Deposits Forecast - YoY Growth"
  ) +
  xlab("") +
  ylab("% YoY") +
  theme_minimal()

tail(yoy_fc,22)
h <- length(best_fc$mean)
hist_last12 <- tail(serie, 12)
yoy_mean <- numeric(h)
yoy_lo95 <- numeric(h)
yoy_hi95 <- numeric(h)
yoy_lo80 <- numeric(h)
yoy_hi80 <- numeric(h)

for(i in 1:h){
  if(i <= 12){
    yoy_mean[i] <-
      (best_fc$mean[i] - hist_last12[i]) * 100
    yoy_lo95[i] <-
      (best_fc$lower[i,"50%"] - hist_last12[i]) * 100
    yoy_hi95[i] <-
      (best_fc$upper[i,"50%"] - hist_last12[i]) * 100
    yoy_lo80[i] <-
      (best_fc$lower[i,"25%"] - hist_last12[i]) * 100
    yoy_hi80[i] <-
      (best_fc$upper[i,"25%"] - hist_last12[i]) * 100
  } else {
    yoy_mean[i] <-
      (best_fc$mean[i] - best_fc$mean[i-12]) * 100
    yoy_lo95[i] <-
      (best_fc$lower[i,"50%"] -
         best_fc$lower[i-12,"50%"]) * 100
    yoy_hi95[i] <-
      (best_fc$upper[i,"50%"] -
         best_fc$upper[i-12,"50%"]) * 100
    yoy_lo80[i] <-
      (best_fc$lower[i,"25%"] -
         best_fc$lower[i-12,"25%"]) * 100
    yoy_hi80[i] <-
      (best_fc$upper[i,"25%"] -
         best_fc$upper[i-12,"25%"]) * 100
  }
}

fc_df <- data.frame(
  date = as.yearmon(time(best_fc$mean)),
  mean = as.numeric(yoy_mean),
  lo95 = as.numeric(yoy_lo95),
  hi95 = as.numeric(yoy_hi95),
  lo80 = as.numeric(yoy_lo80),
  hi80 = as.numeric(yoy_hi80)
)
hist_df <- data.frame(
  date = as.yearmon(time(yoy_hist)),
  value = as.numeric(yoy_hist)
)
forecast_color <- "#2C7FB8"
forecast_color_line <- "blue"
start_plot <- as.yearmon("2020-01")
first_fc_obs <- data.frame(
  date = head(fc_df$date, 1),
  value = head(fc_df$mean, 1)
)
hist_line_df <- rbind(
  hist_df,
  first_fc_obs
)
label_point <- subset(
  fc_df,
  format(as.Date(date), "%Y-%m") == "2026-12"
)
ggplot() +
  geom_hline(
    yintercept = 0,
    color = "grey60",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = subset(fc_df, date >= start_plot),
    aes(x = date, ymin = lo95, ymax = hi95),
    fill = forecast_color,
    alpha = 0.10
  ) +
  geom_ribbon(
    data = subset(fc_df, date >= start_plot),
    aes(x = date, ymin = lo80, ymax = hi80),
    fill = forecast_color,
    alpha = 0.22
  ) +
  geom_line(
    data = subset(hist_line_df, date >= start_plot),
    aes(x = date, y = value),
    linewidth = 1.2,
    color = "black"
  ) +
  geom_line(
    data = subset(fc_df, date >= start_plot),
    aes(x = date, y = mean),
    linewidth = 1.2,
    color = forecast_color_line,
    linetype = "solid"
  ) +
  geom_point(
    data = label_point,
    aes(x = date, y = mean),
    color = forecast_color_line,
    size = 3
  ) +
  geom_text(
    data = label_point,
    aes(
      x = date,
      y = mean,
      label = round(mean,1)
    ),
    vjust = 1.5,
    hjust = 1.02,
    color = forecast_color_line,
    size = 6
  ) +
  labs(
    x = "",
    y = "% YoY"
  ) +
  scale_x_yearmon(
    format = "%Y",
    n = 10
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 12
    )
  )


#################################
#################################
#################################
serie <- ts_cred

results <- data.frame()

for(p in 0:6){
  for(q in 0:6){
    nombre <- paste0("arma",p,q
    )
    fit <- try(
      Arima(
        serie,
        order = c(p,1,q),
        include.drift = TRUE
      ),
      silent = TRUE
    )
    if(!inherits(fit,"try-error")){
      assign(nombre, fit)
      fc <- forecast(fit)
      acc <- accuracy(fc)
      rmse <- acc[1,"RMSE"]
      mape <- acc[1,"MAPE"]
      lb <- Box.test(
        residuals(fit),
        lag = 12,
        type = "Ljung-Box"
      )
      lb_pvalue <- lb$p.value
      complexity <- p + q
      results <- rbind(
        results,
        data.frame(
          modelo = nombre,
          p = p,
          q = q,
          RMSE = rmse,
          MAPE = mape,
          LjungBox = lb_pvalue,
          Complexity = complexity
        )
      )
    }
  }
}

results <- results[
  order(
    -results$LjungBox,
    results$MAPE,
    results$RMSE,
    results$Complexity
  ),
]
results$RMSE_std <- as.numeric(scale(results$RMSE))
results$MAPE_std <- as.numeric(scale(results$MAPE))
results$Complexity_std <- as.numeric(scale(results$Complexity))
results$LB_penalty <- ifelse(
  results$LjungBox < 0.05, 1, 0)
results$Score <-
  results$RMSE_std +
  results$MAPE_std +
  0.3*results$Complexity_std +
  2*results$LB_penalty
results <- results %>%
  arrange(Score)
best_row <- results[1,]

best_model <- get(best_row$modelo)
best_fc <- forecast(
  best_model,
  h = 12,
  level = c(25,50)
)

fc_mean  <- exp(best_fc$mean)
fc_lower <- exp(best_fc$lower)
fc_upper <- exp(best_fc$upper)

hist_levels <- exp(serie)

checkresiduals(best_model)
accuracy(best_model)
autoplot(best_model)

autoplot(best_fc, 
         series = "Forecast") +
  autolayer(
    serie,
    series = c("Observed")
  ) +
  ggtitle(
    paste0(
      "Fan Chart Forecast - ",
      best_row$modelo
    )
  ) +
  xlab("") +
  ylab("Log Deposits") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

yoy_hist <- diff(serie, lag=12)*100
full_series <- ts(
  c(serie,best_fc$mean),
  start = start(serie),
  frequency = 12
)
yoy_full <- diff(
  full_series,
  lag = 12)*100
forecast_start <- length(yoy_hist)

yoy_fc <- window(
  yoy_full,
  start = time(yoy_full)[forecast_start]
)

autoplot(yoy_hist,
         series = "Observed") +
  autolayer(
    yoy_fc,
    series = "Forecast"
  ) +
  ggtitle(
    "Loans Forecast - YoY Growth"
  ) +
  xlab("") +
  ylab("% YoY") +
  theme_minimal()

tail(yoy_fc, 22)


h <- length(best_fc$mean)
hist_last12 <- tail(serie, 12)
yoy_mean <- numeric(h)
yoy_lo95 <- numeric(h)
yoy_hi95 <- numeric(h)
yoy_lo80 <- numeric(h)
yoy_hi80 <- numeric(h)

for(i in 1:h){
  if(i <= 12){
    yoy_mean[i] <-
      (best_fc$mean[i] - hist_last12[i]) * 100
    yoy_lo95[i] <-
      (best_fc$lower[i,"50%"] - hist_last12[i]) * 100
    yoy_hi95[i] <-
      (best_fc$upper[i,"50%"] - hist_last12[i]) * 100
    yoy_lo80[i] <-
      (best_fc$lower[i,"25%"] - hist_last12[i]) * 100
    yoy_hi80[i] <-
      (best_fc$upper[i,"25%"] - hist_last12[i]) * 100
  } else {
    yoy_mean[i] <-
      (best_fc$mean[i] - best_fc$mean[i-12]) * 100
    yoy_lo95[i] <-
      (best_fc$lower[i,"50%"] -
         best_fc$lower[i-12,"50%"]) * 100
    yoy_hi95[i] <-
      (best_fc$upper[i,"50%"] -
         best_fc$upper[i-12,"50%"]) * 100
    yoy_lo80[i] <-
      (best_fc$lower[i,"25%"] -
         best_fc$lower[i-12,"25%"]) * 100
    yoy_hi80[i] <-
      (best_fc$upper[i,"25%"] -
         best_fc$upper[i-12,"25%"]) * 100
  }
}

fc_df <- data.frame(
  date = as.yearmon(time(best_fc$mean)),
  mean = as.numeric(yoy_mean),
  lo95 = as.numeric(yoy_lo95),
  hi95 = as.numeric(yoy_hi95),
  lo80 = as.numeric(yoy_lo80),
  hi80 = as.numeric(yoy_hi80)
)
hist_df <- data.frame(
  date = as.yearmon(time(yoy_hist)),
  value = as.numeric(yoy_hist)
)
forecast_color <- "#2C7FB8"
forecast_color_line <- "blue"
start_plot <- as.yearmon("2020-01")
first_fc_obs <- data.frame(
  date = head(fc_df$date, 1),
  value = head(fc_df$mean, 1)
)
hist_line_df <- rbind(
  hist_df,
  first_fc_obs
)
label_point <- subset(
  fc_df,
  format(as.Date(date), "%Y-%m") == "2026-12"
)
ggplot() +
  geom_hline(
    yintercept = 0,
    color = "grey60",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = subset(fc_df, date >= start_plot),
    aes(x = date, ymin = lo95, ymax = hi95),
    fill = forecast_color,
    alpha = 0.10
  ) +
  geom_ribbon(
    data = subset(fc_df, date >= start_plot),
    aes(x = date, ymin = lo80, ymax = hi80),
    fill = forecast_color,
    alpha = 0.22
  ) +
  geom_line(
    data = subset(hist_line_df, date >= start_plot),
    aes(x = date, y = value),
    linewidth = 1.2,
    color = "black"
  ) +
  geom_line(
    data = subset(fc_df, date >= start_plot),
    aes(x = date, y = mean),
    linewidth = 1.2,
    color = forecast_color_line,
    linetype = "solid"
  ) +
  geom_point(
    data = label_point,
    aes(x = date, y = mean),
    color = forecast_color_line,
    size = 2
  ) +
  geom_text(
    data = label_point,
    aes(
      x = date,
      y = mean,
      label = round(mean,1)
    ),
    vjust = -1.5,
    color = forecast_color_line,
    size = 6
  ) +
  labs(
    x = "",
    y = "% YoY"
  ) +
  scale_x_yearmon(
    format = "%Y",
    n = 10
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 12
    )
  )


