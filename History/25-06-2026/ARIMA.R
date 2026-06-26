# ============================================================================
# ARIMA.R — Modelizacion ARIMA de Inflacion en Bolivia
# Dos metodologias en paralelo:
#   M1 "Directo"     : ARIMA sobre inflacion mensual (var. % del IPC)
#   M2 "Logaritmico" : ARIMA sobre log(IPC), pronostico via exp()
#
# Seleccion por scoring compuesto (RMSE, MAPE, AICc, BIC, HQ, LB, parsimonia)
# ============================================================================

# ---- 0. Configuracion ----

required <- c("readxl","dplyr","lubridate","ggplot2","forecast",
              "tseries","urca","lmtest","scales","openxlsx")
new_pkg  <- required[!(required %in% installed.packages()[,"Package"])]
if (length(new_pkg)) install.packages(new_pkg, repos = "https://cloud.r-project.org")

library(readxl);  library(dplyr);    library(lubridate)
library(ggplot2); library(forecast); library(tseries)
library(urca);    library(lmtest);   library(scales)
library(openxlsx)

tema <- theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        panel.grid.minor = element_blank())
theme_set(tema)

RUTA_BD     <- "Inputs/BD Bruta.xlsx"
RUTA_OUTPUT <- "Outputs/Proyecciones.xlsx"
ANIO_INICIO <- 2023
H_BACKTEST  <- 6
H_FORECAST  <- 7
MAX_P       <- 4
MAX_Q       <- 4

# ---- 1. Carga de datos ----

meses_map <- c("Enero"=1,"Febrero"=2,"Marzo"=3,"Abril"=4,"Mayo"=5,
               "Junio"=6,"Julio"=7,"Agosto"=8,"Septiembre"=9,
               "Octubre"=10,"Noviembre"=11,"Diciembre"=12)

ipc_raw <- read_excel(RUTA_BD, sheet = "IPC")
ipc <- ipc_raw %>%
  mutate(Valor_num = suppressWarnings(as.numeric(Valor))) %>%
  filter(!is.na(Valor_num)) %>%
  mutate(mes_num = meses_map[Mes],
         fecha = as.Date(paste(Año, mes_num, "01", sep = "-")),
         IPC = Valor_num) %>%
  arrange(fecha) %>%
  select(fecha, Año, Mes, mes_num, IPC) %>%
  mutate(log_IPC        = log(IPC),
         inf_mensual    = (IPC / lag(IPC) - 1) * 100,
         inf_interanual = (IPC / lag(IPC, 12) - 1) * 100)

tc_raw <- read_excel(RUTA_BD, sheet = "TC")
names(tc_raw) <- trimws(names(tc_raw))
tc_mensual <- tc_raw %>%
  mutate(fecha_date = as.Date(fecha), anio = year(fecha_date), mes = month(fecha_date)) %>%
  group_by(anio, mes) %>%
  summarise(TC_preferencial = mean(`TC PREFERENCIAL`, na.rm=TRUE),
            TC_digital = mean(`TC-DIGITAL`, na.rm=TRUE),
            TC_oficial = first(`TC-OFICIAL`),
            TC_referencial = mean(`TC-REFERENCIAL`, na.rm=TRUE), .groups="drop") %>%
  mutate(fecha = as.Date(paste(anio, mes, "01", sep = "-")))

ipc_analisis <- ipc %>% filter(year(fecha) >= ANIO_INICIO, !is.na(inf_mensual))

# ---- 2. Series de analisis ----

## M1: Inflacion mensual directa
inf_ts <- ts(ipc_analisis$inf_mensual,
             start = c(year(min(ipc_analisis$fecha)),
                       month(min(ipc_analisis$fecha))),
             frequency = 12)

## M2: log(IPC) — incluye 1 obs mas (no pierde la primera por lag)
ipc_log <- ipc %>% filter(year(fecha) >= ANIO_INICIO, !is.na(IPC))
log_ipc_ts <- ts(ipc_log$log_IPC,
                 start = c(year(min(ipc_log$fecha)),
                           month(min(ipc_log$fecha))),
                 frequency = 12)
dlog_ipc_ts <- diff(log_ipc_ts)

cat("Serie M1 (inf. mensual %): n =", length(inf_ts), "\n")
cat("Serie M2 (log IPC)       : n =", length(log_ipc_ts), "\n")
cat("Serie M2 (Dlog IPC)      : n =", length(dlog_ipc_ts), "\n")

# ---- 3. Analisis preliminar ----

par(mfrow = c(2, 2))
plot(inf_ts, main = "M1: Inflacion mensual (%)", ylab = "Var. %", col = "#2980b9")
abline(h = 0, col = "red", lty = 2)
plot(log_ipc_ts, main = "M2: log(IPC)", ylab = "log(IPC)", col = "#8e44ad")
plot(dlog_ipc_ts, main = expression("M2: "*Delta*"log(IPC)"), ylab = "Dlog", col = "#27ae60")
abline(h = 0, col = "red", lty = 2)
hist(as.numeric(inf_ts), breaks = 15, prob = TRUE,
     main = "M1: Distribucion inf. mensual", col = "#3498db80", border = "white", xlab = "%")
curve(dnorm(x, mean(inf_ts), sd(inf_ts)), add = TRUE, col = "#e74c3c", lwd = 2)
par(mfrow = c(1, 1))

# ---- 4. Pruebas de raiz unitaria ----

cat("\n===== RAIZ UNITARIA =====\n")
cat("\n--- M1: Inflacion mensual ---\n")
cat(sprintf("  ADF  p=%.4f | PP p=%.4f | KPSS p=%.4f\n",
            adf.test(inf_ts)$p.value, pp.test(inf_ts)$p.value,
            kpss.test(inf_ts, null="Level")$p.value))
cat(sprintf("  ndiffs = %d\n", ndiffs(inf_ts, test = "adf")))

cat("\n--- M2: log(IPC) ---\n")
cat(sprintf("  ADF  p=%.4f | PP p=%.4f | KPSS p=%.4f\n",
            adf.test(log_ipc_ts)$p.value, pp.test(log_ipc_ts)$p.value,
            kpss.test(log_ipc_ts, null="Level")$p.value))
cat(sprintf("  ndiffs = %d\n", ndiffs(log_ipc_ts, test = "adf")))

cat("\n--- M2: Dlog(IPC) ---\n")
cat(sprintf("  ADF  p=%.4f | PP p=%.4f | KPSS p=%.4f\n",
            adf.test(dlog_ipc_ts)$p.value, pp.test(dlog_ipc_ts)$p.value,
            kpss.test(dlog_ipc_ts, null="Level")$p.value))

# ---- 5. ACF / PACF ----

par(mfrow = c(2, 2))
Acf(inf_ts, lag.max = 24, main = "ACF — M1: Inf. mensual")
Pacf(inf_ts, lag.max = 24, main = "PACF — M1: Inf. mensual")
Acf(dlog_ipc_ts, lag.max = 24, main = "ACF — M2: Dlog(IPC)")
Pacf(dlog_ipc_ts, lag.max = 24, main = "PACF — M2: Dlog(IPC)")
par(mfrow = c(1, 1))

# ---- 6. Funcion de estimacion + scoring ----

estimar_grilla <- function(serie, max_p = MAX_P, max_q = MAX_Q, etiqueta = "") {
  cat(sprintf("\n--- Estimando grilla %s (n=%d) ---\n", etiqueta, length(serie)))
  res <- data.frame()
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      for (d in 0:1) {
        if (p == 0 && d == 0 && q == 0) next
        fit <- try(Arima(serie, order = c(p, d, q), method = "ML"), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          acc <- accuracy(fit)
          ll  <- as.numeric(logLik(fit))
          k   <- attr(logLik(fit), "df"); n <- fit$nobs
          lb12 <- tryCatch(
            Box.test(residuals(fit), lag = 12, type = "Ljung-Box",
                     fitdf = length(coef(fit)))$p.value,
            error = function(e) NA)
          arc <- fit$model$phi; mac <- fit$model$theta; rok <- TRUE
          if (length(arc) > 0 && any(arc != 0))
            rok <- rok && all(Mod(1 / polyroot(c(1, -arc))) < 1)
          if (length(mac) > 0 && any(mac != 0))
            rok <- rok && all(Mod(1 / polyroot(c(1,  mac))) < 1)
          res <- rbind(res, data.frame(
            Modelo = sprintf("ARIMA(%d,%d,%d)", p, d, q), p=p, d=d, q=q,
            AIC = AIC(fit), BIC = BIC(fit), AICc = fit$aicc,
            HQ = -2*ll + 2*k*log(log(n)),
            RMSE = acc[1,"RMSE"], MAE = acc[1,"MAE"], MAPE = acc[1,"MAPE"],
            LogLik = ll, LB12 = lb12, Raices_ok = rok, Complexity = p + q,
            stringsAsFactors = FALSE))
        }
      }
    }
  }
  res <- res %>%
    mutate(RMSE_std = as.numeric(scale(RMSE)),
           MAPE_std = as.numeric(scale(MAPE)),
           AICc_std = as.numeric(scale(AICc)),
           BIC_std  = as.numeric(scale(BIC)),
           HQ_std   = as.numeric(scale(HQ)),
           Complexity_std = as.numeric(scale(Complexity)),
           LB_pen   = ifelse(is.na(LB12) | LB12 < 0.05, 1, 0),
           R_pen    = ifelse(!Raices_ok, 1, 0),
           Score = RMSE_std + MAPE_std + (AICc_std + BIC_std + HQ_std)/3 +
                   0.3*Complexity_std + 2*LB_pen + 2*R_pen) %>%
    arrange(Score)
  cat(sprintf("  Modelos estimados: %d | Mejor: %s (Score=%.3f)\n",
              nrow(res), res$Modelo[1], res$Score[1]))
  res
}

# ---- 7. Estimacion M1 y M2 ----

cat("\n===== ESTIMACION =====\n")

results_m1 <- estimar_grilla(inf_ts,      etiqueta = "M1: Inf. mensual")
results_m2 <- estimar_grilla(log_ipc_ts,  etiqueta = "M2: log(IPC)")

mejor_m1 <- results_m1[1, ]
mejor_m2 <- results_m2[1, ]

cat("\n--- COMPARACION DE GANADORES ---\n")
cat(sprintf("  M1 (Directo)     : %s  Score=%.3f  AICc=%.2f  RMSE=%.6f  MAPE=%.4f\n",
            mejor_m1$Modelo, mejor_m1$Score, mejor_m1$AICc, mejor_m1$RMSE, mejor_m1$MAPE))
cat(sprintf("  M2 (Logaritmico) : %s  Score=%.3f  AICc=%.2f  RMSE=%.6f  MAPE=%.4f\n",
            mejor_m2$Modelo, mejor_m2$Score, mejor_m2$AICc, mejor_m2$RMSE, mejor_m2$MAPE))

modelo_m1 <- Arima(inf_ts,     order = c(mejor_m1$p, mejor_m1$d, mejor_m1$q), method = "ML")
modelo_m2 <- Arima(log_ipc_ts, order = c(mejor_m2$p, mejor_m2$d, mejor_m2$q), method = "ML")

auto_m1 <- auto.arima(inf_ts, seasonal=TRUE, stepwise=FALSE, approximation=FALSE)
auto_m2 <- auto.arima(log_ipc_ts, seasonal=TRUE, stepwise=FALSE, approximation=FALSE)
cat(sprintf("\n  auto.arima M1: ARIMA(%s)  AIC=%.2f\n",
            paste(arimaorder(auto_m1), collapse=","), AIC(auto_m1)))
cat(sprintf("  auto.arima M2: ARIMA(%s)  AIC=%.2f\n",
            paste(arimaorder(auto_m2), collapse=","), AIC(auto_m2)))

# ---- 8. Diagnostico de ambos ganadores ----

cat("\n===== DIAGNOSTICO =====\n")

for (info in list(list(mod = modelo_m1, nom = paste("M1:", mejor_m1$Modelo)),
                  list(mod = modelo_m2, nom = paste("M2:", mejor_m2$Modelo)))) {
  cat(sprintf("\n--- %s ---\n", info$nom))
  res <- residuals(info$mod)
  coefs <- coef(info$mod); se <- sqrt(diag(vcov(info$mod)))
  z <- coefs / se; pv <- 2*(1-pnorm(abs(z)))
  for (i in seq_along(coefs)) {
    sig <- ifelse(pv[i]<0.01,"***",ifelse(pv[i]<0.05,"**",ifelse(pv[i]<0.1,"*","   ")))
    cat(sprintf("  %-12s est=%9.6f se=%9.6f z=%7.3f p=%.4f %s\n",
                names(coefs)[i], coefs[i], se[i], z[i], pv[i], sig))
  }
  for (h in c(6, 12, 18)) {
    lb <- Box.test(res, lag=h, type="Ljung-Box", fitdf=length(coefs))
    cat(sprintf("  LB h=%2d: p=%.4f -> %s\n", h, lb$p.value,
                ifelse(lb$p.value > 0.05, "OK", "AUTOCORR")))
  }
  jb <- jarque.bera.test(res)
  cat(sprintf("  JB: p=%.4f  SW: p=%.4f\n", jb$p.value, shapiro.test(res)$p.value))
}

# ---- 9. Backtesting ambos metodos ----

cat("\n##########################################################\n")
cat("#                    BACKTESTING                         #\n")
cat("##########################################################\n")

ultimo_mes <- max(ipc_analisis$fecha)
meses_pron <- seq.Date(ultimo_mes %m+% months(1), by = "month", length.out = H_FORECAST)
fechas_test_m1  <- tail(ipc_analisis$fecha, H_BACKTEST)
fechas_train_m1 <- head(ipc_analisis$fecha, length(inf_ts) - H_BACKTEST)

## M1: backtest en niveles de IPC
n1 <- length(inf_ts); nt1 <- n1 - H_BACKTEST
train_m1 <- window(inf_ts, end = time(inf_ts)[nt1])
test_m1  <- window(inf_ts, start = time(inf_ts)[nt1 + 1])

bt_m1 <- Arima(train_m1, order = c(mejor_m1$p, mejor_m1$d, mejor_m1$q), method = "ML")
fc_m1 <- forecast(bt_m1, h = H_BACKTEST, level = c(80, 95))

ultimo_ipc_train <- ipc_analisis$IPC[nt1]
fc_inf_m1 <- as.numeric(fc_m1$mean)
ipc_fc_m1 <- numeric(H_BACKTEST)
ipc_fc_m1[1] <- ultimo_ipc_train * (1 + fc_inf_m1[1] / 100)
for (i in 2:H_BACKTEST) ipc_fc_m1[i] <- ipc_fc_m1[i-1] * (1 + fc_inf_m1[i] / 100)

ipc_obs <- tail(ipc_analisis$IPC, H_BACKTEST)

## M2: backtest en niveles de IPC
fechas_test_m2  <- tail(ipc_log$fecha, H_BACKTEST)
n2 <- length(log_ipc_ts); nt2 <- n2 - H_BACKTEST
train_m2 <- window(log_ipc_ts, end = time(log_ipc_ts)[nt2])
test_m2  <- window(log_ipc_ts, start = time(log_ipc_ts)[nt2 + 1])

bt_m2 <- Arima(train_m2, order = c(mejor_m2$p, mejor_m2$d, mejor_m2$q), method = "ML")
fc_m2 <- forecast(bt_m2, h = H_BACKTEST, level = c(80, 95))
ipc_fc_m2 <- exp(as.numeric(fc_m2$mean))

## Metricas comparadas (ambas en IPC)
err_m1 <- ipc_fc_m1 - ipc_obs;  err_m2 <- ipc_fc_m2 - ipc_obs
rmse_m1 <- sqrt(mean(err_m1^2)); rmse_m2 <- sqrt(mean(err_m2^2))
mape_m1 <- mean(abs(err_m1/ipc_obs))*100; mape_m2 <- mean(abs(err_m2/ipc_obs))*100

lo95_m2 <- exp(as.numeric(fc_m2$lower[,2])); hi95_m2 <- exp(as.numeric(fc_m2$upper[,2]))
dentro_m2 <- ipc_obs >= lo95_m2 & ipc_obs <= hi95_m2

cat("\n--- Backtest M1 (Directo) ---\n")
cat(sprintf("  RMSE IPC = %.4f | MAPE IPC = %.2f%%\n", rmse_m1, mape_m1))

cat("\n--- Backtest M2 (Logaritmico) ---\n")
cat(sprintf("  RMSE IPC = %.4f | MAPE IPC = %.2f%%\n", rmse_m2, mape_m2))
cat(sprintf("  Dentro IC95: %d/%d\n", sum(dentro_m2), H_BACKTEST))

cat("\n--- Comparacion punto a punto (IPC) ---\n")
cat(sprintf("%-12s %10s %10s %10s\n", "Fecha", "IPC obs", "M1 pron", "M2 pron"))
for (i in 1:H_BACKTEST)
  cat(sprintf("%-12s %10.2f %10.2f %10.2f\n",
              format(fechas_test_m1[i], "%b %Y"), ipc_obs[i], ipc_fc_m1[i], ipc_fc_m2[i]))

# ---- 10. Pronostico final ambos metodos ----

cat("\n##########################################################\n")
cat("#                 PRONOSTICO FINAL                       #\n")
cat("##########################################################\n")

## M1: pronostico
fc_final_m1 <- forecast(modelo_m1, h = H_FORECAST, level = c(80, 95))
ultimo_ipc <- tail(ipc_analisis$IPC, 1)
ipc_proy_m1 <- numeric(H_FORECAST)
fc_inf <- as.numeric(fc_final_m1$mean)
ipc_proy_m1[1] <- ultimo_ipc * (1 + fc_inf[1] / 100)
for (i in 2:H_FORECAST) ipc_proy_m1[i] <- ipc_proy_m1[i-1] * (1 + fc_inf[i] / 100)

## M2: pronostico
fc_final_m2 <- forecast(modelo_m2, h = H_FORECAST, level = c(80, 95))
ipc_proy_m2     <- exp(as.numeric(fc_final_m2$mean))
ipc_proy_m2_lo95 <- exp(as.numeric(fc_final_m2$lower[, 2]))
ipc_proy_m2_hi95 <- exp(as.numeric(fc_final_m2$upper[, 2]))
ipc_proy_m2_lo80 <- exp(as.numeric(fc_final_m2$lower[, 1]))
ipc_proy_m2_hi80 <- exp(as.numeric(fc_final_m2$upper[, 1]))

## Inflacion interanual
ipc_hist_all <- ipc$IPC[!is.na(ipc$IPC)]
fechas_hist  <- ipc$fecha[!is.na(ipc$IPC)]

derivar_ia <- function(ipc_proy) {
  ipc_c <- c(ipc_hist_all, ipc_proy); fc <- c(fechas_hist, meses_pron)
  sapply(seq_along(meses_pron), function(i) {
    (ipc_proy[i] / ipc_c[which.min(abs(fc - (meses_pron[i] %m-% months(12))))] - 1) * 100
  })
}

ia_m1 <- derivar_ia(ipc_proy_m1)
ia_m2 <- derivar_ia(ipc_proy_m2)

cat("\n--- Proyeccion comparada (IPC) ---\n")
cat(sprintf("%-15s %10s %10s %12s %12s\n",
            "Fecha", "M1 IPC", "M2 IPC", "M1 Inf.ia%", "M2 Inf.ia%"))
for (i in 1:H_FORECAST)
  cat(sprintf("%-15s %10.2f %10.2f %12.2f %12.2f\n",
              format(meses_pron[i], "%B %Y"),
              ipc_proy_m1[i], ipc_proy_m2[i], ia_m1[i], ia_m2[i]))

# ---- 11. Graficos finales ----

## Fan chart — enlazando historico con proyeccion
ipc_obs_df <- ipc %>% filter(year(fecha) >= ANIO_INICIO, !is.na(IPC)) %>%
  select(fecha, IPC)
ci_fan <- data.frame(fecha = meses_pron,
                     lo80 = ipc_proy_m2_lo80, hi80 = ipc_proy_m2_hi80,
                     lo95 = ipc_proy_m2_lo95, hi95 = ipc_proy_m2_hi95)

pto_enlace <- tail(ipc_obs_df, 1)
m1_linea <- data.frame(fecha = c(pto_enlace$fecha, meses_pron),
                       IPC   = c(pto_enlace$IPC,   ipc_proy_m1))
m2_linea <- data.frame(fecha = c(pto_enlace$fecha, meses_pron),
                       IPC   = c(pto_enlace$IPC,   ipc_proy_m2))

print(ggplot() +
  geom_ribbon(data = ci_fan, aes(fecha, ymin=lo95, ymax=hi95), fill="#2980b9", alpha=0.10) +
  geom_ribbon(data = ci_fan, aes(fecha, ymin=lo80, ymax=hi80), fill="#2980b9", alpha=0.22) +
  geom_line(data = ipc_obs_df, aes(fecha, IPC), color="#2c3e50", linewidth=0.9) +
  geom_line(data = m2_linea, aes(fecha, IPC), color="blue", linewidth=0.9) +
  geom_line(data = m1_linea, aes(fecha, IPC), color="#e74c3c", linewidth=0.9, linetype="dashed") +
  labs(title = "IPC proyectado — M1 (rojo) vs M2 (azul)",
       subtitle = sprintf("M1: %s | M2: %s", mejor_m1$Modelo, mejor_m2$Modelo),
       x = NULL, y = "IPC") +
  scale_x_date(date_breaks="4 months", date_labels="%b %Y") +
  theme(axis.text.x = element_text(angle=45, hjust=1)))

## Inflacion interanual — enlazando historico con proyeccion
ia_h <- ipc %>% filter(!is.na(inf_interanual), year(fecha) >= ANIO_INICIO) %>%
  select(fecha, yoy = inf_interanual)

ult_ia <- tail(ia_h, 1)
ia_m1_linea <- data.frame(fecha = c(ult_ia$fecha, meses_pron),
                          yoy   = c(ult_ia$yoy,   ia_m1), tipo = "M1 Directo")
ia_m2_linea <- data.frame(fecha = c(ult_ia$fecha, meses_pron),
                          yoy   = c(ult_ia$yoy,   ia_m2), tipo = "M2 Logaritmico")

ia_all <- bind_rows(ia_h %>% mutate(tipo = "Observado"), ia_m1_linea, ia_m2_linea)

print(ggplot(ia_all, aes(fecha, yoy, color=tipo, linetype=tipo)) +
  geom_line(linewidth=0.9) + geom_point(size=1.2) +
  scale_color_manual(values=c(Observado="black", `M1 Directo`="#e74c3c", `M2 Logaritmico`="blue")) +
  scale_linetype_manual(values=c(Observado="solid", `M1 Directo`="dashed", `M2 Logaritmico`="solid")) +
  labs(title="Inflacion interanual — M1 vs M2", x=NULL, y="% interanual", color=NULL, linetype=NULL) +
  scale_x_date(date_breaks="4 months", date_labels="%b %Y") +
  theme(axis.text.x=element_text(angle=45,hjust=1)))

# ---- 12. Exportar ----

cat("\n===== EXPORTACION =====\n")

mn <- c("Enero","Febrero","Marzo","Abril","Mayo","Junio",
        "Julio","Agosto","Septiembre","Octubre","Noviembre","Diciembre")

hoja_proy <- data.frame(
  Anio = year(meses_pron), Mes = mn[month(meses_pron)],
  M1_IPC = round(ipc_proy_m1, 4), M2_IPC = round(ipc_proy_m2, 4),
  M1_Inf_ia = round(ia_m1, 2), M2_Inf_ia = round(ia_m2, 2),
  M2_IC95_inf = round(ipc_proy_m2_lo95, 2), M2_IC95_sup = round(ipc_proy_m2_hi95, 2))

hoja_bt <- data.frame(
  Fecha = format(fechas_test_m1, "%b %Y"), IPC_obs = round(ipc_obs, 2),
  M1_pron = round(ipc_fc_m1, 2), M2_pron = round(ipc_fc_m2, 2),
  M1_err = round(err_m1, 2), M2_err = round(err_m2, 2))

hoja_scoring_m1 <- results_m1 %>%
  select(Modelo,p,d,q,AICc,BIC,HQ,RMSE,MAPE,LB12,Score) %>%
  mutate(across(where(is.numeric), ~round(.,4)), Metodo = "M1 Directo")
hoja_scoring_m2 <- results_m2 %>%
  select(Modelo,p,d,q,AICc,BIC,HQ,RMSE,MAPE,LB12,Score) %>%
  mutate(across(where(is.numeric), ~round(.,4)), Metodo = "M2 Logaritmico")

hoja_comp <- data.frame(
  Criterio = c("Modelo ganador","Score","AICc","BIC","HQ","RMSE (train)","MAPE (train)",
               "RMSE backtest (IPC)","MAPE backtest (IPC)"),
  M1_Directo = c(mejor_m1$Modelo, round(mejor_m1$Score,3), round(mejor_m1$AICc,2),
                 round(mejor_m1$BIC,2), round(mejor_m1$HQ,2),
                 round(mejor_m1$RMSE,6), round(mejor_m1$MAPE,4),
                 round(rmse_m1,4), paste0(round(mape_m1,2),"%")),
  M2_Logaritmico = c(mejor_m2$Modelo, round(mejor_m2$Score,3), round(mejor_m2$AICc,2),
                     round(mejor_m2$BIC,2), round(mejor_m2$HQ,2),
                     round(mejor_m2$RMSE,6), round(mejor_m2$MAPE,4),
                     round(rmse_m2,4), paste0(round(mape_m2,2),"%")))

## Serie completa: historico + proyectado
hoja_serie <- bind_rows(
  ipc %>%
    filter(year(fecha) >= ANIO_INICIO, !is.na(IPC)) %>%
    transmute(Anio = Año, Mes = Mes, IPC = round(IPC, 4),
              M1_IPC = round(IPC, 4), M2_IPC = round(IPC, 4),
              Inf_interanual = round(inf_interanual, 2),
              Tipo = "HISTORICO"),
  data.frame(
    Anio = year(meses_pron), Mes = mn[month(meses_pron)],
    IPC = NA,
    M1_IPC = round(ipc_proy_m1, 4), M2_IPC = round(ipc_proy_m2, 4),
    Inf_interanual = NA,
    Tipo = "PROYECCION")
) %>%
  mutate(
    M1_Inf_ia = round(ia_m1[match(paste(Anio, Mes), paste(year(meses_pron), mn[month(meses_pron)]))], 2),
    M2_Inf_ia = round(ia_m2[match(paste(Anio, Mes), paste(year(meses_pron), mn[month(meses_pron)]))], 2)
  ) %>%
  mutate(
    M1_Inf_ia = ifelse(Tipo == "HISTORICO", Inf_interanual, M1_Inf_ia),
    M2_Inf_ia = ifelse(Tipo == "HISTORICO", Inf_interanual, M2_Inf_ia)
  ) %>%
  select(Anio, Mes, IPC, M1_IPC, M2_IPC, M1_Inf_ia, M2_Inf_ia, Tipo)

wb <- createWorkbook()
addWorksheet(wb,"Serie_Completa"); writeData(wb,"Serie_Completa", hoja_serie)
addWorksheet(wb,"Proyecciones");   writeData(wb,"Proyecciones",   hoja_proy)
addWorksheet(wb,"Backtesting");    writeData(wb,"Backtesting",    hoja_bt)
addWorksheet(wb,"Scoring_M1");     writeData(wb,"Scoring_M1",     hoja_scoring_m1)
addWorksheet(wb,"Scoring_M2");     writeData(wb,"Scoring_M2",     hoja_scoring_m2)
addWorksheet(wb,"Comparacion");    writeData(wb,"Comparacion",    hoja_comp)
saveWorkbook(wb, RUTA_OUTPUT, overwrite = TRUE)
cat(sprintf("Exportado: %s\n", RUTA_OUTPUT))
cat("  Hoja 'Serie_Completa': historico + proyectado (M1 y M2)\n")

# ---- 13. Resumen ----

cat("\n==========================================================\n")
cat("                    RESUMEN FINAL\n")
cat("==========================================================\n")
cat(sprintf("  %-22s %-20s %-20s\n", "", "M1 DIRECTO", "M2 LOGARITMICO"))
cat(sprintf("  %-22s %-20s %-20s\n", "Modelo",     mejor_m1$Modelo, mejor_m2$Modelo))
cat(sprintf("  %-22s %-20.3f %-20.3f\n", "Score",  mejor_m1$Score,  mejor_m2$Score))
cat(sprintf("  %-22s %-20.2f %-20.2f\n", "AICc",   mejor_m1$AICc,   mejor_m2$AICc))
cat(sprintf("  %-22s %-20.4f %-20.4f\n", "BT RMSE (IPC)", rmse_m1, rmse_m2))
cat(sprintf("  %-22s %-20s %-20s\n", "BT MAPE (IPC)",
            paste0(round(mape_m1,2),"%"), paste0(round(mape_m2,2),"%")))
cat(sprintf("  %-22s %-20.2f %-20.2f\n", "IPC dic 2026",
            tail(ipc_proy_m1,1), tail(ipc_proy_m2,1)))
cat(sprintf("  %-22s %-20.2f %-20.2f\n", "Inf.ia dic 2026",
            tail(ia_m1,1), tail(ia_m2,1)))
cat("==========================================================\n")
