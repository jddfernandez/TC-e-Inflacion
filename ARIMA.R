# ============================================================================
# ARIMA.R — Modelizacion ARIMA de Inflacion en Bolivia
# Metodologia de Box-Jenkins (script interactivo)
#
# IMPORTANTE: Este script NUNCA modifica BD Bruta.xlsx.
#   Los datos historicos se leen como solo-lectura.
#   Las proyecciones se escriben en Outputs/Proyecciones.xlsx
#
# Ejecutar seccion por seccion en RStudio (Ctrl+Enter) para validar
# cada etapa.
# ============================================================================

# ---- 0. Configuracion ----

required <- c("readxl","dplyr","lubridate","ggplot2","forecast",
              "tseries","urca","lmtest","scales","openxlsx")
new_pkg  <- required[!(required %in% installed.packages()[,"Package"])]
if (length(new_pkg)) install.packages(new_pkg, repos = "https://cloud.r-project.org")

library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(forecast)
library(tseries)
library(urca)
library(lmtest)
library(scales)
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
H_FORECAST  <- 7   # jun-dic 2026

# ---- 1. Carga de datos (solo lectura) ----

meses_map <- c("Enero"=1,"Febrero"=2,"Marzo"=3,"Abril"=4,"Mayo"=5,
               "Junio"=6,"Julio"=7,"Agosto"=8,"Septiembre"=9,
               "Octubre"=10,"Noviembre"=11,"Diciembre"=12)

ipc_raw <- read_excel(RUTA_BD, sheet = "IPC")

ipc <- ipc_raw %>%
  mutate(Valor_num = suppressWarnings(as.numeric(Valor))) %>%
  filter(!is.na(Valor_num)) %>%
  mutate(
    mes_num = meses_map[Mes],
    fecha   = as.Date(paste(Año, mes_num, "01", sep = "-")),
    IPC     = Valor_num
  ) %>%
  arrange(fecha) %>%
  select(fecha, Año, Mes, mes_num, IPC)

cat("IPC historico cargado:", nrow(ipc), "obs —",
    format(min(ipc$fecha), "%b %Y"), "a",
    format(max(ipc$fecha), "%b %Y"), "\n")
cat("(BD Bruta.xlsx se usa en modo SOLO LECTURA)\n")

tc_raw <- read_excel(RUTA_BD, sheet = "TC")
names(tc_raw) <- trimws(names(tc_raw))

tc_mensual <- tc_raw %>%
  mutate(fecha_date = as.Date(fecha),
         anio = year(fecha_date),
         mes  = month(fecha_date)) %>%
  group_by(anio, mes) %>%
  summarise(
    TC_preferencial = mean(`TC PREFERENCIAL`, na.rm = TRUE),
    TC_digital      = mean(`TC-DIGITAL`,      na.rm = TRUE),
    TC_oficial      = first(`TC-OFICIAL`),
    TC_referencial  = mean(`TC-REFERENCIAL`,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(fecha = as.Date(paste(anio, mes, "01", sep = "-")))

cat("TC mensualizado:", nrow(tc_mensual), "obs\n")

# ---- 2. Transformacion de datos ----

ipc <- ipc %>%
  mutate(
    inf_mensual    = (IPC / lag(IPC)     - 1) * 100,
    inf_interanual = (IPC / lag(IPC, 12) - 1) * 100
  )

bd_mensual <- ipc %>%
  filter(year(fecha) >= ANIO_INICIO) %>%
  left_join(tc_mensual %>% select(fecha, starts_with("TC_")), by = "fecha")

cat("\nBase mensual consolidada:\n")
print(bd_mensual, n = Inf)

ipc_analisis <- ipc %>% filter(year(fecha) >= ANIO_INICIO, !is.na(inf_mensual))

inf_ts <- ts(ipc_analisis$inf_mensual,
             start = c(year(min(ipc_analisis$fecha)),
                       month(min(ipc_analisis$fecha))),
             frequency = 12)

cat("\nSerie completa: n =", length(inf_ts),
    "| Media =", round(mean(inf_ts), 4),
    "| DE =", round(sd(inf_ts), 4), "\n")

# ---- 3. Analisis preliminar ----

p_ipc <- ipc %>%
  filter(year(fecha) >= ANIO_INICIO) %>%
  ggplot(aes(fecha, IPC)) +
  geom_line(color = "#2c3e50", linewidth = 0.8) +
  geom_point(color = "#2c3e50", size = 1.2) +
  labs(title = "IPC — Bolivia (base 2017 = 100)", x = NULL, y = "IPC") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_ipc)

p_inf <- ipc_analisis %>%
  ggplot(aes(fecha, inf_mensual)) +
  geom_line(color = "#2980b9", linewidth = 0.7) +
  geom_point(color = "#2980b9", size = 1.2) +
  geom_hline(yintercept = mean(ipc_analisis$inf_mensual),
             linetype = "dashed", color = "#e74c3c") +
  labs(title = "Inflacion mensual (%)", x = NULL, y = "Var. mensual (%)") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_inf)

p_ia <- ipc %>%
  filter(year(fecha) >= ANIO_INICIO, !is.na(inf_interanual)) %>%
  ggplot(aes(fecha, inf_interanual)) +
  geom_line(color = "#27ae60", linewidth = 0.7) +
  geom_point(color = "#27ae60", size = 1.2) +
  labs(title = "Inflacion interanual (%)", x = NULL, y = "Var. 12 meses (%)") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_ia)

x_cent    <- as.numeric(inf_ts) - mean(inf_ts)
asimetria <- sum(x_cent^3) / length(inf_ts) / (sd(inf_ts)^3)
curtosis  <- sum(x_cent^4) / length(inf_ts) / (sd(inf_ts)^4)

cat("\n===== Estadisticas descriptivas — Inflacion mensual =====\n")
cat(sprintf("  Media        : %8.4f\n", mean(inf_ts)))
cat(sprintf("  Mediana      : %8.4f\n", median(inf_ts)))
cat(sprintf("  Desv. Est.   : %8.4f\n", sd(inf_ts)))
cat(sprintf("  Minimo       : %8.4f\n", min(inf_ts)))
cat(sprintf("  Maximo       : %8.4f\n", max(inf_ts)))
cat(sprintf("  Asimetria    : %8.4f\n", asimetria))
cat(sprintf("  Curtosis     : %8.4f\n", curtosis))
cat(sprintf("  Observaciones: %d\n",    length(inf_ts)))

p_hist <- data.frame(x = as.numeric(inf_ts)) %>%
  ggplot(aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15,
                 fill = "#3498db", alpha = 0.6, color = "white") +
  geom_density(color = "#e74c3c", linewidth = 1) +
  labs(title = "Distribucion de la inflacion mensual",
       x = "Var. mensual (%)", y = "Densidad")
print(p_hist)

# ---- 4. Pruebas de raiz unitaria ----

ipc_nivel_ts <- ts(ipc$IPC[ipc$Año >= ANIO_INICIO & !is.na(ipc$IPC)],
                   start = c(ANIO_INICIO, 1), frequency = 12)

cat("\n===== PRUEBAS DE RAIZ UNITARIA =====\n\n")

adf_nivel <- adf.test(ipc_nivel_ts)
adf_inf   <- adf.test(inf_ts)

cat("--- ADF (H0: raiz unitaria) ---\n")
cat(sprintf("  IPC nivel   : stat = %7.4f  p = %.4f  -> %s\n",
            adf_nivel$statistic, adf_nivel$p.value,
            ifelse(adf_nivel$p.value < 0.05, "ESTACIONARIA", "NO estacionaria")))
cat(sprintf("  Inf. mensual: stat = %7.4f  p = %.4f  -> %s\n",
            adf_inf$statistic, adf_inf$p.value,
            ifelse(adf_inf$p.value < 0.05, "ESTACIONARIA", "NO estacionaria")))

pp_nivel <- pp.test(ipc_nivel_ts)
pp_inf   <- pp.test(inf_ts)

cat("\n--- Phillips-Perron (H0: raiz unitaria) ---\n")
cat(sprintf("  IPC nivel   : stat = %7.4f  p = %.4f  -> %s\n",
            pp_nivel$statistic, pp_nivel$p.value,
            ifelse(pp_nivel$p.value < 0.05, "ESTACIONARIA", "NO estacionaria")))
cat(sprintf("  Inf. mensual: stat = %7.4f  p = %.4f  -> %s\n",
            pp_inf$statistic, pp_inf$p.value,
            ifelse(pp_inf$p.value < 0.05, "ESTACIONARIA", "NO estacionaria")))

kpss_nivel <- kpss.test(ipc_nivel_ts, null = "Level")
kpss_inf   <- kpss.test(inf_ts, null = "Level")

cat("\n--- KPSS (H0: estacionaria) ---\n")
cat(sprintf("  IPC nivel   : stat = %7.4f  p = %.4f  -> %s\n",
            kpss_nivel$statistic, kpss_nivel$p.value,
            ifelse(kpss_nivel$p.value < 0.05, "NO estacionaria", "ESTACIONARIA")))
cat(sprintf("  Inf. mensual: stat = %7.4f  p = %.4f  -> %s\n",
            kpss_inf$statistic, kpss_inf$p.value,
            ifelse(kpss_inf$p.value < 0.05, "NO estacionaria", "ESTACIONARIA")))

cat("\n--- ADF con seleccion de rezagos (urca) ---\n")
ur_nivel <- ur.df(ipc_nivel_ts, type = "trend", selectlags = "AIC")
ur_inf   <- ur.df(inf_ts, type = "drift", selectlags = "AIC")
cat("IPC nivel:\n"); print(summary(ur_nivel))
cat("Inf. mensual:\n"); print(summary(ur_inf))

d_recomendado <- ndiffs(inf_ts, test = "adf")
D_recomendado <- nsdiffs(inf_ts, test = "ocsb")
cat(sprintf("\nndiffs (ADF) = %d  |  nsdiffs (OCSB) = %d\n",
            d_recomendado, D_recomendado))

if (d_recomendado > 0) {
  inf_diff <- diff(inf_ts, differences = d_recomendado)
  cat("ADF sobre serie diferenciada: p =",
      round(adf.test(inf_diff)$p.value, 4), "\n")
  par(mfrow = c(1, 2))
  plot(inf_ts, main = "Serie original", ylab = "Inf. mensual (%)")
  plot(inf_diff, main = paste0("Diferenciada (d=", d_recomendado, ")"),
       ylab = "Inf. diferenciada")
  par(mfrow = c(1, 1))
}

# ---- 5. Identificacion (ACF / PACF) ----

cat("\n===== IDENTIFICACION DEL MODELO =====\n")

acf_vals  <- Acf(inf_ts, lag.max = 24, plot = FALSE)
pacf_vals <- Pacf(inf_ts, lag.max = 24, plot = FALSE)

par(mfrow = c(1, 2))
plot(acf_vals,  main = "ACF — Inflacion mensual")
plot(pacf_vals, main = "PACF — Inflacion mensual")
par(mfrow = c(1, 1))

lim_ic   <- qnorm(0.975) / sqrt(length(inf_ts))
acf_sig  <- which(abs(acf_vals$acf[-1]) > lim_ic)
pacf_sig <- which(abs(pacf_vals$acf)     > lim_ic)

cat("Banda de confianza 95%: +/-", round(lim_ic, 4), "\n")
cat("Rezagos significativos ACF :",
    if (length(acf_sig))  paste(acf_sig,  collapse = ", ") else "Ninguno", "\n")
cat("Rezagos significativos PACF:",
    if (length(pacf_sig)) paste(pacf_sig, collapse = ", ") else "Ninguno", "\n")

# ---- 6. Estimacion de modelos (serie completa) ----
#
# Dado que ndiffs=0 y las pruebas PP/KPSS/ADF-urca apoyan que la inflacion
# mensual es estacionaria, la seleccion principal se realiza sobre modelos
# con d=0.  Los modelos con d=1 se estiman aparte como analisis de
# sensibilidad y NO compiten en la seleccion del modelo ganador.

cat("\n===== ESTIMACION DE MODELOS ARIMA =====\n")

## Bloque principal: modelos estacionarios (d=0)
ordenes_principal <- list(
  c(1,0,0), c(2,0,0), c(3,0,0), c(4,0,0),
  c(0,0,1), c(0,0,2), c(0,0,3),
  c(1,0,1), c(2,0,1), c(1,0,2), c(2,0,2),
  c(3,0,1), c(1,0,3)
)

## Bloque de sensibilidad: modelos con diferenciación (d=1)
ordenes_sensibilidad <- list(
  c(0,1,1), c(1,1,0), c(1,1,1),
  c(2,1,0), c(2,1,1), c(0,1,2), c(1,1,2)
)

ajustar_arima <- function(serie, orden) {
  tryCatch({
    mod <- Arima(serie, order = orden, method = "ML")
    list(modelo = mod, orden = orden,
         aic = AIC(mod), bic = BIC(mod), aicc = mod$aicc,
         loglik = as.numeric(logLik(mod)),
         npar = length(coef(mod)) + 1, converged = TRUE)
  }, error = function(e) {
    list(modelo = NULL, orden = orden,
         aic = NA, bic = NA, aicc = NA,
         loglik = NA, npar = NA, converged = FALSE)
  })
}

diagnostico_modelo <- function(resultado) {
  if (!resultado$converged) return(NULL)
  mod      <- resultado$modelo
  res      <- residuals(mod)
  n_params <- length(coef(mod))
  ord      <- resultado$orden

  lb_pvals <- sapply(c(6, 12, 18), function(h) {
    tryCatch(Box.test(res, lag = h, type = "Ljung-Box",
                      fitdf = n_params)$p.value,
             error = function(e) NA)
  })
  jb_pval <- tryCatch(jarque.bera.test(res)$p.value, error = function(e) NA)

  ar_c <- mod$model$phi;  ma_c <- mod$model$theta
  raices_ok <- TRUE
  if (length(ar_c) > 0 && any(ar_c != 0))
    raices_ok <- raices_ok && all(Mod(1 / polyroot(c(1, -ar_c))) < 1)
  if (length(ma_c) > 0 && any(ma_c != 0))
    raices_ok <- raices_ok && all(Mod(1 / polyroot(c(1,  ma_c))) < 1)

  lb_ok <- all(lb_pvals > 0.05, na.rm = TRUE)

  data.frame(
    Modelo = sprintf("ARIMA(%d,%d,%d)", ord[1], ord[2], ord[3]),
    p = ord[1], d = ord[2], q = ord[3],
    AIC = round(resultado$aic, 2), BIC = round(resultado$bic, 2),
    LB_ok = ifelse(lb_ok, "Si", "No"),
    JB_pval = round(jb_pval, 3),
    Raices_ok = ifelse(raices_ok, "Si", "No"),
    Diagnostico = ifelse(lb_ok & raices_ok, "APROBADO", "NO APROBADO"),
    stringsAsFactors = FALSE
  )
}

## 6a. Seleccion principal (d=0)
cat("\n--- SELECCION PRINCIPAL (d=0, modelos estacionarios) ---\n")
cat("Candidatos:", length(ordenes_principal), "\n")

res_principal  <- lapply(ordenes_principal, function(o) ajustar_arima(inf_ts, o))
eval_principal <- do.call(rbind, lapply(res_principal, diagnostico_modelo))
eval_principal <- eval_principal[order(eval_principal$AIC), ]

cat("\n")
print(eval_principal %>% select(Modelo, AIC, BIC, LB_ok, JB_pval, Raices_ok, Diagnostico),
      row.names = FALSE)

aprobados <- eval_principal %>% filter(Diagnostico == "APROBADO")
if (nrow(aprobados) > 0) {
  mejor_nombre <- aprobados$Modelo[1]
  mejor_p <- aprobados$p[1]; mejor_d <- aprobados$d[1]; mejor_q <- aprobados$q[1]
} else {
  mejor_nombre <- eval_principal$Modelo[1]
  mejor_p <- eval_principal$p[1]; mejor_d <- eval_principal$d[1]; mejor_q <- eval_principal$q[1]
}
cat(sprintf("\n>>> MODELO GANADOR (d=0): %s\n", mejor_nombre))

## 6b. Sensibilidad (d=1)
cat("\n--- SENSIBILIDAD (d=1, modelos con diferenciacion) ---\n")
cat("Los modelos con d=1 modelan CAMBIOS en la inflacion mensual,\n")
cat("no la inflacion en nivel. Se reportan como referencia.\n")
cat("Candidatos:", length(ordenes_sensibilidad), "\n\n")

res_sensibilidad  <- lapply(ordenes_sensibilidad, function(o) ajustar_arima(inf_ts, o))
eval_sensibilidad <- do.call(rbind, lapply(res_sensibilidad, diagnostico_modelo))
eval_sensibilidad <- eval_sensibilidad[order(eval_sensibilidad$AIC), ]

print(eval_sensibilidad %>% select(Modelo, AIC, BIC, LB_ok, JB_pval, Raices_ok, Diagnostico),
      row.names = FALSE)

cat(sprintf("\n  Mejor d=1: %s (AIC=%.2f) — solo como referencia\n",
            eval_sensibilidad$Modelo[1], eval_sensibilidad$AIC[1]))

## 6c. auto.arima como benchmark

auto_fit   <- auto.arima(inf_ts, seasonal = TRUE, stepwise = FALSE,
                         approximation = FALSE, trace = FALSE)
auto_orden <- arimaorder(auto_fit)
cat(sprintf("\n  auto.arima: ARIMA(%d,%d,%d)  AIC=%.2f  (d=%d)\n",
            auto_orden[1], auto_orden[2], auto_orden[3], AIC(auto_fit), auto_orden[2]))


# ---- 6d. Sensibilidad con muestra comun de 40 observaciones ----
#
# Esta prueba recorta una observacion inicial de la serie completa.
# El objetivo es verificar si la seleccion del modelo cambia cuando todos
# los candidatos se estiman sobre una ventana comun mas corta.
#
# IMPORTANTE:
# - Esta NO reemplaza la seleccion principal.
# - Sirve como prueba de robustez.
# - El AIC/BIC de esta tabla debe compararse internamente dentro de esta muestra,
#   no contra los AIC/BIC de la muestra completa.

cat("\n--- SENSIBILIDAD: MUESTRA COMUN DE 40 OBSERVACIONES ---\n")

if (length(inf_ts) < 40) {
  stop("La serie tiene menos de 40 observaciones. No se puede hacer sensibilidad n=40.")
}

# Tomar las ultimas 40 observaciones disponibles
# En tu caso: probablemente feb-2023 a may-2026 si inf_ts tiene 41 obs.
inf_ts_40 <- tail(inf_ts, 40)

# Reconstruir ts preservando frecuencia mensual y fecha de inicio correcta
fecha_inicio_40 <- ipc_analisis$fecha[length(ipc_analisis$fecha) - 40 + 1]

inf_ts_40 <- ts(
  as.numeric(inf_ts_40),
  start = c(year(fecha_inicio_40), month(fecha_inicio_40)),
  frequency = 12
)

cat(sprintf("Serie sensibilidad n=40: inicio = %s | fin = %s | n = %d\n",
            format(fecha_inicio_40, "%b %Y"),
            format(max(ipc_analisis$fecha), "%b %Y"),
            length(inf_ts_40)))

# Todos los modelos: d=0 + d=1
ordenes_todos_40 <- c(ordenes_principal, ordenes_sensibilidad)

res_40 <- lapply(ordenes_todos_40, function(o) ajustar_arima(inf_ts_40, o))

eval_40 <- do.call(rbind, lapply(res_40, diagnostico_modelo))
eval_40 <- eval_40[order(eval_40$AIC), ]

cat("\nEvaluacion de TODOS los modelos sobre muestra comun n=40:\n")
print(eval_40 %>%
        select(Modelo, AIC, BIC, LB_ok, JB_pval, Raices_ok, Diagnostico),
      row.names = FALSE)

aprobados_40 <- eval_40 %>% filter(Diagnostico == "APROBADO")

if (nrow(aprobados_40) > 0) {
  mejor_40 <- aprobados_40[1, ]
} else {
  mejor_40 <- eval_40[1, ]
}

cat(sprintf("\n  Mejor modelo en sensibilidad n=40: %s (AIC=%.2f, BIC=%.2f)\n",
            mejor_40$Modelo, mejor_40$AIC, mejor_40$BIC))

cat(sprintf("  Modelo principal seleccionado con muestra completa: %s\n", mejor_nombre))

if (mejor_40$Modelo == mejor_nombre) {
  cat("  Resultado: ROBUSTO — el modelo principal tambien gana en n=40.\n")
} else {
  cat("  Resultado: SENSIBLE — el ganador cambia con n=40. Revisar con cautela.\n")
}

# ---- 6e. Sensibilidad con muestra compatible AR(4) (n=37) ----
#
# El AR(4) consume 4 observaciones iniciales para estimar.  Para que la
# comparacion de AIC/BIC sea estrictamente justa entre todos los modelos,
# se recorta la serie a las ultimas 37 observaciones — el numero efectivo
# que usa un AR(4) — y se re-estima todo.

cat("\n--- SENSIBILIDAD: MUESTRA COMPATIBLE AR(4) (n=37) ---\n")

N_COMPAT <- 37

if (length(inf_ts) < N_COMPAT) {
  cat("  Serie insuficiente para n=37. Saltando.\n")
} else {
  inf_ts_37 <- tail(inf_ts, N_COMPAT)
  fecha_inicio_37 <- ipc_analisis$fecha[length(ipc_analisis$fecha) - N_COMPAT + 1]
  inf_ts_37 <- ts(as.numeric(inf_ts_37),
                  start = c(year(fecha_inicio_37), month(fecha_inicio_37)),
                  frequency = 12)

  cat(sprintf("Serie sensibilidad n=%d: inicio = %s | fin = %s\n",
              N_COMPAT, format(fecha_inicio_37, "%b %Y"),
              format(max(ipc_analisis$fecha), "%b %Y")))

  ordenes_todos_37 <- c(ordenes_principal, ordenes_sensibilidad)
  res_37  <- lapply(ordenes_todos_37, function(o) ajustar_arima(inf_ts_37, o))
  eval_37 <- do.call(rbind, lapply(res_37, diagnostico_modelo))
  eval_37 <- eval_37[order(eval_37$AIC), ]

  cat("\nEvaluacion de TODOS los modelos sobre muestra n=37:\n")
  print(eval_37 %>%
          select(Modelo, AIC, BIC, LB_ok, JB_pval, Raices_ok, Diagnostico),
        row.names = FALSE)

  aprobados_37 <- eval_37 %>% filter(Diagnostico == "APROBADO")
  mejor_37 <- if (nrow(aprobados_37) > 0) aprobados_37[1, ] else eval_37[1, ]

  cat(sprintf("\n  Mejor modelo en n=%d: %s (AIC=%.2f, BIC=%.2f)\n",
              N_COMPAT, mejor_37$Modelo, mejor_37$AIC, mejor_37$BIC))
  cat(sprintf("  Modelo principal (muestra completa): %s\n", mejor_nombre))

  if (mejor_37$Modelo == mejor_nombre) {
    cat("  Resultado: ROBUSTO — el modelo principal tambien gana con informacion homogenea.\n")
  } else {
    cat("  Resultado: SENSIBLE — el ganador cambia con n=37. Revisar con cautela.\n")
  }
}

## Tabla combinada para referencia
eval_todos <- bind_rows(
  eval_principal %>% mutate(Bloque = "Principal d=0, muestra completa"),
  eval_sensibilidad %>% mutate(Bloque = "Sensibilidad d=1, muestra completa"),
  eval_40 %>% mutate(Bloque = "Sensibilidad todos, n=40"),
  eval_37 %>% mutate(Bloque = "Sensibilidad todos, n=37 (AR4-compatible)")
)

eval_todos <- eval_todos[order(eval_todos$Bloque, eval_todos$AIC), ]

# ---- 7. Diagnostico del modelo ganador ----

cat("\n===== DIAGNOSTICO DEL MODELO GANADOR =====\n")

modelo_ganador <- Arima(inf_ts, order = c(mejor_p, mejor_d, mejor_q), method = "ML")
residuos <- residuals(modelo_ganador)

cat("\n--- Coeficientes ---\n")
coefs   <- coef(modelo_ganador)
se_vals <- sqrt(diag(vcov(modelo_ganador)))
z_vals  <- coefs / se_vals
p_vals  <- 2 * (1 - pnorm(abs(z_vals)))
for (i in seq_along(coefs)) {
  sig <- ifelse(p_vals[i] < 0.01, "***",
         ifelse(p_vals[i] < 0.05, "**",
         ifelse(p_vals[i] < 0.10, "*", "   ")))
  cat(sprintf("  %-12s  est = %9.6f  se = %9.6f  z = %7.3f  p = %.4f %s\n",
              names(coefs)[i], coefs[i], se_vals[i], z_vals[i], p_vals[i], sig))
}

par(mfrow = c(2, 2))
plot(residuos, main = "Residuos del modelo", ylab = "Residuos", col = "#2c3e50")
abline(h = 0, col = "red", lty = 2)
hist(residuos, breaks = 15, probability = TRUE, main = "Histograma de residuos",
     col = "#3498db80", border = "white", xlab = "Residuos")
curve(dnorm(x, mean(residuos), sd(residuos)), add = TRUE, col = "#e74c3c", lwd = 2)
Acf(residuos,  main = "ACF de residuos",  lag.max = 20)
Pacf(residuos, main = "PACF de residuos", lag.max = 20)
par(mfrow = c(1, 1))

qqnorm(residuos, main = "Q-Q Plot de residuos", col = "#2c3e50", pch = 16)
qqline(residuos, col = "#e74c3c", lwd = 2)

cat("\n--- Ljung-Box ---\n")
for (h in c(6, 12, 18, 24)) {
  lb <- Box.test(residuos, lag = h, type = "Ljung-Box",
                 fitdf = length(coef(modelo_ganador)))
  cat(sprintf("  h=%2d : Q=%8.4f  p=%.4f  -> %s\n",
              h, lb$statistic, lb$p.value,
              ifelse(lb$p.value > 0.05, "Ruido blanco", "AUTOCORRELACION")))
}

jb <- jarque.bera.test(residuos)
sw <- shapiro.test(residuos)
cat(sprintf("\n--- Normalidad ---\n  Jarque-Bera : p=%.4f -> %s\n",
            jb$p.value, ifelse(jb$p.value > 0.05, "Normal", "NO normal")))
cat(sprintf("  Shapiro-Wilk: p=%.4f -> %s\n",
            sw$p.value, ifelse(sw$p.value > 0.05, "Normal", "NO normal")))

ar_coefs <- modelo_ganador$model$phi
ma_coefs <- modelo_ganador$model$theta
tiene_ar <- length(ar_coefs) > 0 && any(ar_coefs != 0)
tiene_ma <- length(ma_coefs) > 0 && any(ma_coefs != 0)
if (tiene_ar || tiene_ma) {
  theta_c <- seq(0, 2 * pi, length.out = 200)
  n_plots <- tiene_ar + tiene_ma
  if (n_plots == 2) par(mfrow = c(1, 2))
  if (tiene_ar) {
    ar_inv <- 1 / polyroot(c(1, -ar_coefs))
    plot(Re(ar_inv), Im(ar_inv), xlim = c(-1.5,1.5), ylim = c(-1.5,1.5),
         pch = 19, col = "#e74c3c", cex = 2, asp = 1,
         main = "Raices inversas AR", xlab = "Real", ylab = "Imaginario")
    lines(cos(theta_c), sin(theta_c), col = "grey50", lwd = 2)
    abline(h = 0, v = 0, lty = 3, col = "grey70")
    cat("\nRaices inversas AR:", round(Mod(ar_inv), 4), "-> OK:", all(Mod(ar_inv) < 1), "\n")
  }
  if (tiene_ma) {
    ma_inv <- 1 / polyroot(c(1, ma_coefs))
    plot(Re(ma_inv), Im(ma_inv), xlim = c(-1.5,1.5), ylim = c(-1.5,1.5),
         pch = 19, col = "#2980b9", cex = 2, asp = 1,
         main = "Raices inversas MA", xlab = "Real", ylab = "Imaginario")
    lines(cos(theta_c), sin(theta_c), col = "grey50", lwd = 2)
    abline(h = 0, v = 0, lty = 3, col = "grey70")
    cat("Raices inversas MA:", round(Mod(ma_inv), 4), "-> OK:", all(Mod(ma_inv) < 1), "\n")
  }
  par(mfrow = c(1, 1))
}

checkresiduals(modelo_ganador)

# ============================================================================
# ---- 8. BACKTESTING (ultimos 6 meses) ----
# ============================================================================

cat("\n")
cat("##########################################################\n")
cat("#              BACKTESTING — ULTIMOS 6 MESES             #\n")
cat("##########################################################\n")

## 8a. Separar entrenamiento / prueba
n_total <- length(inf_ts)
n_train <- n_total - H_BACKTEST

train_ts <- window(inf_ts, end = time(inf_ts)[n_train])
test_ts  <- window(inf_ts, start = time(inf_ts)[n_train + 1])

fechas_test  <- tail(ipc_analisis$fecha, H_BACKTEST)
fechas_train <- head(ipc_analisis$fecha, n_train)

cat(sprintf("\n  Entrenamiento : %s a %s  (n = %d)\n",
            format(min(fechas_train), "%b %Y"),
            format(max(fechas_train), "%b %Y"), length(train_ts)))
cat(sprintf("  Prueba        : %s a %s  (n = %d)\n",
            format(min(fechas_test), "%b %Y"),
            format(max(fechas_test), "%b %Y"), length(test_ts)))

## 8b. Estimar modelo en datos de entrenamiento
modelo_bt <- Arima(train_ts, order = c(mejor_p, mejor_d, mejor_q), method = "ML")
fc_bt     <- forecast(modelo_bt, h = H_BACKTEST, level = c(80, 95))

## 8c. Comparacion punto a punto
actual_vals   <- as.numeric(test_ts)
forecast_vals <- as.numeric(fc_bt$mean)
error_vals    <- forecast_vals - actual_vals
pct_error     <- (error_vals / abs(actual_vals)) * 100
lo95 <- as.numeric(fc_bt$lower[, 2])
hi95 <- as.numeric(fc_bt$upper[, 2])
dentro_ic95 <- actual_vals >= lo95 & actual_vals <= hi95

cat("\n--- Comparacion punto a punto ---\n")
cat(sprintf("%-15s %10s %10s %10s %10s %6s\n",
            "Fecha", "Observado", "Pronostic", "Error", "Error(%)", "IC95?"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (i in seq_along(actual_vals)) {
  cat(sprintf("%-15s %10.4f %10.4f %+10.4f %+9.2f%%   %s\n",
              format(fechas_test[i], "%b %Y"),
              actual_vals[i], forecast_vals[i],
              error_vals[i], pct_error[i],
              ifelse(dentro_ic95[i], "Si", "NO")))
}

## 8d. Metricas de error
rmse <- sqrt(mean(error_vals^2))
mae  <- mean(abs(error_vals))
mape <- mean(abs(pct_error))
me   <- mean(error_vals)
pct_dentro <- mean(dentro_ic95) * 100

cat(sprintf("\n--- Metricas de error ---\n"))
cat(sprintf("  RMSE                        : %8.4f pp\n", rmse))
cat(sprintf("  MAE                         : %8.4f pp\n", mae))
cat(sprintf("  MAPE                        : %8.2f%%\n",  mape))
cat(sprintf("  Error medio (sesgo)         : %+8.4f pp\n", me))
cat(sprintf("  Obs. dentro del IC 95%%      : %d/%d (%.0f%%)\n",
            sum(dentro_ic95), H_BACKTEST, pct_dentro))

bt_aceptable <- pct_dentro >= 50 && rmse < 2 * sd(train_ts)
cat(sprintf("\n  >>> Backtesting: %s\n",
            ifelse(bt_aceptable,
                   "ACEPTABLE — el modelo captura razonablemente la dinamica",
                   "ADVERTENCIA — revisar especificacion del modelo")))

## 8e. Grafico 1: Backtest — pronostico vs observado
bt_df <- data.frame(
  fecha     = rep(fechas_test, 2),
  valor     = c(actual_vals, forecast_vals),
  tipo      = rep(c("Observado", "Pronostico backtest"), each = H_BACKTEST)
)

bt_ci <- data.frame(
  fecha = fechas_test,
  lo80  = as.numeric(fc_bt$lower[, 1]),
  hi80  = as.numeric(fc_bt$upper[, 1]),
  lo95  = lo95,
  hi95  = hi95
)

p_bt1 <- ggplot() +
  geom_ribbon(data = bt_ci, aes(x = fecha, ymin = lo95, ymax = hi95),
              fill = "#e74c3c", alpha = 0.10) +
  geom_ribbon(data = bt_ci, aes(x = fecha, ymin = lo80, ymax = hi80),
              fill = "#e74c3c", alpha = 0.20) +
  geom_line(data = data.frame(fecha = fechas_train,
                               valor = as.numeric(train_ts)),
            aes(fecha, valor), color = "#95a5a6", linewidth = 0.5) +
  geom_line(data = bt_df, aes(fecha, valor, color = tipo, linetype = tipo),
            linewidth = 0.9) +
  geom_point(data = bt_df, aes(fecha, valor, color = tipo), size = 2.5) +
  geom_segment(data = data.frame(fecha = fechas_test,
                                  obs = actual_vals, pron = forecast_vals),
               aes(x = fecha, xend = fecha, y = obs, yend = pron),
               color = "#7f8c8d", linetype = "dotted", linewidth = 0.6) +
  scale_color_manual(values = c("Observado" = "#2c3e50",
                                "Pronostico backtest" = "#e74c3c")) +
  scale_linetype_manual(values = c("Observado" = "solid",
                                   "Pronostico backtest" = "dashed")) +
  labs(title = sprintf("Backtesting — %s (ultimos %d meses)", mejor_nombre, H_BACKTEST),
       subtitle = sprintf("RMSE=%.4f | MAE=%.4f | MAPE=%.1f%% | Dentro IC95: %.0f%%",
                          rmse, mae, mape, pct_dentro),
       x = NULL, y = "Inflacion mensual (%)", color = NULL, linetype = NULL) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
print(p_bt1)

## 8f. Grafico 2: Barras de error del backtest
p_bt2 <- data.frame(fecha = fechas_test, error = error_vals) %>%
  ggplot(aes(x = fecha, y = error, fill = error > 0)) +
  geom_col(width = 20, alpha = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = c(-rmse, rmse), linetype = "dashed", color = "#e74c3c") +
  annotate("text", x = min(fechas_test), y = rmse + 0.05,
           label = paste0("RMSE = ", round(rmse, 3)),
           hjust = 0, size = 3.5, color = "#e74c3c") +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2980b9")) +
  labs(title = "Error del backtest (pronostico - observado)",
       x = NULL, y = "Error (pp)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_bt2)

# ============================================================================
# ---- 9. Pronostico final (solo si el backtest es aceptable) ----
# ============================================================================

cat("\n")
cat("##########################################################\n")
cat("#                 PRONOSTICO FINAL                       #\n")
cat("##########################################################\n")

if (!bt_aceptable) {
  cat("\n  ADVERTENCIA: El backtest no fue del todo satisfactorio.\n")
  cat("  Se procede con la proyeccion, pero interpretar con cautela.\n")
}

pronostico <- forecast(modelo_ganador, h = H_FORECAST, level = c(80, 95))

ultimo_mes <- max(ipc_analisis$fecha)
meses_pron <- seq.Date(ultimo_mes %m+% months(1), by = "month", length.out = H_FORECAST)

cat("\n--- Pronostico de inflacion mensual ---\n")
cat(sprintf("%-15s %8s %10s %10s %10s %10s\n",
            "Fecha", "Punto", "IC80_inf", "IC80_sup", "IC95_inf", "IC95_sup"))
for (i in 1:H_FORECAST) {
  cat(sprintf("%-15s %8.4f %10.4f %10.4f %10.4f %10.4f\n",
              format(meses_pron[i], "%B %Y"),
              pronostico$mean[i],
              pronostico$lower[i, 1], pronostico$upper[i, 1],
              pronostico$lower[i, 2], pronostico$upper[i, 2]))
}

ultimo_ipc <- tail(ipc_analisis$IPC, 1)
ipc_proy   <- numeric(H_FORECAST)
ipc_proy[1] <- ultimo_ipc * (1 + pronostico$mean[1] / 100)
for (i in 2:H_FORECAST) {
  ipc_proy[i] <- ipc_proy[i - 1] * (1 + pronostico$mean[i] / 100)
}

ipc_completo    <- c(ipc$IPC[!is.na(ipc$IPC)], ipc_proy)
fechas_completo <- c(ipc$fecha[!is.na(ipc$IPC)], meses_pron)

inf_interanual_proy <- sapply(seq_along(meses_pron), function(i) {
  fecha_12 <- meses_pron[i] %m-% months(12)
  idx_12   <- which.min(abs(fechas_completo - fecha_12))
  (ipc_proy[i] / ipc_completo[idx_12] - 1) * 100
})

cat("\n--- IPC e inflacion proyectados ---\n")
cat(sprintf("%-15s %12s %12s %16s\n",
            "Fecha", "IPC proy.", "Inf.mens.(%)", "Inf.interan.(%)"))
for (i in 1:H_FORECAST) {
  cat(sprintf("%-15s %12.4f %12.4f %16.4f\n",
              format(meses_pron[i], "%B %Y"),
              ipc_proy[i], as.numeric(pronostico$mean[i]),
              inf_interanual_proy[i]))
}

# ---- 10. Graficos finales ----

## 10a. Grafico integrado: historico + backtest + proyeccion
historico_df <- ipc_analisis %>%
  select(fecha, inf_mensual) %>%
  mutate(tipo = "Historico")

backtest_pron_df <- data.frame(
  fecha       = fechas_test,
  inf_mensual = forecast_vals,
  tipo        = "Backtest (pronostico)"
)

backtest_obs_df <- data.frame(
  fecha       = fechas_test,
  inf_mensual = actual_vals,
  tipo        = "Backtest (observado)"
)

proyeccion_df <- data.frame(
  fecha       = meses_pron,
  inf_mensual = as.numeric(pronostico$mean),
  tipo        = "Proyeccion"
)

ci_proy <- data.frame(
  fecha = meses_pron,
  lo95  = as.numeric(pronostico$lower[, 2]),
  hi95  = as.numeric(pronostico$upper[, 2])
)

plot_todo <- bind_rows(historico_df, backtest_obs_df, backtest_pron_df, proyeccion_df) %>%
  mutate(tipo = factor(tipo, levels = c("Historico", "Backtest (observado)",
                                         "Backtest (pronostico)", "Proyeccion")))

colores <- c("Historico"              = "#2c3e50",
             "Backtest (observado)"   = "#27ae60",
             "Backtest (pronostico)"  = "#e67e22",
             "Proyeccion"             = "#e74c3c")
lineas  <- c("Historico"              = "solid",
             "Backtest (observado)"   = "solid",
             "Backtest (pronostico)"  = "dashed",
             "Proyeccion"             = "dashed")

p_final <- ggplot() +
  geom_ribbon(data = ci_proy, aes(x = fecha, ymin = lo95, ymax = hi95),
              fill = "#e74c3c", alpha = 0.12) +
  geom_vline(xintercept = min(fechas_test), linetype = "dotted",
             color = "#7f8c8d", linewidth = 0.5) +
  geom_vline(xintercept = max(ipc_analisis$fecha) + 15, linetype = "dotted",
             color = "#7f8c8d", linewidth = 0.5) +
  annotate("text", x = min(fechas_test) - 20,
           y = max(plot_todo$inf_mensual) * 0.95,
           label = "Entrenamiento", hjust = 1, size = 3, color = "#7f8c8d") +
  annotate("text", x = min(fechas_test) + 45,
           y = max(plot_todo$inf_mensual) * 0.95,
           label = "Backtest", hjust = 0.5, size = 3, color = "#7f8c8d") +
  annotate("text", x = max(ipc_analisis$fecha) + 60,
           y = max(plot_todo$inf_mensual) * 0.95,
           label = "Proyeccion", hjust = 0, size = 3, color = "#7f8c8d") +
  geom_line(data = plot_todo, aes(fecha, inf_mensual, color = tipo, linetype = tipo),
            linewidth = 0.8) +
  geom_point(data = plot_todo, aes(fecha, inf_mensual, color = tipo), size = 1.8) +
  geom_segment(data = data.frame(fecha = fechas_test,
                                  obs = actual_vals, pron = forecast_vals),
               aes(x = fecha, xend = fecha, y = obs, yend = pron),
               color = "#7f8c8d", linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = lineas) +
  labs(title = sprintf("Inflacion mensual — %s", mejor_nombre),
       subtitle = sprintf("Historico + Backtest (RMSE=%.3f, MAPE=%.1f%%) + Proyeccion a dic 2026",
                          rmse, mape),
       x = NULL, y = "Variacion mensual (%)", color = NULL, linetype = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
print(p_final)

## 10b. IPC observado + proyectado
ipc_obs_df <- ipc %>%
  filter(year(fecha) >= ANIO_INICIO, !is.na(IPC)) %>%
  select(fecha, IPC) %>%
  mutate(tipo = "Historico")
ipc_proy_df <- data.frame(fecha = meses_pron, IPC = ipc_proy, tipo = "Proyectado")
ipc_all <- bind_rows(ipc_obs_df, ipc_proy_df)

p_ipc_final <- ggplot(ipc_all, aes(fecha, IPC, color = tipo, linetype = tipo)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  geom_vline(xintercept = max(ipc_analisis$fecha) + 15,
             linetype = "dotted", color = "#7f8c8d") +
  scale_color_manual(values = c(Historico = "#2c3e50", Proyectado = "#e74c3c")) +
  scale_linetype_manual(values = c(Historico = "solid", Proyectado = "dashed")) +
  labs(title = "IPC historico y proyectado", subtitle = "Base 2017 = 100",
       x = NULL, y = "IPC", color = NULL, linetype = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
print(p_ipc_final)

# ---- 11. Exportar resultados (archivo separado) ----

cat("\n===== EXPORTACION DE RESULTADOS =====\n")

meses_nombres <- c("Enero","Febrero","Marzo","Abril","Mayo","Junio",
                   "Julio","Agosto","Septiembre","Octubre","Noviembre","Diciembre")

## Hoja 1: Proyecciones
hoja_proy <- data.frame(
  Anio  = year(meses_pron),
  Mes   = meses_nombres[month(meses_pron)],
  IPC_proyectado         = round(ipc_proy, 4),
  Inf_mensual_pct        = round(as.numeric(pronostico$mean), 4),
  Inf_interanual_pct     = round(inf_interanual_proy, 4),
  IC95_inf               = round(as.numeric(pronostico$lower[, 2]), 4),
  IC95_sup               = round(as.numeric(pronostico$upper[, 2]), 4),
  Tipo                   = "PROYECCION"
)

## Hoja 2: Backtesting
hoja_bt <- data.frame(
  Anio      = year(fechas_test),
  Mes       = meses_nombres[month(fechas_test)],
  Observado = round(actual_vals, 4),
  Pronostico = round(forecast_vals, 4),
  Error      = round(error_vals, 4),
  Error_pct  = round(pct_error, 2),
  IC95_inf   = round(lo95, 4),
  IC95_sup   = round(hi95, 4),
  Dentro_IC95 = ifelse(dentro_ic95, "Si", "No")
)

## Hoja 3: Metricas
hoja_metricas <- data.frame(
  Metrica = c("Modelo", "AIC", "BIC",
              "RMSE (backtest)", "MAE (backtest)", "MAPE (backtest)",
              "Error medio (backtest)", "Obs. dentro IC95",
              "Ljung-Box", "Jarque-Bera"),
  Valor = c(
    mejor_nombre,
    round(AIC(modelo_ganador), 2),
    round(BIC(modelo_ganador), 2),
    round(rmse, 4),
    round(mae, 4),
    paste0(round(mape, 2), "%"),
    round(me, 4),
    paste0(sum(dentro_ic95), "/", H_BACKTEST),
    ifelse(all(sapply(c(6,12,18), function(h)
      Box.test(residuos, lag=h, type="Ljung-Box",
               fitdf=length(coef(modelo_ganador)))$p.value > 0.05)),
      "APROBADO", "NO aprobado"),
    ifelse(jb$p.value > 0.05, "APROBADO",
           paste0("NO aprobado (p=", round(jb$p.value, 4), ")"))
  )
)

## Hoja 4: Datos historicos (base completa para referencia)
hoja_hist <- ipc_analisis %>%
  select(fecha, Año, Mes, IPC, inf_mensual, inf_interanual) %>%
  mutate(Tipo = "HISTORICO",
         inf_mensual    = round(inf_mensual, 4),
         inf_interanual = round(inf_interanual, 4),
         IPC            = round(IPC, 4))

wb_out <- createWorkbook()
addWorksheet(wb_out, "Proyecciones")
addWorksheet(wb_out, "Backtesting")
addWorksheet(wb_out, "Metricas")
addWorksheet(wb_out, "Historico")

writeData(wb_out, "Proyecciones", hoja_proy)
writeData(wb_out, "Backtesting",  hoja_bt)
writeData(wb_out, "Metricas",     hoja_metricas)
writeData(wb_out, "Historico",    hoja_hist)

saveWorkbook(wb_out, RUTA_OUTPUT, overwrite = TRUE)
cat(sprintf("Resultados exportados a: %s\n", RUTA_OUTPUT))
cat("  Hoja 'Proyecciones' : IPC e inflacion proyectados (jun-dic 2026)\n")
cat("  Hoja 'Backtesting'  : Comparacion obs vs pronostico (ultimos 6 meses)\n")
cat("  Hoja 'Metricas'     : Indicadores de bondad del modelo\n")
cat("  Hoja 'Historico'    : Serie historica completa\n")

# ---- 12. Resumen final ----

cat("\n")
cat("==========================================================\n")
cat("                    RESUMEN FINAL\n")
cat("==========================================================\n")
cat(sprintf("  Modelo ganador  : %s\n", mejor_nombre))
cat(sprintf("  AIC / BIC       : %.2f / %.2f\n", AIC(modelo_ganador), BIC(modelo_ganador)))
cat(sprintf("  Ljung-Box       : %s\n",
            ifelse(all(sapply(c(6,12,18), function(h)
              Box.test(residuos, lag=h, type="Ljung-Box",
                       fitdf=length(coef(modelo_ganador)))$p.value > 0.05)),
              "APROBADO", "NO aprobado")))
cat(sprintf("  Jarque-Bera     : p=%.4f\n", jb$p.value))
cat("  ---\n")
cat(sprintf("  Backtest RMSE   : %.4f pp\n", rmse))
cat(sprintf("  Backtest MAE    : %.4f pp\n", mae))
cat(sprintf("  Backtest MAPE   : %.2f%%\n",  mape))
cat(sprintf("  Dentro IC95     : %d/%d\n", sum(dentro_ic95), H_BACKTEST))
cat("  ---\n")
cat(sprintf("  Inf. mensual promedio proyectada : %.4f%%\n", mean(pronostico$mean)))
cat(sprintf("  IPC dic 2026 proyectado          : %.4f\n", tail(ipc_proy, 1)))
if (length(inf_interanual_proy) > 0)
  cat(sprintf("  Inf. interanual dic 2026         : %.2f%%\n",
              tail(inf_interanual_proy, 1)))
cat("==========================================================\n")
cat("BD Bruta.xlsx NO fue modificada.\n")
cat(sprintf("Resultados en: %s\n", RUTA_OUTPUT))
cat("Fin del script.\n")
