# ============================================================================
# ARIMA.R — Funciones auxiliares para modelización ARIMA de inflación
# Utilizado por: ARIMA_Inflacion.Rmd (documento knitr)
# ============================================================================

# Para generar el documento completo, ejecutar en R:
#   rmarkdown::render("ARIMA_Inflacion.Rmd")

# --- Instalación de paquetes (ejecutar una vez) ---
instalar_paquetes <- function() {
  required <- c("readxl","dplyr","lubridate","ggplot2","forecast",
                 "tseries","urca","lmtest","kableExtra","scales",
                 "gridExtra","rmarkdown","knitr")
  new_pkg  <- required[!(required %in% installed.packages()[,"Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg, repos = "https://cran.r-project.org",
                     dependencies = TRUE)
  }
  cat("Todos los paquetes instalados.\n")
}

# --- Carga de datos IPC ---
cargar_ipc <- function(path = "Inputs/BD Bruta.xlsx") {
  library(readxl)
  library(dplyr)

  meses_map <- c("Enero"=1,"Febrero"=2,"Marzo"=3,"Abril"=4,"Mayo"=5,
                 "Junio"=6,"Julio"=7,"Agosto"=8,"Septiembre"=9,
                 "Octubre"=10,"Noviembre"=11,"Diciembre"=12)

  ipc_raw <- read_excel(path, sheet = "IPC")

  ipc <- ipc_raw %>%
    mutate(Valor_num = suppressWarnings(as.numeric(Valor))) %>%
    filter(!is.na(Valor_num)) %>%
    mutate(
      mes_num = meses_map[Mes],
      fecha   = as.Date(paste(Año, mes_num, "01", sep = "-")),
      IPC     = Valor_num
    ) %>%
    arrange(fecha) %>%
    select(fecha, Año, Mes, mes_num, IPC) %>%
    mutate(
      inf_mensual    = (IPC / lag(IPC)     - 1) * 100,
      inf_interanual = (IPC / lag(IPC, 12) - 1) * 100
    )

  return(ipc)
}

# --- Carga y mensualización de TC ---
cargar_tc_mensual <- function(path = "Inputs/BD Bruta.xlsx") {
  library(readxl)
  library(dplyr)
  library(lubridate)

  tc_raw <- read_excel(path, sheet = "TC")
  names(tc_raw) <- trimws(names(tc_raw))

  tc_mensual <- tc_raw %>%
    mutate(
      fecha_date = as.Date(fecha),
      año_tc     = year(fecha_date),
      mes_tc     = month(fecha_date)
    ) %>%
    group_by(año_tc, mes_tc) %>%
    summarise(
      TC_preferencial = mean(`TC PREFERENCIAL`, na.rm = TRUE),
      TC_digital      = mean(`TC-DIGITAL`,      na.rm = TRUE),
      TC_oficial      = first(`TC-OFICIAL`),
      TC_referencial  = mean(`TC-REFERENCIAL`,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(fecha = as.Date(paste(año_tc, mes_tc, "01", sep = "-")))

  return(tc_mensual)
}

# --- Ajustar modelo ARIMA con manejo de errores ---
ajustar_arima <- function(serie, orden, metodo = "ML") {
  tryCatch({
    modelo <- forecast::Arima(serie, order = orden, method = metodo)
    list(
      modelo    = modelo,
      orden     = orden,
      aic       = AIC(modelo),
      bic       = BIC(modelo),
      aicc      = modelo$aicc,
      loglik    = as.numeric(logLik(modelo)),
      converged = TRUE
    )
  }, error = function(e) {
    list(
      modelo    = NULL,
      orden     = orden,
      aic       = NA, bic = NA, aicc = NA, loglik = NA,
      converged = FALSE
    )
  })
}

# --- Evaluar diagnóstico de un modelo ---
diagnostico_arima <- function(modelo, rezagos = c(6, 12, 18)) {
  res      <- residuals(modelo)
  n_params <- length(coef(modelo))

  lb_pvals <- sapply(rezagos, function(h) {
    tryCatch(
      Box.test(res, lag = h, type = "Ljung-Box", fitdf = n_params)$p.value,
      error = function(e) NA
    )
  })

  jb <- tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL)

  ar_c <- modelo$model$phi
  ma_c <- modelo$model$theta
  raices_ok <- TRUE
  if (length(ar_c) > 0 && any(ar_c != 0))
    raices_ok <- raices_ok && all(Mod(1/polyroot(c(1, -ar_c))) < 1)
  if (length(ma_c) > 0 && any(ma_c != 0))
    raices_ok <- raices_ok && all(Mod(1/polyroot(c(1, ma_c))) < 1)

  list(
    lb_pvals   = lb_pvals,
    lb_ok      = all(lb_pvals > 0.05, na.rm = TRUE),
    jb_pval    = if (!is.null(jb)) jb$p.value else NA,
    jb_ok      = if (!is.null(jb)) jb$p.value > 0.05 else NA,
    raices_ok  = raices_ok
  )
}

# --- Renderizar el documento ---
renderizar <- function() {
  rmarkdown::render("ARIMA_Inflacion.Rmd", output_format = "html_document")
  cat("Documento generado: ARIMA_Inflacion.html\n")
}
