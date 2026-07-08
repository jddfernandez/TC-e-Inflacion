%% ARIMA_M2.m — Modelo Logaritmico (M2) de Inflacion en Bolivia
%
% Replica en MATLAB del modelo ganador M2 (R): ARIMA(p,d,q) sobre log(IPC).
% A diferencia de code_P2.m / OLS.m (identificacion Hannan-Rissanen via
% MCO manual), aqui se usa Econometrics Toolbox (estimate/forecast/infer)
% para estimacion por maxima verosimilitud, manteniendo el mismo espiritu
% metodologico: grilla sistematica + seleccion por criterios de informacion
% (AICc, BIC, HQC) combinados con precision (RMSE, MAPE) y diagnostico de
% residuos (Ljung-Box, raices invertidas).
%
% Requiere: matlab workspace ARIMA.mat (variables BDBruta_IPC, BDBruta_TC)
% en la misma carpeta que este script.

clc; clear; close all;
PROJ_ROOT = 'C:\Users\Juande\Documents\Scripts Python\TC e π';
scriptDir = fullfile(PROJ_ROOT, 'matlab scripts');
cd(scriptDir);

%% 0. Configuracion
ANIO_INICIO = 2023;
H_BACKTEST  = 6;
H_FORECAST  = 7;
MAX_P = 4;
MAX_Q = 4;
RUTA_OUTPUT = fullfile(PROJ_ROOT, 'Outputs', 'Proyecciones_MATLAB.xlsx');
addpath 'C:\Users\Juande\Documents\Scripts Python\TC e π\matlab scripts';

%% 1. Carga de datos (workspace pre-importado)
load('matlab workspace ARIMA.mat');

ipc_tbl = sortrows(BDBruta_IPC(BDBruta_IPC.A_o >= ANIO_INICIO, :), 'Time');
fechas  = ipc_tbl.Time;
IPC     = ipc_tbl.Valor;
logIPC  = log(IPC);
n       = length(logIPC);

fprintf('Serie M2 cargada: %s a %s (n=%d)\n', ...
    datestr(fechas(1),'mmm yyyy'), datestr(fechas(end),'mmm yyyy'), n);

%% 2. Analisis preliminar

dlogIPC = diff(logIPC);

figure('Name','Analisis preliminar','Position',[100 100 900 700]);
subplot(2,2,1); plot(fechas, IPC, '-o', 'Color',[0.17 0.24 0.31]);
title('IPC (base 2017=100)'); grid on;
subplot(2,2,2); plot(fechas, logIPC, '-o', 'Color',[0.56 0.27 0.68]);
title('log(IPC)'); grid on;
subplot(2,2,3); plot(fechas(2:end), dlogIPC, '-o', 'Color',[0.16 0.50 0.73]);
yline(0,'r--'); title('\Delta log(IPC)'); grid on;
subplot(2,2,4); histogram(dlogIPC, 12, 'FaceColor',[0.20 0.60 0.86]);
title('Distribucion \Delta log(IPC)'); grid on;

fprintf('\n--- Estadisticas Dlog(IPC) ---\n');
fprintf('  Media   : %.6f (~%.4f%% mensual)\n', mean(dlogIPC), mean(dlogIPC)*100);
fprintf('  DE      : %.6f\n', std(dlogIPC));
fprintf('  Asimetria: %.4f\n', skewness(dlogIPC));
fprintf('  Curtosis : %.4f\n', kurtosis(dlogIPC));

%% 3. Pruebas de raiz unitaria

fprintf('\n===== PRUEBAS DE RAIZ UNITARIA =====\n');
fprintf('  H0 (ADF / PP / ZA) : serie tiene raiz unitaria\n');
fprintf('  H0 (KPSS)          : serie es estacionaria\n');
fprintf('  ZA: H1 = estacionaria con un quiebre estructural unico en fecha desconocida\n\n');

% --- ADF, Phillips-Perron, KPSS ---
[~, p_adf_log]   = adftest(logIPC);
[~, p_pp_log]    = pptest(logIPC);
[~, p_kpss_log]  = kpsstest(logIPC);

[~, p_adf_dlog]  = adftest(dlogIPC);
[~, p_pp_dlog]   = pptest(dlogIPC);
[~, p_kpss_dlog] = kpsstest(dlogIPC);

% --- Zivot-Andrews (zatest.m) ---
% Modelo 'ARD' = quiebre en intercepto y tendencia (caso mas general).

za_ok = false;

try
    [h_za_log,  p_za_log,  stat_za_log,  cv_za_log,  bd_log]  = zatest(logIPC,  'Model','ARD');
    [h_za_dlog, p_za_dlog, stat_za_dlog, cv_za_dlog, bd_dlog] = zatest(dlogIPC, 'Model','ARD');

    za_ok = true;

catch ME
    warning('zatest no pudo calcularse: %s', ME.message);

    p_za_log  = NaN;
    p_za_dlog = NaN;

    bd_log  = NaN;
    bd_dlog = NaN;
end

% --- Tabla de resultados ---
fprintf('  %-24s  log(IPC)                     Dlog(IPC)\n', 'Prueba');
fprintf('  %s\n', repmat('-',1,82));

fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n', 'ADF', ...
    p_adf_log,  ternary(p_adf_log<0.05,  '[Estacionaria]  ', '[Raiz unitaria]'), ...
    p_adf_dlog, ternary(p_adf_dlog<0.05, '[Estacionaria]',   '[Raiz unitaria]'));

fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n', 'Phillips-Perron', ...
    p_pp_log,  ternary(p_pp_log<0.05,  '[Estacionaria]  ', '[Raiz unitaria]'), ...
    p_pp_dlog, ternary(p_pp_dlog<0.05, '[Estacionaria]',   '[Raiz unitaria]'));

fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n', 'KPSS (H0: estac.)', ...
    p_kpss_log,  ternary(p_kpss_log<0.05,  '[Raiz unitaria] ', '[Estacionaria] '), ...
    p_kpss_dlog, ternary(p_kpss_dlog<0.05, '[Raiz unitaria]',  '[Estacionaria]'));

if za_ok
    fprintf('  %-24s  tau=%.4f p≈%.4f %-14s  tau=%.4f p≈%.4f %s\n', ...
        'Zivot-Andrews (*)', ...
        stat_za_log, ...
        p_za_log, ...
        ternary(h_za_log, '[Estac.+quiebre]', '[Raiz unitaria]'), ...
        stat_za_dlog, ...
        p_za_dlog, ...
        ternary(h_za_dlog, '[Estac.+quiebre]', '[Raiz unitaria]'));

    % Mapeo de fechas de quiebre
    % dlogIPC(i) = logIPC(i+1)-logIPC(i)
    % Entonces, el quiebre de dlogIPC en i se reporta como fechas(i+1).

    try
        fd_log  = fechas(bd_log);
        fd_dlog = fechas(min(bd_dlog + 1, length(fechas)));

        fprintf('  Fecha de quiebre: log(IPC) = %-10s  |  Dlog(IPC) = %s\n', ...
            datestr(fd_log,'mmm yyyy'), datestr(fd_dlog,'mmm yyyy'));

        fprintf('  CV 5%% (modelo ARD): log(IPC)=%.3f | Dlog(IPC)=%.3f\n', ...
            cv_za_log, cv_za_dlog);
    catch
    end

    fprintf('  (*) Rechazo implica estacionariedad condicionada al quiebre detectado\n');
    fprintf('      p≈ es una aproximacion; la decision principal usa valores criticos.\n');

else
    fprintf('  %-24s  (zatest no disponible en este entorno)\n', 'Zivot-Andrews');
end

%% 4. Identificacion ACF / PACF

figure('Name','ACF/PACF','Position',[100 100 800 350]);
subplot(1,2,1); autocorr(dlogIPC, 'NumLags', 24); title('ACF — \Delta log(IPC)');
subplot(1,2,2); parcorr(dlogIPC, 'NumLags', 24); title('PACF — \Delta log(IPC)');

%% 5. Grilla SARIMA/ARIMA + scoring compuesto

fprintf('\n===== ESTIMACION (grilla p=0:%d, q=0:%d, d=0:1) =====\n', MAX_P, MAX_Q);

Modelo = strings(0,1); pV=[]; dV=[]; qV=[];
AICv=[]; BICv=[]; AICcv=[]; HQv=[]; RMSEv=[]; MAPEv=[]; LogLikv=[];
LB12v=[]; RaicesOk=logical([]); Complexity=[];

for p = 0:MAX_P
  for q = 0:MAX_Q
    for d = 0:1
      if p==0 && d==0 && q==0, continue; end

      if d == 0
        constVal = NaN;   % se estima una media (analogo a include.mean=TRUE)
      else
        constVal = 0;     % sin drift, analogo al default de forecast::Arima con d>0
      end

      Mdl = arima('ARLags',1:p,'D',d,'MALags',1:q,'Constant',constVal);
      try
        [EstMdl, ~, logL] = estimate(Mdl, logIPC, 'Display','off');
      catch
        continue
      end

      R = summarize(EstMdl);
      k_full = R.NumEstimatedParameters;      % incluye Variance
      k_coef = max(k_full - 1, 0);            % excluye Variance (~ R coef())
      n_eff  = n - d;                         % observaciones efectivas

      res = infer(EstMdl, logIPC);
      rmse = sqrt(mean(res.^2));
      mape = mean(abs(res ./ logIPC)) * 100;

      aicc = R.AIC + (2*k_full*(k_full+1)) / max(n_eff - k_full - 1, 1);
      hq   = -2*logL + 2*k_full*log(log(n_eff));

      dof_lb = max(12 - k_coef, 1);
      try
        [~, p_lb] = lbqtest(res, 'Lags', 12, 'DoF', dof_lb);
      catch
        p_lb = NaN;
      end

      raices_ok = true;
      if p > 0
        arc = cell2mat(EstMdl.AR);
        rt  = roots([-fliplr(arc) 1]);
        raices_ok = raices_ok && all(abs(1./rt) < 1);
      end
      if q > 0
        mac = cell2mat(EstMdl.MA);
        rt  = roots([fliplr(mac) 1]);
        raices_ok = raices_ok && all(abs(1./rt) < 1);
      end

      Modelo(end+1,1)   = sprintf("ARIMA(%d,%d,%d)",p,d,q); %#ok<*SAGROW>
      pV(end+1,1)=p; dV(end+1,1)=d; qV(end+1,1)=q;
      AICv(end+1,1)=R.AIC; BICv(end+1,1)=R.BIC; AICcv(end+1,1)=aicc; HQv(end+1,1)=hq;
      RMSEv(end+1,1)=rmse; MAPEv(end+1,1)=mape; LogLikv(end+1,1)=logL;
      LB12v(end+1,1)=p_lb; RaicesOk(end+1,1)=raices_ok; Complexity(end+1,1)=p+q;
    end
  end
end

results = table(Modelo,pV,dV,qV,AICv,BICv,AICcv,HQv,RMSEv,MAPEv,LogLikv,LB12v,RaicesOk,Complexity, ...
    'VariableNames',{'Modelo','p','d','q','AIC','BIC','AICc','HQ','RMSE','MAPE','LogLik','LB12','Raices_ok','Complexity'});

fprintf('Modelos estimados: %d\n', height(results));

%% 6. Scoring compuesto (identico a la metodologia R)

results.RMSE_std = zscore(results.RMSE);
results.MAPE_std = zscore(results.MAPE);
results.AICc_std = zscore(results.AICc);
results.BIC_std  = zscore(results.BIC);
results.HQ_std   = zscore(results.HQ);
results.Cx_std   = zscore(results.Complexity);
results.LB_pen   = double(isnan(results.LB12) | results.LB12 < 0.05);
results.R_pen    = double(~results.Raices_ok);

results.Score = results.RMSE_std + results.MAPE_std + ...
    (results.AICc_std + results.BIC_std + results.HQ_std)/3 + ...
    0.3*results.Cx_std + 2*results.LB_pen + 2*results.R_pen;

results = sortrows(results, 'Score');

fprintf('\n--- Top 10 modelos por Score compuesto ---\n');
disp(results(1:min(10,height(results)), {'Modelo','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));

mejor = results(1,:);
fprintf('\n>>> MODELO GANADOR: %s (Score=%.4f)\n', mejor.Modelo, mejor.Score);

%% 7. Diagnostico del modelo ganador

p_g = mejor.p; d_g = mejor.d; q_g = mejor.q;
constVal_g = 0; if d_g==0, constVal_g = NaN; end
Mdl_g = arima('ARLags',1:p_g,'D',d_g,'MALags',1:q_g,'Constant',constVal_g);
[EstMdl_g, ~, logL_g] = estimate(Mdl_g, logIPC, 'Display','off');
Rg = summarize(EstMdl_g);

fprintf('\n===== DIAGNOSTICO: %s =====\n', mejor.Modelo);
disp(Rg.Table);

res_g = infer(EstMdl_g, logIPC);

figure('Name','Diagnostico residuos','Position',[100 100 900 700]);
subplot(2,2,1); plot(fechas, res_g, '-o', 'Color',[0.17 0.24 0.31]); yline(0,'r--');
title('Residuos'); grid on;
subplot(2,2,2); histogram(res_g, 12, 'FaceColor',[0.20 0.60 0.86]); title('Histograma de residuos');
subplot(2,2,3); autocorr(res_g, 'NumLags', 20); title('ACF de residuos');
subplot(2,2,4); qqplot(res_g); title('Q-Q Plot');

fprintf('\n--- Ljung-Box ---\n');
nCoef_g = Rg.NumEstimatedParameters - 1;
for h = [6 12 18 24]
  dof_h = max(h - nCoef_g, 1);
  [~, p_h] = lbqtest(res_g, 'Lags', h, 'DoF', dof_h);
  fprintf('  h=%2d : p=%.4f -> %s\n', h, p_h, ternary(p_h>0.05,"Ruido blanco","AUTOCORRELACION"));
end

[~, p_jb] = jbtest(res_g);
fprintf('Jarque-Bera: p=%.4f -> %s\n', p_jb, ternary(p_jb>0.05,"Normal","No normal"));

%% 8. Backtesting (ultimos H_BACKTEST meses)

fprintf('\n##########################################################\n');
fprintf('#              BACKTESTING — ULTIMOS %d MESES             #\n', H_BACKTEST);
fprintf('##########################################################\n');

nt = n - H_BACKTEST;
train_y = logIPC(1:nt);
fechas_test = fechas(nt+1:end);
ipc_obs = IPC(nt+1:end);

Mdl_bt = arima('ARLags',1:p_g,'D',d_g,'MALags',1:q_g,'Constant',constVal_g);
EstMdl_bt = estimate(Mdl_bt, train_y, 'Display','off');
[Ybt, YMSEbt] = forecast(EstMdl_bt, H_BACKTEST, 'Y0', train_y);

ipc_fc_bt   = exp(Ybt);
ipc_lo95_bt = exp(Ybt - norminv(0.975)*sqrt(YMSEbt));
ipc_hi95_bt = exp(Ybt + norminv(0.975)*sqrt(YMSEbt));

err_bt  = ipc_fc_bt - ipc_obs;
rmse_bt = sqrt(mean(err_bt.^2));
mape_bt = mean(abs(err_bt ./ ipc_obs)) * 100;
dentro_bt = (ipc_obs >= ipc_lo95_bt) & (ipc_obs <= ipc_hi95_bt);

fprintf('\n--- Comparacion en niveles de IPC ---\n');
for i = 1:H_BACKTEST
  fprintf('  %-10s obs=%8.2f pron=%8.2f err=%+7.2f  %s\n', ...
      datestr(fechas_test(i),'mmm yyyy'), ipc_obs(i), ipc_fc_bt(i), err_bt(i), ...
      ternary(dentro_bt(i),"dentro IC95","FUERA IC95"));
end
fprintf('\nRMSE = %.4f | MAPE = %.2f%% | Dentro IC95: %d/%d\n', ...
    rmse_bt, mape_bt, sum(dentro_bt), H_BACKTEST);

figure('Name','Backtesting','Position',[100 100 800 450]);
plot(fechas(1:nt), IPC(1:nt), '-', 'Color',[0.74 0.76 0.78], 'LineWidth',1); hold on;
plot(fechas_test, ipc_obs, '-o', 'Color',[0.17 0.24 0.31], 'LineWidth',1.5);
plot(fechas_test, ipc_fc_bt, '--o', 'Color',[0.91 0.30 0.24], 'LineWidth',1.5);
fill([fechas_test; flipud(fechas_test)], [ipc_lo95_bt; flipud(ipc_hi95_bt)], ...
    [0.91 0.30 0.24], 'FaceAlpha',0.15, 'EdgeColor','none');
legend({'Entrenamiento','Observado','Pronostico','IC 95%'}, 'Location','northwest');
title(sprintf('Backtesting IPC — %s (RMSE=%.2f, MAPE=%.2f%%)', mejor.Modelo, rmse_bt, mape_bt));
grid on; hold off;

%% 9. Pronostico final

fprintf('\n##########################################################\n');
fprintf('#                 PRONOSTICO FINAL                       #\n');
fprintf('##########################################################\n');

[Yf, YMSEf] = forecast(EstMdl_g, H_FORECAST, 'Y0', logIPC);

ipc_proy     = exp(Yf);
ipc_lo80     = exp(Yf - norminv(0.90)*sqrt(YMSEf));
ipc_hi80     = exp(Yf + norminv(0.90)*sqrt(YMSEf));
ipc_lo95     = exp(Yf - norminv(0.975)*sqrt(YMSEf));
ipc_hi95     = exp(Yf + norminv(0.975)*sqrt(YMSEf));

fecha_ultima = fechas(end);
fechas_pron  = dateshift(fecha_ultima + calmonths(1:H_FORECAST)', 'end', 'month');

% Inflacion mensual (encadenada)
ipc_serie_completa = [IPC; ipc_proy];
inf_mensual_proy = (ipc_proy ./ [IPC(end); ipc_proy(1:end-1)] - 1) * 100;

% Inflacion interanual (base 12 meses atras, dentro de la propia muestra)
inf_ia_proy = nan(H_FORECAST,1);
for i = 1:H_FORECAST
  fecha_base = fechas_pron(i) - calmonths(12);
  [~, idx12] = min(abs(fechas - fecha_base));
  inf_ia_proy(i) = (ipc_proy(i) / IPC(idx12) - 1) * 100;
end

meses_es = {'Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio', ...
            'Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

fprintf('\n--- Proyeccion IPC e inflacion ---\n');
for i = 1:H_FORECAST
  fprintf('  %-15s IPC=%8.2f  Inf.mens=%7.4f%%  Inf.ia=%7.2f%%\n', ...
      sprintf('%s %d', meses_es{month(fechas_pron(i))}, year(fechas_pron(i))), ...
      ipc_proy(i), inf_mensual_proy(i), inf_ia_proy(i));
end

figure('Name','Fan Chart','Position',[100 100 800 450]);
plot(fechas, IPC, '-', 'Color',[0.17 0.24 0.31], 'LineWidth',1.3); hold on;
fechas_fan = [fecha_ultima; fechas_pron];
ipc_fan    = [IPC(end); ipc_proy];
plot(fechas_fan, ipc_fan, '-', 'Color',[0 0 1], 'LineWidth',1.5);
fill([fechas_pron; flipud(fechas_pron)], [ipc_lo95; flipud(ipc_hi95)], ...
    [0 0 1], 'FaceAlpha',0.10, 'EdgeColor','none');
fill([fechas_pron; flipud(fechas_pron)], [ipc_lo80; flipud(ipc_hi80)], ...
    [0 0 1], 'FaceAlpha',0.20, 'EdgeColor','none');
legend({'Historico','Proyectado','IC 95%','IC 80%'}, 'Location','northwest');
title(sprintf('Fan Chart — IPC proyectado (%s)', mejor.Modelo));
grid on; hold off;

%% 10. Exportar resultados

if ~exist(fullfile(scriptDir,'..','Outputs'), 'dir')
  mkdir(fullfile(scriptDir,'..','Outputs'));
end

T_proy = table(year(fechas_pron), string(meses_es(month(fechas_pron)))', ...
    round(ipc_proy,4), round(inf_mensual_proy,4), round(inf_ia_proy,2), ...
    round(ipc_lo95,4), round(ipc_hi95,4), ...
    'VariableNames', {'Anio','Mes','IPC_proyectado','Inf_mensual_pct','Inf_interanual_pct','IC95_inf','IC95_sup'});

T_bt = table(string(datestr(fechas_test,'mmm yyyy')), round(ipc_obs,2), ...
    round(ipc_fc_bt,2), round(err_bt,2), ...
    'VariableNames', {'Fecha','IPC_observado','IPC_pronostico','Error'});

T_scoring = results(:, {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});

T_metricas = table( ...
    ["Modelo ganador";"Score";"AICc";"BIC";"HQ";"RMSE (train)";"MAPE (train)"; ...
     "RMSE backtest (IPC)";"MAPE backtest (IPC)";"Dentro IC95"], ...
    [string(mejor.Modelo); string(round(mejor.Score,4)); string(round(mejor.AICc,2)); ...
     string(round(mejor.BIC,2)); string(round(mejor.HQ,2)); string(round(mejor.RMSE,6)); ...
     string(round(mejor.MAPE,4)); string(round(rmse_bt,4)); ...
     string(round(mape_bt,2))+"%"; string(sum(dentro_bt))+"/"+string(H_BACKTEST)], ...
    'VariableNames', {'Metrica','Valor'});

writetable(T_proy,     RUTA_OUTPUT, 'Sheet','Proyeccion_MATLAB');
writetable(T_bt,       RUTA_OUTPUT, 'Sheet','Backtesting_MATLAB');
writetable(T_scoring,  RUTA_OUTPUT, 'Sheet','Scoring_MATLAB');
writetable(T_metricas, RUTA_OUTPUT, 'Sheet','Metricas_MATLAB');

fprintf('\nResultados exportados a: %s\n', RUTA_OUTPUT);

%% 11. Resumen final

fprintf('\n==========================================================\n');
fprintf('                    RESUMEN FINAL\n');
fprintf('==========================================================\n');
fprintf('  Modelo ganador  : %s\n', mejor.Modelo);
fprintf('  Score compuesto : %.4f\n', mejor.Score);
fprintf('  AICc / BIC / HQ : %.2f / %.2f / %.2f\n', mejor.AICc, mejor.BIC, mejor.HQ);
fprintf('  Backtest RMSE   : %.4f (IPC)\n', rmse_bt);
fprintf('  Backtest MAPE   : %.2f%%\n', mape_bt);
fprintf('  IPC dic 2026    : %.2f\n', ipc_proy(end));
fprintf('  Inf. ia dic 2026: %.2f%%\n', inf_ia_proy(end));
fprintf('==========================================================\n');

%% Funciones auxiliares

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end