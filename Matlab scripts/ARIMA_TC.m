%% ARIMA_TC.m — Modelo Logaritmico (M2) del Tipo de Cambio Digital en Bolivia
%
% Replica en MATLAB de ARIMA_TC.Rmd: ARIMA(p,d,q) sobre log(TC digital),
% mismo espiritu metodologico que ARIMA_M2.m (grilla sistematica + seleccion
% por criterios de informacion combinados con precision y diagnostico de
% residuos), aplicado ahora al TC digital en vez de la inflacion.
%
% Requiere: matlab workspace ARIMA.mat (variable BDBruta_TC) en la misma
% carpeta que este script. Ejecutar actualizar_workspace.m primero si
% BD Bruta.xlsx cambio.

clc; clear; close all;
PROJ_ROOT = 'C:\Users\Juande\Documents\Scripts Python\TC e π';
scriptDir = fullfile(PROJ_ROOT, 'matlab scripts');
cd(scriptDir);

%% 0. Configuracion
% El TC digital solo cotiza realmente desde oct-2023; oct-2023 (13/31 dias)
% y el mes en curso suelen venir incompletos, asi que se usa la ventana de
% meses calendario completos: nov-2023 a jun-2026 (mismo ultimo mes cerrado
% que los modelos de inflacion).
FECHA_INICIO = datetime(2023,11,1);
FECHA_FIN    = datetime(2026,6,1);
H_BACKTEST  = 6;
H_FORECAST  = 7;
MAX_P = 4;
MAX_Q = 4;
RUTA_OUTPUT = fullfile(PROJ_ROOT, 'Outputs', 'Proyecciones_MATLAB.xlsx');

%% 1. Carga y mensualizacion del TC digital

load('matlab workspace ARIMA.mat');   % BDBruta_TC

tc_m = retime(BDBruta_TC, 'monthly', 'mean');
tc_m = tc_m(tc_m.Time >= FECHA_INICIO & tc_m.Time <= FECHA_FIN, :);
tc_m = sortrows(tc_m, 'Time');

fechas = tc_m.Time;
TC     = tc_m.TC_DIGITAL;
n      = length(TC);
logTC  = log(TC);

fprintf('Serie TC digital: %s a %s (n=%d)\n', ...
    datestr(fechas(1),'mmm yyyy'), datestr(fechas(end),'mmm yyyy'), n);
fprintf('  Nota: n=%d es una muestra corta frente a los modelos de inflacion\n', n);
fprintf('  (el TC digital solo tiene cotizacion real desde oct-2023).\n');

%% 2. Analisis preliminar

varTC   = (TC(2:end)./TC(1:end-1) - 1) * 100;   % M1: variacion mensual (%)
dlogTC  = diff(logTC);                           % M2: Dlog(TC)

figure('Name','Analisis preliminar TC','Position',[100 100 900 700]);
subplot(2,2,1); plot(fechas, TC, '-o', 'Color',[0.09 0.63 0.52]);
title('TC digital (Bs/USD)'); grid on;
subplot(2,2,2); plot(fechas, logTC, '-o', 'Color',[0.56 0.27 0.68]);
title('log(TC digital)'); grid on;
subplot(2,2,3); plot(fechas(2:end), dlogTC, '-o', 'Color',[0.16 0.50 0.73]);
yline(0,'r--'); title('\Delta log(TC digital)'); grid on;
subplot(2,2,4); histogram(varTC, 10, 'FaceColor',[0.91 0.30 0.24]);
title('M1: Distribucion var. mensual (%)'); grid on;

%% 3. Pruebas de raiz unitaria

fprintf('\n===== PRUEBAS DE RAIZ UNITARIA =====\n');
[~, p_adf_var]  = adftest(varTC);
[~, p_pp_var]   = pptest(varTC);
[~, p_kpss_var] = kpsstest(varTC);
[~, p_adf_log]  = adftest(logTC);
[~, p_pp_log]   = pptest(logTC);
[~, p_kpss_log] = kpsstest(logTC);
[~, p_adf_dlog]  = adftest(dlogTC);
[~, p_pp_dlog]   = pptest(dlogTC);
[~, p_kpss_dlog] = kpsstest(dlogTC);

za_ok = false;
try
    [~, p_za_var,  stat_za_var]  = zatest(varTC,  'Model','LS');
    [~, p_za_log,  stat_za_log]  = zatest(logTC,  'Model','LS');
    [~, p_za_dlog, stat_za_dlog] = zatest(dlogTC, 'Model','LS');
    za_ok = true;
catch ME
    warning('zatest no pudo calcularse: %s', ME.message);
end

fprintf('  %-20s  M1: Var.mensual   M2: log(TC)      Dlog(TC)\n','Prueba');
fprintf('  ADF     p=%.4f  p=%.4f  p=%.4f\n', p_adf_var, p_adf_log, p_adf_dlog);
fprintf('  PP      p=%.4f  p=%.4f  p=%.4f\n', p_pp_var,  p_pp_log,  p_pp_dlog);
fprintf('  KPSS    p=%.4f  p=%.4f  p=%.4f\n', p_kpss_var, p_kpss_log, p_kpss_dlog);
if za_ok
    fprintf('  ZA      p~%.4f  p~%.4f  p~%.4f\n', p_za_var, p_za_log, p_za_dlog);
end

%% 4. Grilla ARIMA + scoring compuesto (identico a ARIMA_M2.m)

fprintf('\n===== ESTIMACION M1 (grilla p=0:%d, q=0:%d, d=0:1) =====\n', MAX_P, MAX_Q);
results_m1 = estimar_grilla_tc(varTC, n, MAX_P, MAX_Q);
fprintf('\n===== ESTIMACION M2 (grilla p=0:%d, q=0:%d, d=0:1) =====\n', MAX_P, MAX_Q);
results_m2 = estimar_grilla_tc(logTC, n, MAX_P, MAX_Q);

mejor_m2 = results_m2(1,:);
fprintf('\n--- Top 10 M2 ---\n');
disp(results_m2(1:min(10,height(results_m2)), {'Modelo','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));

% M1 se fuerza al mismo p,q que M2 pero con d=0 (mismo criterio que
% ARIMA_Inflacion.Rmd / ARIMA_M2.m: la variacion mensual ya es una primera
% diferencia del log(TC)).
p_g = mejor_m2.p; q_g = mejor_m2.q; d_g = mejor_m2.d;
Mdl_m2 = arima('ARLags',1:p_g,'D',d_g,'MALags',1:q_g,'Constant',ternary(d_g==0,NaN,0));
[EstMdl_m2, ~, logL_m2] = estimate(Mdl_m2, logTC, 'Display','off');
Mdl_m1 = arima('ARLags',1:p_g,'D',0,'MALags',1:q_g,'Constant',NaN);
[EstMdl_m1, ~, logL_m1] = estimate(Mdl_m1, varTC, 'Display','off');

orden_m2 = sprintf('ARIMA(%d,%d,%d)', p_g, d_g, q_g);
orden_m1 = sprintf('ARIMA(%d,0,%d)',  p_g, q_g);
fprintf('\n>>> M2 GANADOR: %s (Score=%.4f)\n', orden_m2, mejor_m2.Score);
fprintf('>>> M1 FORZADO: %s (mismo AR/MA que M2, d=0)\n', orden_m1);

%% 5. Diagnostico de ambos ganadores

res_m1 = infer(EstMdl_m1, varTC);
res_m2 = infer(EstMdl_m2, logTC);

figure('Name','Diagnostico TC','Position',[100 100 900 700]);
subplot(2,4,1); plot(fechas(2:end), res_m1, '-o', 'Color',[0.91 0.30 0.24]); yline(0,'r--');
title(sprintf('M1: %s Residuos',orden_m1)); grid on;
subplot(2,4,2); histogram(res_m1, 10, 'FaceColor',[0.91 0.30 0.24]); title('M1: Histograma');
subplot(2,4,3); autocorr(res_m1, 'NumLags', 15); title('M1: ACF');
subplot(2,4,4); qqplot(res_m1); title('M1: Q-Q');
subplot(2,4,5); plot(fechas, res_m2, '-o', 'Color',[0.16 0.50 0.73]); yline(0,'r--');
title(sprintf('M2: %s Residuos',orden_m2)); grid on;
subplot(2,4,6); histogram(res_m2, 10, 'FaceColor',[0.16 0.50 0.73]); title('M2: Histograma');
subplot(2,4,7); autocorr(res_m2, 'NumLags', 15); title('M2: ACF');
subplot(2,4,8); qqplot(res_m2); title('M2: Q-Q');

nc_m1 = length(EstMdl_m1.AR) + length(EstMdl_m1.MA);
nc_m2 = length(EstMdl_m2.AR) + length(EstMdl_m2.MA);
fprintf('\n--- Ljung-Box ---\n');
for h = [6 12]
  [~,p_h1] = lbqtest(res_m1,'Lags',h,'DoF',max(h-nc_m1,1));
  [~,p_h2] = lbqtest(res_m2,'Lags',h,'DoF',max(h-nc_m2,1));
  fprintf('  h=%2d : M1 p=%.4f [%s]   M2 p=%.4f [%s]\n', h, ...
      p_h1, ternary(p_h1>0.05,'OK','AUTO'), p_h2, ternary(p_h2>0.05,'OK','AUTO'));
end

%% 6. Backtesting (ultimos H_BACKTEST meses)

fprintf('\n##########################################################\n');
fprintf('#              BACKTESTING — ULTIMOS %d MESES             #\n', H_BACKTEST);
fprintf('##########################################################\n');

nt1 = length(varTC) - H_BACKTEST;
train_m1 = varTC(1:nt1);
Mdl_bt1 = arima('ARLags',1:p_g,'D',0,'MALags',1:q_g,'Constant',NaN);
EstMdl_bt1 = estimate(Mdl_bt1, train_m1, 'Display','off');
Ybt1 = forecast(EstMdl_bt1, H_BACKTEST, 'Y0', train_m1);

nt2 = n - H_BACKTEST;
train_m2 = logTC(1:nt2);
Mdl_bt2 = arima('ARLags',1:p_g,'D',d_g,'MALags',1:q_g,'Constant',ternary(d_g==0,NaN,0));
EstMdl_bt2 = estimate(Mdl_bt2, train_m2, 'Display','off');
[Ybt2, YMSEbt2] = forecast(EstMdl_bt2, H_BACKTEST, 'Y0', train_m2);

tc_obs = TC(nt2+1:end);
fechas_test = fechas(nt2+1:end);

% varTC(i) = pct.change de TC(i) a TC(i+1), asi que entrenar hasta
% varTC(nt1) usa informacion hasta TC(nt1+1) — ese es el ancla correcta
% para encadenar el primer pronostico fuera de muestra (no TC(nt1)).
tc_fc1 = zeros(H_BACKTEST,1);
tc_pre = TC(nt1+1);
tc_fc1(1) = tc_pre * (1 + Ybt1(1)/100);
for i = 2:H_BACKTEST, tc_fc1(i) = tc_fc1(i-1) * (1 + Ybt1(i)/100); end

tc_fc2 = exp(Ybt2);

err1 = tc_fc1 - tc_obs; err2 = tc_fc2 - tc_obs;
rmse1 = sqrt(mean(err1.^2)); rmse2 = sqrt(mean(err2.^2));
mape1 = mean(abs(err1./tc_obs))*100; mape2 = mean(abs(err2./tc_obs))*100;

fprintf('\n  Fecha       Obs      M1 pron    M2 pron\n');
for i = 1:H_BACKTEST
  fprintf('  %-10s  %7.4f   %7.4f    %7.4f\n', datestr(fechas_test(i),'mmm yyyy'), tc_obs(i), tc_fc1(i), tc_fc2(i));
end
fprintf('\nRMSE: M1=%.4f | M2=%.4f\nMAPE: M1=%.2f%% | M2=%.2f%%\n', rmse1, rmse2, mape1, mape2);

%% 7. Pronostico final

fprintf('\n##########################################################\n');
fprintf('#                 PRONOSTICO FINAL                       #\n');
fprintf('##########################################################\n');

[Yf1] = forecast(EstMdl_m1, H_FORECAST, 'Y0', varTC);
[Yf2, YMSEf2] = forecast(EstMdl_m2, H_FORECAST, 'Y0', logTC);

tc_proy2 = exp(Yf2);
tc_lo95 = exp(Yf2 - norminv(0.975)*sqrt(YMSEf2));
tc_hi95 = exp(Yf2 + norminv(0.975)*sqrt(YMSEf2));

tc_proy1 = zeros(H_FORECAST,1);
tc_proy1(1) = TC(end) * (1 + Yf1(1)/100);
for i = 2:H_FORECAST, tc_proy1(i) = tc_proy1(i-1) * (1 + Yf1(i)/100); end

fecha_ultima = fechas(end);
fechas_pron  = dateshift(fecha_ultima + calmonths(1:H_FORECAST)', 'end', 'month');
meses_es = {'Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio', ...
            'Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

fprintf('\n  Mes              M1 TC       M2 TC\n');
for i = 1:H_FORECAST
  fprintf('  %-15s  %7.4f     %7.4f\n', ...
      sprintf('%s %d', meses_es{month(fechas_pron(i))}, year(fechas_pron(i))), tc_proy1(i), tc_proy2(i));
end

figure('Name','Fan Chart TC','Position',[100 100 800 450]);
plot(fechas, TC, '-', 'Color',[0.17 0.24 0.31], 'LineWidth',1.3); hold on;
fill([fechas_pron; flipud(fechas_pron)], [tc_lo95; flipud(tc_hi95)], [0 0 1], 'FaceAlpha',0.12,'EdgeColor','none');
plot([fecha_ultima;fechas_pron], [TC(end);tc_proy2], '-',  'Color',[0 0 1],        'LineWidth',1.5);
plot([fecha_ultima;fechas_pron], [TC(end);tc_proy1], '--', 'Color',[0.91 0.3 0.24],'LineWidth',1.5);
legend({'Historico','IC95 M2','M2 Logaritmico','M1 Directo'}, 'Location','northwest');
title(sprintf('Fan Chart — TC digital proyectado (M1:%s | M2:%s)', orden_m1, orden_m2));
grid on; hold off;

%% 8. Exportar resultados

if ~exist(fullfile(scriptDir,'..','Outputs'), 'dir')
  mkdir(fullfile(scriptDir,'..','Outputs'));
end

T_proy = table(year(fechas_pron), string(meses_es(month(fechas_pron)))', ...
    round(tc_proy1,4), round(tc_proy2,4), round(tc_lo95,4), round(tc_hi95,4), ...
    'VariableNames', {'Anio','Mes','M1_TC','M2_TC','M2_IC95_inf','M2_IC95_sup'});

T_bt = table(string(datestr(fechas_test,'mmm yyyy')), round(tc_obs,4), ...
    round(tc_fc1,4), round(tc_fc2,4), round(err1,4), round(err2,4), ...
    'VariableNames', {'Fecha','TC_observado','M1_pronostico','M2_pronostico','M1_error','M2_error'});

T_sc_m1 = results_m1(:, {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});
T_sc_m2 = results_m2(:, {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});

writetable(T_proy,   RUTA_OUTPUT, 'Sheet','TC_ARIMA_Proyeccion');
writetable(T_bt,     RUTA_OUTPUT, 'Sheet','TC_ARIMA_Backtesting');
writetable(T_sc_m1,  RUTA_OUTPUT, 'Sheet','TC_ARIMA_Scoring_M1');
writetable(T_sc_m2,  RUTA_OUTPUT, 'Sheet','TC_ARIMA_Scoring_M2');

fprintf('\nResultados exportados a: %s\n', RUTA_OUTPUT);

%% 9. Resumen final

fprintf('\n==========================================================\n');
fprintf('                    RESUMEN FINAL\n');
fprintf('==========================================================\n');
fprintf('  %-22s %-20s %-20s\n', '', 'M1 DIRECTO', 'M2 LOGARITMICO');
fprintf('  %-22s %-20s %-20s\n', 'Modelo', orden_m1, orden_m2);
fprintf('  %-22s %-20.4f %-20.4f\n', 'BT RMSE (TC)', rmse1, rmse2);
fprintf('  %-22s %-20s %-20s\n', 'BT MAPE (TC)', sprintf('%.2f%%',mape1), sprintf('%.2f%%',mape2));
fprintf('  %-22s %-20.4f %-20.4f\n', 'TC dic 2026', tc_proy1(6), tc_proy2(6));
fprintf('==========================================================\n');

%% =========================================================================
%% FUNCIONES LOCALES
%% =========================================================================

function results = estimar_grilla_tc(y, n, max_p, max_q)
Modelo = strings(0,1); pV=[]; dV=[]; qV=[];
AICv=[]; BICv=[]; AICcv=[]; HQv=[]; RMSEv=[]; MAPEv=[]; LogLikv=[];
LB12v=[]; RaicesOk=logical([]); Complexity=[];

for p = 0:max_p
  for q = 0:max_q
    for d = 0:1
      if p==0 && d==0 && q==0, continue; end
      try
        cst = NaN; if d ~= 0, cst = 0; end
        Mdl = arima('ARLags',1:p,'D',d,'MALags',1:q,'Constant',cst);
        [EstMdl, ~, logL] = estimate(Mdl, y, 'Display','off');
      catch
        continue
      end

      R = summarize(EstMdl);
      k_full = R.NumEstimatedParameters;
      k_coef = max(k_full - 1, 0);
      n_eff  = n - d;

      res = infer(EstMdl, y);
      rmse = sqrt(mean(res.^2));
      mape = mean(abs(res ./ y)) * 100;

      aicc = R.AIC + (2*k_full*(k_full+1)) / max(n_eff - k_full - 1, 1);
      hq   = -2*logL + 2*k_full*log(log(n_eff));

      lags_lb = max(min(12, n_eff-k_coef-2), k_coef+1);
      dof_lb  = max(lags_lb - k_coef, 1);
      try
        [~, p_lb] = lbqtest(res, 'Lags', lags_lb, 'DoF', dof_lb);
      catch
        p_lb = NaN;
      end

      raices_ok = true;
      if p > 0
        arc = cell2mat(EstMdl.AR); rt = roots([-fliplr(arc) 1]); raices_ok = raices_ok && all(abs(1./rt) < 1);
      end
      if q > 0
        mac = cell2mat(EstMdl.MA); rt = roots([fliplr(mac) 1]); raices_ok = raices_ok && all(abs(1./rt) < 1);
      end

      Modelo(end+1,1) = sprintf("ARIMA(%d,%d,%d)",p,d,q); %#ok<*SAGROW>
      pV(end+1,1)=p; dV(end+1,1)=d; qV(end+1,1)=q;
      AICv(end+1,1)=R.AIC; BICv(end+1,1)=R.BIC; AICcv(end+1,1)=aicc; HQv(end+1,1)=hq;
      RMSEv(end+1,1)=rmse; MAPEv(end+1,1)=mape; LogLikv(end+1,1)=logL;
      LB12v(end+1,1)=p_lb; RaicesOk(end+1,1)=raices_ok; Complexity(end+1,1)=p+q;
    end
  end
end

results = table(Modelo,pV,dV,qV,AICv,BICv,AICcv,HQv,RMSEv,MAPEv,LogLikv,LB12v,RaicesOk,Complexity, ...
    'VariableNames',{'Modelo','p','d','q','AIC','BIC','AICc','HQ','RMSE','MAPE','LogLik','LB12','Raices_ok','Complexity'});

results.RMSE_std = zscore(results.RMSE); results.MAPE_std = zscore(results.MAPE);
results.AICc_std = zscore(results.AICc); results.BIC_std  = zscore(results.BIC);
results.HQ_std   = zscore(results.HQ);   results.Cx_std   = zscore(results.Complexity);
results.LB_pen   = double(isnan(results.LB12) | results.LB12 < 0.05);
results.R_pen    = double(~results.Raices_ok);
results.Score = results.RMSE_std + results.MAPE_std + ...
    (results.AICc_std + results.BIC_std + results.HQ_std)/3 + ...
    0.3*results.Cx_std + 2*results.LB_pen + 2*results.R_pen;
results = sortrows(results, 'Score');
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
