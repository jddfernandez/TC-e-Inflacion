%% ARIMAX_M2.m — ARIMA con Tipo de Cambio Digital como Regresor Exogeno
%
% Replica en MATLAB de ARIMAX_Inflacion.Rmd: extiende el modelo ganador M2
% (ARIMA_M2.m, log(IPC) sin regresor) agregando el TC digital como exogena
% via regARIMA (regresion con errores ARIMA). Se estiman dos grillas en
% paralelo — Base (arima, sin X) y ARIMAX (regARIMA, X=TC digital) — con el
% mismo scoring compuesto que el resto del proyecto, y la eleccion final
% entre ambas se decide por backtesting (fuera de muestra), no por score.
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
load('matlab workspace ARIMA.mat');   % BDBruta_IPC, BDBruta_TC

ipc_tbl = sortrows(BDBruta_IPC(BDBruta_IPC.A_o >= ANIO_INICIO, :), 'Time');
fechas  = ipc_tbl.Time;
IPC     = ipc_tbl.Valor;
logIPC  = log(IPC);
n       = length(logIPC);

fprintf('Serie ARIMAX cargada: %s a %s (n=%d)\n', ...
    datestr(fechas(1),'mmm yyyy'), datestr(fechas(end),'mmm yyyy'), n);

% --- TC digital mensual (media de TC diario por mes) ---
% Antes de oct-2023 el TC digital no existia (mercado paralelo/cripto);
% se interpola/retropropaga igual que en ARIMAX_Inflacion.Rmd y en
% SARIMA_M2.m, con piso de 6.96 si aun quedara NaN.
tc_ok = false;
try
    tc_m = retime(BDBruta_TC, 'monthly', 'mean');
    tc_m = tc_m(tc_m.Time >= fechas(1) & tc_m.Time <= fechas(end), :);
    tc_dig_v = nan(n,1);
    for i = 1:n
        tgt = tc_m.Time(year(tc_m.Time)==year(fechas(i)) & month(tc_m.Time)==month(fechas(i)));
        if ~isempty(tgt)
            idx = find(tc_m.Time == tgt, 1);
            tc_dig_v(i) = tc_m.TC_DIGITAL(idx);
        end
    end
    tc_dig_v = fillmissing(tc_dig_v, 'linear');
    tc_dig_v = fillmissing(tc_dig_v, 'previous');
    tc_dig_v(isnan(tc_dig_v)) = 6.96;
    tc_ok = true;
catch ME
    warning('No se pudo procesar TC digital: %s — ARIMAX omitido.', ME.message);
end
if ~tc_ok
    error('ARIMAX_M2:sinTC', 'No se pudo construir el regresor TC digital; revisar hoja TC de BD Bruta.xlsx.');
end

fprintf('  TC digital: media=%.4f  min=%.4f  max=%.4f\n', mean(tc_dig_v), min(tc_dig_v), max(tc_dig_v));
fprintf('  Correlacion log(IPC) vs TC digital: r=%.4f\n', corr(logIPC, tc_dig_v));

%% 2. Analisis preliminar

dlogIPC = diff(logIPC);

figure('Name','Datos ARIMAX','Position',[50 50 1000 700]);
subplot(2,2,1); plot(fechas, IPC, '-o', 'Color',[0.17 0.24 0.31]); title('IPC (base 2017=100)'); grid on;
subplot(2,2,2); plot(fechas, tc_dig_v, '-o', 'Color',[0.09 0.63 0.52]); title('TC digital (mensualizado)'); grid on;
subplot(2,2,3); plot(fechas(2:end), dlogIPC, '-o', 'Color',[0.16 0.50 0.73]); yline(0,'r--');
title('\Delta log(IPC)'); grid on;
subplot(2,2,4); scatter(tc_dig_v, logIPC, 24, [0.56 0.27 0.68], 'filled');
title(sprintf('log(IPC) vs TC digital (r=%.3f)', corr(logIPC,tc_dig_v))); xlabel('TC digital'); ylabel('log(IPC)'); grid on;

%% 3. Grilla de estimacion: Base (sin X) vs ARIMAX (X = TC digital)

fprintf('\n===== ESTIMACION (grilla p=0:%d, q=0:%d, d=0:1) =====\n', MAX_P, MAX_Q);

fprintf('Modelo Base (ARIMA, sin regresor)...\n');
[res_Base, mdl_Base] = estimar_arimax(logIPC, n, MAX_P, MAX_Q, [], 'Base');
fprintf('  Ganador Base  : %s (Score=%.4f)\n', mdl_Base.nombre, mdl_Base.score);

fprintf('Modelo ARIMAX (+TC digital)...\n');
[res_ARIMAX, mdl_ARIMAX] = estimar_arimax(logIPC, n, MAX_P, MAX_Q, tc_dig_v, 'ARIMAX');
fprintf('  Ganador ARIMAX: %s (Score=%.4f)\n', mdl_ARIMAX.nombre, mdl_ARIMAX.score);

fprintf('\n--- Top 10 Base ---\n');
disp(res_Base(1:min(10,height(res_Base)), {'Modelo','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));
fprintf('\n--- Top 10 ARIMAX ---\n');
disp(res_ARIMAX(1:min(10,height(res_ARIMAX)), {'Modelo','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));

%% 4. Coeficiente del regresor TC digital

Tw = summarize(mdl_ARIMAX.EstMdl).Table;
idx_beta = find(startsWith(string(Tw.Properties.RowNames), 'Beta'));
fprintf('\n===== COEFICIENTE TC DIGITAL (ganador ARIMAX: %s) =====\n', mdl_ARIMAX.nombre);
if ~isempty(idx_beta)
    disp(Tw(idx_beta,:));
    tc_sig = Tw.PValue(idx_beta(1)) < 0.05;
    fprintf('  %s\n', ternary(tc_sig, 'Significativo al 5%.', 'NO significativo al 5% -> tratar como exploratorio hasta ver el backtest.'));
else
    fprintf('  No se pudo ubicar el coeficiente Beta en la tabla de resumen.\n');
end

%% 5. Diagnostico de residuos (ganadores Base y ARIMAX)

figure('Name','Diagnostico ARIMAX','Position',[50 50 1000 500]);
col_base = [0.16 0.50 0.73]; col_arimax = [0.56 0.27 0.68];
subplot(2,3,1); plot(fechas(end-length(mdl_Base.res)+1:end), mdl_Base.res, '-', 'Color',col_base); yline(0,'r--');
title('Base: Residuos'); grid on;
subplot(2,3,2); autocorr(mdl_Base.res,'NumLags',20); title('Base: ACF');
subplot(2,3,3); qqplot(mdl_Base.res); title('Base: Q-Q');
subplot(2,3,4); plot(fechas(end-length(mdl_ARIMAX.res)+1:end), mdl_ARIMAX.res, '-', 'Color',col_arimax); yline(0,'r--');
title('ARIMAX: Residuos'); grid on;
subplot(2,3,5); autocorr(mdl_ARIMAX.res,'NumLags',20); title('ARIMAX: ACF');
subplot(2,3,6); qqplot(mdl_ARIMAX.res); title('ARIMAX: Q-Q');

fprintf('\n===== LJUNG-BOX Y JARQUE-BERA (ganadores) =====\n');
for h = [6 12 18]
    [~,p_lb_b] = lbqtest(mdl_Base.res,  'Lags',h,'DoF',max(h-mdl_Base.n_coef,1));
    [~,p_lb_x] = lbqtest(mdl_ARIMAX.res,'Lags',h,'DoF',max(h-mdl_ARIMAX.n_coef,1));
    fprintf('  LB h=%2d : Base p=%.4f [%s]   ARIMAX p=%.4f [%s]\n', h, ...
        p_lb_b, ternary(p_lb_b>0.05,'OK','AUTO'), p_lb_x, ternary(p_lb_x>0.05,'OK','AUTO'));
end
[~,p_jb_b] = jbtest(mdl_Base.res); [~,p_jb_x] = jbtest(mdl_ARIMAX.res);
fprintf('  Jarque-Bera: Base p=%.4f [%s]   ARIMAX p=%.4f [%s]\n', ...
    p_jb_b, ternary(p_jb_b>0.05,'Normal','No-normal'), p_jb_x, ternary(p_jb_x>0.05,'Normal','No-normal'));

%% 6. Backtesting — criterio decisivo para elegir entre Base y ARIMAX

fprintf('\n##########################################################\n');
fprintf('#   BACKTESTING — ULTIMOS %d MESES (criterio decisivo)   #\n', H_BACKTEST);
fprintf('##########################################################\n');

nt = n - H_BACKTEST;
fechas_test = fechas(nt+1:end);
ipc_obs = IPC(nt+1:end);

ipc_bt_base   = backtest_arimax(logIPC, nt, n, mdl_Base,   [],       H_BACKTEST, 'Base');
ipc_bt_arimax = backtest_arimax(logIPC, nt, n, mdl_ARIMAX, tc_dig_v, H_BACKTEST, 'ARIMAX');

err_base   = ipc_bt_base   - ipc_obs;
err_arimax = ipc_bt_arimax - ipc_obs;
rmse_base   = sqrt(mean(err_base.^2));   mape_base   = mean(abs(err_base./ipc_obs))*100;
rmse_arimax = sqrt(mean(err_arimax.^2)); mape_arimax = mean(abs(err_arimax./ipc_obs))*100;

fprintf('\n  Fecha       IPC obs     Base pron   ARIMAX pron\n');
for i = 1:H_BACKTEST
    fprintf('  %-10s  %8.2f    %8.2f     %8.2f\n', ...
        datestr(fechas_test(i),'mmm yyyy'), ipc_obs(i), ipc_bt_base(i), ipc_bt_arimax(i));
end
fprintf('\n  RMSE (IPC): Base=%.4f | ARIMAX=%.4f\n', rmse_base, rmse_arimax);
fprintf('  MAPE (%%)  : Base=%.2f%% | ARIMAX=%.2f%%\n', mape_base, mape_arimax);

gana_arimax = rmse_arimax < rmse_base;
if gana_arimax
    mdl_final = mdl_ARIMAX; etiqueta_final = 'ARIMAX (+TC digital)';
    fprintf('\n>>> GANADOR (backtest): ARIMAX — el TC digital mejora el pronostico fuera de muestra.\n');
else
    mdl_final = mdl_Base; etiqueta_final = 'ARIMA base (sin exogena)';
    fprintf('\n>>> GANADOR (backtest): ARIMA base — el TC digital no mejora el pronostico fuera de muestra;\n');
    fprintf('    ARIMAX queda como escenario exploratorio / sensibilidad al tipo de cambio.\n');
end

figure('Name','Backtesting ARIMAX','Position',[50 50 900 500]);
plot(fechas(1:nt), IPC(1:nt), '-', 'Color',[0.74 0.76 0.78], 'LineWidth',1); hold on;
plot(fechas_test, ipc_obs,       '-o',  'Color',[0.17 0.24 0.31], 'LineWidth',1.5);
plot(fechas_test, ipc_bt_base,   '--o', 'Color',col_base,        'LineWidth',1.3);
plot(fechas_test, ipc_bt_arimax, '--o', 'Color',col_arimax,      'LineWidth',1.3);
legend({'Entrenamiento','Observado','Base','ARIMAX'}, 'Location','northwest');
title(sprintf('Backtesting IPC — Base (RMSE=%.2f) vs ARIMAX (RMSE=%.2f)', rmse_base, rmse_arimax));
grid on; hold off;

%% 7. Pronostico final (ambos, y el elegido con IC)

fprintf('\n##########################################################\n');
fprintf('#                 PRONOSTICO FINAL                       #\n');
fprintf('##########################################################\n');

fecha_ultima = fechas(end);
fechas_pron  = dateshift(fecha_ultima + calmonths(1:H_FORECAST)', 'end', 'month');
X_fut = repmat(tc_dig_v(end), H_FORECAST, 1);   % TC digital futuro = ultimo valor observado

[ipc_pron_base,   lo95_base,   hi95_base]   = forecast_arimax(logIPC, n, mdl_Base,   [],    H_FORECAST);
[ipc_pron_arimax, lo95_arimax, hi95_arimax] = forecast_arimax(logIPC, n, mdl_ARIMAX, X_fut, H_FORECAST);

if gana_arimax
    ipc_pron_final = ipc_pron_arimax; lo95_final = lo95_arimax; hi95_final = hi95_arimax;
else
    ipc_pron_final = ipc_pron_base;   lo95_final = lo95_base;   hi95_final = hi95_base;
end

inf_mens_base   = (ipc_pron_base   ./ [IPC(end); ipc_pron_base(1:end-1)]   - 1) * 100;
inf_mens_arimax = (ipc_pron_arimax ./ [IPC(end); ipc_pron_arimax(1:end-1)] - 1) * 100;

inf_ia_base = nan(H_FORECAST,1); inf_ia_arimax = nan(H_FORECAST,1);
for i = 1:H_FORECAST
    f_base = fechas_pron(i) - calmonths(12);
    [~,ix] = min(abs(fechas - f_base));
    inf_ia_base(i)   = (ipc_pron_base(i)   / IPC(ix) - 1) * 100;
    inf_ia_arimax(i) = (ipc_pron_arimax(i) / IPC(ix) - 1) * 100;
end

meses_es = {'Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio', ...
            'Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

fprintf('\n  %-15s  Base (mens%%/ia%%)        ARIMAX (mens%%/ia%%)\n', 'Mes');
for i = 1:H_FORECAST
    fprintf('  %-15s  %7.4f%% / %6.2f%%      %7.4f%% / %6.2f%%\n', ...
        sprintf('%s %d', meses_es{month(fechas_pron(i))}, year(fechas_pron(i))), ...
        inf_mens_base(i), inf_ia_base(i), inf_mens_arimax(i), inf_ia_arimax(i));
end

figure('Name','Pronostico ARIMAX','Position',[50 50 900 500]);
plot(fechas, IPC, '-', 'Color',[0.17 0.24 0.31], 'LineWidth',1.3); hold on;
fill([fechas_pron; flipud(fechas_pron)], [lo95_final; flipud(hi95_final)], ...
    col_arimax, 'FaceAlpha',0.12, 'EdgeColor','none');
plot([fecha_ultima; fechas_pron], [IPC(end); ipc_pron_base],   '-', 'Color',col_base,   'LineWidth',1.5);
plot([fecha_ultima; fechas_pron], [IPC(end); ipc_pron_arimax], '-', 'Color',col_arimax, 'LineWidth',1.5);
legend({'Historico', sprintf('IC95 %s',etiqueta_final), 'Base','ARIMAX'}, 'Location','northwest');
title(sprintf('Pronostico IPC — Elegido: %s', etiqueta_final));
grid on; hold off;

%% 8. Contrafactual vs. Post-bloqueos (usa el modelo elegido por backtest)

esc_pb = struct('fechas', [datetime(2026,7,1),datetime(2026,8,1),datetime(2026,9,1), ...
                            datetime(2026,10,1),datetime(2026,11,1),datetime(2026,12,1)], ...
                'ia',     [10.1178, 10.5425, 11.6078, 12.2512, 13.0547, 13.6194]);

if gana_arimax, ia_final_all = inf_ia_arimax; else, ia_final_all = inf_ia_base; end
idx_jul  = fechas_pron >= datetime(2026,7,1);
ia_final_fut_all = ia_final_all(idx_jul);
n_cf = min(length(ia_final_fut_all), length(esc_pb.ia));
ia_final_fut = ia_final_fut_all(1:n_cf);

fprintf('\n===== CONTRAFACTUAL (%s) vs. POST-BLOQUEOS =====\n', etiqueta_final);
fprintf('  Mes         Contrafactual   Post-bloqueos  Dif (pp)\n');
fprintf('  %s\n', repmat('-',1,55));
for i = 1:n_cf
    fprintf('  %-12s  %8.2f%%       %8.2f%%      %+.2f\n', ...
        meses_es{month(esc_pb.fechas(i))}, ia_final_fut(i), esc_pb.ia(i), esc_pb.ia(i) - ia_final_fut(i));
end
fprintf('  Brecha dic 2026: %.1f pp\n', esc_pb.ia(end) - ia_final_fut(end));

%% 9. Exportar resultados

if ~exist(fullfile(scriptDir,'..','Outputs'), 'dir')
    mkdir(fullfile(scriptDir,'..','Outputs'));
end

mes_pron_str = string(meses_es(month(fechas_pron)))';
T_proy = table(year(fechas_pron), mes_pron_str, ...
    round(inf_mens_base,4),   round(inf_ia_base,2), ...
    round(inf_mens_arimax,4), round(inf_ia_arimax,2), ...
    'VariableNames',{'Anio','Mes','Base_mens','Base_ia','ARIMAX_mens','ARIMAX_ia'});
writetable(T_proy, RUTA_OUTPUT, 'Sheet','ARIMAX_Proyeccion');

T_bt = table(string(datestr(fechas_test,'mmm yyyy')), round(ipc_obs,2), ...
    round(ipc_bt_base,2), round(ipc_bt_arimax,2), round(err_base,2), round(err_arimax,2), ...
    'VariableNames', {'Fecha','IPC_observado','Base_pronostico','ARIMAX_pronostico','Base_error','ARIMAX_error'});
writetable(T_bt, RUTA_OUTPUT, 'Sheet','ARIMAX_Backtesting');

T_sc_base   = res_Base(:,   {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});
T_sc_arimax = res_ARIMAX(:, {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});
writetable(T_sc_base,   RUTA_OUTPUT, 'Sheet','ARIMAX_Scoring_Base');
writetable(T_sc_arimax, RUTA_OUTPUT, 'Sheet','ARIMAX_Scoring_ARIMAX');

T_metricas = table( ...
    ["Ganador Base";"Ganador ARIMAX";"RMSE backtest Base";"RMSE backtest ARIMAX"; ...
     "MAPE backtest Base";"MAPE backtest ARIMAX";"Modelo elegido"], ...
    [string(mdl_Base.nombre); string(mdl_ARIMAX.nombre); string(round(rmse_base,4)); string(round(rmse_arimax,4)); ...
     string(round(mape_base,2))+"%"; string(round(mape_arimax,2))+"%"; string(etiqueta_final)], ...
    'VariableNames', {'Metrica','Valor'});
writetable(T_metricas, RUTA_OUTPUT, 'Sheet','ARIMAX_Metricas');

fprintf('\nResultados exportados a: %s\n', RUTA_OUTPUT);

%% 10. Resumen final

fprintf('\n==========================================================\n');
fprintf('                    RESUMEN FINAL\n');
fprintf('==========================================================\n');
fprintf('  Ganador Base    : %s (Score=%.4f)\n', mdl_Base.nombre, mdl_Base.score);
fprintf('  Ganador ARIMAX  : %s (Score=%.4f)\n', mdl_ARIMAX.nombre, mdl_ARIMAX.score);
fprintf('  Backtest RMSE   : Base=%.4f | ARIMAX=%.4f\n', rmse_base, rmse_arimax);
fprintf('  Modelo elegido  : %s\n', etiqueta_final);
fprintf('  IPC dic 2026    : %.2f\n', ipc_pron_final(end));
fprintf('==========================================================\n');

%% =========================================================================
%% FUNCIONES LOCALES
%% =========================================================================

function [results, mejor] = estimar_arimax(y, n, max_p, max_q, X, ~)
% Grid search ARIMA(p,d,q) sobre y, con X opcional (regresion con errores
% ARIMA via regARIMA) o [] para el ARIMA base (via arima). d se recorre
% 0:1 dentro de la misma grilla (igual criterio que ARIMA_M2.m).

use_xreg = ~isempty(X);

Modelo = strings(0,1); pV=[]; dV=[]; qV=[];
AICv=[]; BICv=[]; AICcv=[]; HQv=[]; RMSEv=[]; MAPEv=[]; LogLikv=[];
LB12v=[]; RaicesOk=logical([]); Complexity=[];

for p = 0:max_p
  for q = 0:max_q
    for d = 0:1
      if p==0 && d==0 && q==0, continue; end
      try
        cst = NaN; if d ~= 0, cst = 0; end
        if ~use_xreg
          Mdl = arima('ARLags',1:p,'D',d,'MALags',1:q,'Constant',cst);
          [EstMdl,~,logL] = estimate(Mdl, y, 'Display','off');
          res = infer(EstMdl, y);
        else
          Mdl = regARIMA('ARLags',1:p,'D',d,'MALags',1:q,'Intercept',cst);
          [EstMdl,~,logL] = estimate(Mdl, y, 'X',X, 'Display','off');
          res = infer(EstMdl, y, 'X',X);
        end
      catch
        continue
      end

      R = summarize(EstMdl);
      k_full = R.NumEstimatedParameters;
      k_coef = max(k_full - 1, 0);
      n_eff  = n - d;

      rmse = sqrt(mean(res.^2));
      mape = mean(abs(res ./ y)) * 100;
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
        arc = cell2mat(EstMdl.AR); rt = roots([-fliplr(arc) 1]); raices_ok = raices_ok && all(abs(1./rt) < 1);
      end
      if q > 0
        mac = cell2mat(EstMdl.MA); rt = roots([fliplr(mac) 1]); raices_ok = raices_ok && all(abs(1./rt) < 1);
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

results.RMSE_std = zscore(results.RMSE); results.MAPE_std = zscore(results.MAPE);
results.AICc_std = zscore(results.AICc); results.BIC_std  = zscore(results.BIC);
results.HQ_std   = zscore(results.HQ);   results.Cx_std   = zscore(results.Complexity);
results.LB_pen   = double(isnan(results.LB12) | results.LB12 < 0.05);
results.R_pen    = double(~results.Raices_ok);
results.Score = results.RMSE_std + results.MAPE_std + ...
    (results.AICc_std + results.BIC_std + results.HQ_std)/3 + ...
    0.3*results.Cx_std + 2*results.LB_pen + 2*results.R_pen;
results = sortrows(results, 'Score');

bst = results(1,:);
mejor.p = bst.p; mejor.d = bst.d; mejor.q = bst.q;
mejor.nombre = bst.Modelo{1}; mejor.score = bst.Score; mejor.use_xreg = use_xreg;

% Re-estimar ganador en muestra completa para guardar EstMdl y residuos
cst_w = NaN; if mejor.d ~= 0, cst_w = 0; end
if ~use_xreg
    Mdl_w = arima('ARLags',1:mejor.p,'D',mejor.d,'MALags',1:mejor.q,'Constant',cst_w);
    [mejor.EstMdl,~] = estimate(Mdl_w, y, 'Display','off');
    mejor.res = infer(mejor.EstMdl, y);
else
    Mdl_w = regARIMA('ARLags',1:mejor.p,'D',mejor.d,'MALags',1:mejor.q,'Intercept',cst_w);
    [mejor.EstMdl,~] = estimate(Mdl_w, y, 'X',X, 'Display','off');
    mejor.res = infer(mejor.EstMdl, y, 'X',X);
end
mejor.n_coef = max(summarize(mejor.EstMdl).NumEstimatedParameters - 1, 0);
mejor.y_fit  = y;
end

% -------------------------------------------------------------------------
function ipc_bt = backtest_arimax(y, nt, n, mm, X, H_BT, etiqueta)
% Re-estima en muestra de entrenamiento y genera pronostico de backtest.
cst = NaN; if mm.d ~= 0, cst = 0; end
y_tr = y(1:nt);
X_tr = []; X_te = [];
if mm.use_xreg, X_tr = X(1:nt,:); X_te = X(nt+1:n,:); end
try
    if ~mm.use_xreg
        Mdl_bt = arima('ARLags',1:mm.p,'D',mm.d,'MALags',1:mm.q,'Constant',cst);
        EstMdl_bt = estimate(Mdl_bt, y_tr, 'Display','off');
        Ybt = forecast(EstMdl_bt, H_BT, 'Y0', y_tr);
    else
        Mdl_bt = regARIMA('ARLags',1:mm.p,'D',mm.d,'MALags',1:mm.q,'Intercept',cst);
        EstMdl_bt = estimate(Mdl_bt, y_tr, 'X',X_tr, 'Display','off');
        Ybt = forecast(EstMdl_bt, H_BT, 'Y0', y_tr, 'XF', X_te);
    end
    ipc_bt = exp(Ybt);
catch ME
    warning('backtest_arimax:fallo', 'Backtest de %s fallo (%s): %s', etiqueta, ME.identifier, ME.message);
    ipc_bt = nan(H_BT,1);
end
end

% -------------------------------------------------------------------------
function [ipc_p, lo95, hi95] = forecast_arimax(y, n, mm, X_fut, H)
% Pronostico H-pasos desde el modelo ganador (ya estimado en muestra
% completa dentro de estimar_arimax) con IC 95%.
try
    if ~mm.use_xreg
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit);
    else
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit, 'XF', X_fut);
    end
    ipc_p = exp(Yf);
    lo95  = exp(Yf - norminv(0.975)*sqrt(YMSEf));
    hi95  = exp(Yf + norminv(0.975)*sqrt(YMSEf));
catch ME
    warning('forecast_arimax:fallo', 'Pronostico fallo (%s): %s', ME.identifier, ME.message);
    ipc_p = nan(H,1); lo95 = nan(H,1); hi95 = nan(H,1);
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
