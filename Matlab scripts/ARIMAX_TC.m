%% ARIMAX_TC.m — TC Digital con Inflacion Interanual como Regresor Exogeno
%
% Replica en MATLAB de ARIMAX_TC.Rmd — espejo de ARIMAX_M2.m con los roles
% invertidos: el TC digital es la variable a explicar, y la inflacion
% interanual (via regARIMA) es el regresor exogeno. Se estiman dos grillas
% en paralelo — Base (arima, sin X) y ARIMAX (regARIMA, X=inflacion) — y la
% eleccion final se decide por backtesting (fuera de muestra), no por score.
%
% Requiere: matlab workspace ARIMA.mat (BDBruta_IPC, BDBruta_TC).

clc; clear; close all;
PROJ_ROOT = 'C:\Users\Juande\Documents\Scripts Python\TC e π';
scriptDir = fullfile(PROJ_ROOT, 'matlab scripts');
cd(scriptDir);

%% 0. Configuracion
FECHA_INICIO = datetime(2023,11,1);
FECHA_FIN    = datetime(2026,6,1);
H_BACKTEST  = 6;
H_FORECAST  = 7;
MAX_P = 4;
MAX_Q = 4;
RUTA_OUTPUT = fullfile(PROJ_ROOT, 'Outputs', 'Proyecciones_MATLAB.xlsx');
addpath 'C:\Users\Juande\Documents\Scripts Python\TC e π\matlab scripts';

%% 1. Carga y mensualizacion

load('matlab workspace ARIMA.mat');   % BDBruta_IPC, BDBruta_TC

tc_m = retime(BDBruta_TC, 'monthly', 'mean');
tc_m = tc_m(tc_m.Time >= FECHA_INICIO & tc_m.Time <= FECHA_FIN, :);
tc_m = sortrows(tc_m, 'Time');

fechas = tc_m.Time;
TC     = tc_m.TC_DIGITAL;
n      = length(TC);
logTC  = log(TC);

fprintf('Serie TC digital cargada: %s a %s (n=%d)\n', ...
    datestr(fechas(1),'mmm yyyy'), datestr(fechas(end),'mmm yyyy'), n);

% --- Inflacion interanual alineada a las fechas del TC digital ---
ipc_tbl = sortrows(BDBruta_IPC, 'Time');
inf_ia = nan(n,1);
for i = 1:n
    f_base = fechas(i) - calmonths(12);
    idx_now  = find(ipc_tbl.Time == fechas(i), 1);
    idx_base = find(ipc_tbl.Time == f_base, 1);
    if ~isempty(idx_now) && ~isempty(idx_base)
        inf_ia(i) = (ipc_tbl.Valor(idx_now) / ipc_tbl.Valor(idx_base) - 1) * 100;
    end
end

fprintf('  Correlacion log(TC digital) vs inflacion interanual: r=%.4f\n', corr(logTC, inf_ia));

%% 2. Analisis preliminar

dlogTC = diff(logTC);

figure('Name','Datos ARIMAX TC','Position',[50 50 1000 700]);
subplot(2,2,1); plot(fechas, TC, '-o', 'Color',[0.09 0.63 0.52]); title('TC digital (Bs/USD)'); grid on;
subplot(2,2,2); plot(fechas, inf_ia, '-o', 'Color',[0.75 0.22 0.17]); title('Inflacion interanual (%)'); grid on;
subplot(2,2,3); plot(fechas(2:end), dlogTC, '-o', 'Color',[0.16 0.50 0.73]); yline(0,'r--');
title('\Delta log(TC digital)'); grid on;
subplot(2,2,4); scatter(inf_ia, logTC, 24, [0.56 0.27 0.68], 'filled');
title(sprintf('log(TC) vs inflacion ia (r=%.3f)', corr(logTC,inf_ia))); xlabel('Inflacion ia (%)'); ylabel('log(TC)'); grid on;

%% 3. Grilla de estimacion: Base (sin X) vs ARIMAX (X = inflacion interanual)

fprintf('\n===== ESTIMACION (grilla p=0:%d, q=0:%d, d=0:1) =====\n', MAX_P, MAX_Q);

fprintf('Modelo Base (ARIMA, sin regresor)...\n');
[res_Base, mdl_Base] = estimar_arimax_tc(logTC, n, MAX_P, MAX_Q, [], 'Base');
fprintf('  Ganador Base  : %s (Score=%.4f)\n', mdl_Base.nombre, mdl_Base.score);

fprintf('Modelo ARIMAX (+inflacion interanual)...\n');
[res_ARIMAX, mdl_ARIMAX] = estimar_arimax_tc(logTC, n, MAX_P, MAX_Q, inf_ia, 'ARIMAX');
fprintf('  Ganador ARIMAX: %s (Score=%.4f)\n', mdl_ARIMAX.nombre, mdl_ARIMAX.score);

fprintf('\n--- Top 10 Base ---\n');
disp(res_Base(1:min(10,height(res_Base)), {'Modelo','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));
fprintf('\n--- Top 10 ARIMAX ---\n');
disp(res_ARIMAX(1:min(10,height(res_ARIMAX)), {'Modelo','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));

%% 4. Coeficiente del regresor de inflacion

Tw = summarize(mdl_ARIMAX.EstMdl).Table;
idx_beta = find(startsWith(string(Tw.Properties.RowNames), 'Beta'));
fprintf('\n===== COEFICIENTE INFLACION INTERANUAL (ganador ARIMAX: %s) =====\n', mdl_ARIMAX.nombre);
if ~isempty(idx_beta)
    disp(Tw(idx_beta,:));
    ia_sig = Tw.PValue(idx_beta(1)) < 0.05;
    fprintf('  %s\n', ternary(ia_sig, 'Significativo al 5%.', 'NO significativo al 5% -> tratar como exploratorio hasta ver el backtest.'));
else
    fprintf('  No se pudo ubicar el coeficiente Beta en la tabla de resumen.\n');
end

%% 5. Diagnostico de residuos

figure('Name','Diagnostico ARIMAX TC','Position',[50 50 1000 500]);
col_base = [0.16 0.50 0.73]; col_arimax = [0.56 0.27 0.68];
subplot(2,3,1); plot(fechas(end-length(mdl_Base.res)+1:end), mdl_Base.res, '-', 'Color',col_base); yline(0,'r--');
title('Base: Residuos'); grid on;
subplot(2,3,2); autocorr(mdl_Base.res,'NumLags',min(15,length(mdl_Base.res)-2)); title('Base: ACF');
subplot(2,3,3); qqplot(mdl_Base.res); title('Base: Q-Q');
subplot(2,3,4); plot(fechas(end-length(mdl_ARIMAX.res)+1:end), mdl_ARIMAX.res, '-', 'Color',col_arimax); yline(0,'r--');
title('ARIMAX: Residuos'); grid on;
subplot(2,3,5); autocorr(mdl_ARIMAX.res,'NumLags',min(15,length(mdl_ARIMAX.res)-2)); title('ARIMAX: ACF');
subplot(2,3,6); qqplot(mdl_ARIMAX.res); title('ARIMAX: Q-Q');

fprintf('\n===== LJUNG-BOX Y JARQUE-BERA (ganadores) =====\n');
for h = [6 12]
    try
      [~,p_lb_b] = lbqtest(mdl_Base.res,  'Lags',min(h,length(mdl_Base.res)-mdl_Base.n_coef-1),  'DoF',max(h-mdl_Base.n_coef,1));
    catch, p_lb_b = NaN; end
    try
      [~,p_lb_x] = lbqtest(mdl_ARIMAX.res,'Lags',min(h,length(mdl_ARIMAX.res)-mdl_ARIMAX.n_coef-1),'DoF',max(h-mdl_ARIMAX.n_coef,1));
    catch, p_lb_x = NaN; end
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
tc_obs = TC(nt+1:end);

tc_bt_base   = backtest_arimax_tc(logTC, nt, n, mdl_Base,   [],     H_BACKTEST, 'Base');
tc_bt_arimax = backtest_arimax_tc(logTC, nt, n, mdl_ARIMAX, inf_ia, H_BACKTEST, 'ARIMAX');

err_base   = tc_bt_base   - tc_obs;
err_arimax = tc_bt_arimax - tc_obs;
rmse_base   = sqrt(mean(err_base.^2));   mape_base   = mean(abs(err_base./tc_obs))*100;
rmse_arimax = sqrt(mean(err_arimax.^2)); mape_arimax = mean(abs(err_arimax./tc_obs))*100;

fprintf('\n  Fecha       TC obs   Base pron   ARIMAX pron\n');
for i = 1:H_BACKTEST
    fprintf('  %-10s  %7.4f    %7.4f     %7.4f\n', ...
        datestr(fechas_test(i),'mmm yyyy'), tc_obs(i), tc_bt_base(i), tc_bt_arimax(i));
end
fprintf('\n  RMSE (TC): Base=%.4f | ARIMAX=%.4f\n', rmse_base, rmse_arimax);
fprintf('  MAPE (%%) : Base=%.2f%% | ARIMAX=%.2f%%\n', mape_base, mape_arimax);

gana_arimax = rmse_arimax < rmse_base;
if gana_arimax
    etiqueta_final = 'ARIMAX (+inflacion)';
    fprintf('\n>>> GANADOR (backtest): ARIMAX — la inflacion mejora el pronostico fuera de muestra.\n');
else
    etiqueta_final = 'ARIMA base (sin exogena)';
    fprintf('\n>>> GANADOR (backtest): ARIMA base — la inflacion no mejora el pronostico fuera de muestra;\n');
    fprintf('    ARIMAX queda como escenario exploratorio / sensibilidad a la inflacion.\n');
end

figure('Name','Backtesting ARIMAX TC','Position',[50 50 900 500]);
plot(fechas(1:nt), TC(1:nt), '-', 'Color',[0.74 0.76 0.78], 'LineWidth',1); hold on;
plot(fechas_test, tc_obs,      '-o',  'Color',[0.17 0.24 0.31], 'LineWidth',1.5);
plot(fechas_test, tc_bt_base,  '--o', 'Color',col_base,         'LineWidth',1.3);
plot(fechas_test, tc_bt_arimax,'--o', 'Color',col_arimax,       'LineWidth',1.3);
legend({'Entrenamiento','Observado','Base','ARIMAX'}, 'Location','northwest');
title(sprintf('Backtesting TC digital — Base (RMSE=%.4f) vs ARIMAX (RMSE=%.4f)', rmse_base, rmse_arimax));
grid on; hold off;

%% 7. Pronostico final

fprintf('\n##########################################################\n');
fprintf('#                 PRONOSTICO FINAL                       #\n');
fprintf('##########################################################\n');

fecha_ultima = fechas(end);
fechas_pron  = dateshift(fecha_ultima + calmonths(1:H_FORECAST)', 'end', 'month');
X_fut = repmat(inf_ia(end), H_FORECAST, 1);   % inflacion futura = ultimo valor observado

[tc_pron_base,   lo95_base,   hi95_base]   = forecast_arimax_tc(logTC, n, mdl_Base,   [],    H_FORECAST);
[tc_pron_arimax, lo95_arimax, hi95_arimax] = forecast_arimax_tc(logTC, n, mdl_ARIMAX, X_fut, H_FORECAST);

meses_es = {'Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio', ...
            'Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

fprintf('\n  Mes              Base TC     ARIMAX TC\n');
for i = 1:H_FORECAST
  fprintf('  %-15s  %8.4f    %8.4f\n', ...
      sprintf('%s %d', meses_es{month(fechas_pron(i))}, year(fechas_pron(i))), tc_pron_base(i), tc_pron_arimax(i));
end

if gana_arimax
    tc_pron_final = tc_pron_arimax; lo95_final = lo95_arimax; hi95_final = hi95_arimax;
else
    tc_pron_final = tc_pron_base;   lo95_final = lo95_base;   hi95_final = hi95_base;
end

figure('Name','Pronostico ARIMAX TC','Position',[50 50 900 500]);
plot(fechas, TC, '-', 'Color',[0.17 0.24 0.31], 'LineWidth',1.3); hold on;
fill([fechas_pron; flipud(fechas_pron)], [lo95_final; flipud(hi95_final)], ...
    col_arimax, 'FaceAlpha',0.12, 'EdgeColor','none');
plot([fecha_ultima; fechas_pron], [TC(end); tc_pron_base],   '-', 'Color',col_base,   'LineWidth',1.5);
plot([fecha_ultima; fechas_pron], [TC(end); tc_pron_arimax], '-', 'Color',col_arimax, 'LineWidth',1.5);
legend({'Historico', sprintf('IC95 %s',etiqueta_final), 'Base','ARIMAX'}, 'Location','northwest');
title(sprintf('Pronostico TC digital — Elegido: %s', etiqueta_final));
grid on; hold off;

%% 8. Exportar resultados

if ~exist(fullfile(scriptDir,'..','Outputs'), 'dir')
    mkdir(fullfile(scriptDir,'..','Outputs'));
end

mes_pron_str = string(meses_es(month(fechas_pron)))';
T_proy = table(year(fechas_pron), mes_pron_str, ...
    round(tc_pron_base,4), round(tc_pron_arimax,4), ...
    'VariableNames',{'Anio','Mes','Base_TC','ARIMAX_TC'});
writetable(T_proy, RUTA_OUTPUT, 'Sheet','TC_ARIMAX_Proyeccion');

T_bt = table(string(datestr(fechas_test,'mmm yyyy')), round(tc_obs,4), ...
    round(tc_bt_base,4), round(tc_bt_arimax,4), round(err_base,4), round(err_arimax,4), ...
    'VariableNames', {'Fecha','TC_observado','Base_pronostico','ARIMAX_pronostico','Base_error','ARIMAX_error'});
writetable(T_bt, RUTA_OUTPUT, 'Sheet','TC_ARIMAX_Backtesting');

T_sc_base   = res_Base(:,   {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});
T_sc_arimax = res_ARIMAX(:, {'Modelo','p','d','q','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'});
writetable(T_sc_base,   RUTA_OUTPUT, 'Sheet','TC_ARIMAX_Scoring_Base');
writetable(T_sc_arimax, RUTA_OUTPUT, 'Sheet','TC_ARIMAX_Scoring_ARIMAX');

T_metricas = table( ...
    ["Ganador Base";"Ganador ARIMAX";"RMSE backtest Base";"RMSE backtest ARIMAX"; ...
     "MAPE backtest Base";"MAPE backtest ARIMAX";"Modelo elegido"], ...
    [string(mdl_Base.nombre); string(mdl_ARIMAX.nombre); string(round(rmse_base,4)); string(round(rmse_arimax,4)); ...
     string(round(mape_base,2))+"%"; string(round(mape_arimax,2))+"%"; string(etiqueta_final)], ...
    'VariableNames', {'Metrica','Valor'});
writetable(T_metricas, RUTA_OUTPUT, 'Sheet','TC_ARIMAX_Metricas');

fprintf('\nResultados exportados a: %s\n', RUTA_OUTPUT);

%% 9. Resumen final

fprintf('\n==========================================================\n');
fprintf('                    RESUMEN FINAL\n');
fprintf('==========================================================\n');
fprintf('  Ganador Base    : %s (Score=%.4f)\n', mdl_Base.nombre, mdl_Base.score);
fprintf('  Ganador ARIMAX  : %s (Score=%.4f)\n', mdl_ARIMAX.nombre, mdl_ARIMAX.score);
fprintf('  Backtest RMSE   : Base=%.4f | ARIMAX=%.4f\n', rmse_base, rmse_arimax);
fprintf('  Modelo elegido  : %s\n', etiqueta_final);
fprintf('  TC dic 2026     : %.4f\n', tc_pron_final(6));
fprintf('==========================================================\n');

%% =========================================================================
%% FUNCIONES LOCALES
%% =========================================================================

function [results, mejor] = estimar_arimax_tc(y, n, max_p, max_q, X, ~)
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

      lags_lb = max(min(12, n_eff-k_coef-2), k_coef+1);
      try
        [~, p_lb] = lbqtest(res, 'Lags', lags_lb, 'DoF', max(lags_lb-k_coef,1));
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
function tc_bt = backtest_arimax_tc(y, nt, n, mm, X, H_BT, etiqueta)
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
    tc_bt = exp(Ybt);
catch ME
    warning('backtest_arimax_tc:fallo', 'Backtest de %s fallo (%s): %s', etiqueta, ME.identifier, ME.message);
    tc_bt = nan(H_BT,1);
end
end

% -------------------------------------------------------------------------
function [tc_p, lo95, hi95] = forecast_arimax_tc(y, n, mm, X_fut, H)
try
    if ~mm.use_xreg
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit);
    else
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit, 'XF', X_fut);
    end
    tc_p = exp(Yf);
    lo95 = exp(Yf - norminv(0.975)*sqrt(YMSEf));
    hi95 = exp(Yf + norminv(0.975)*sqrt(YMSEf));
catch ME
    warning('forecast_arimax_tc:fallo', 'Pronostico fallo (%s): %s', ME.identifier, ME.message);
    tc_p = nan(H,1); lo95 = nan(H,1); hi95 = nan(H,1);
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
