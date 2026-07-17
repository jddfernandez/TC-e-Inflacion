%% SARIMA_TC.m — Modelos Estacionales para el Tipo de Cambio Digital
%
% Replica en MATLAB de SARIMA_TC.Rmd — espejo de SARIMA_M2.m con los roles
% invertidos: el TC digital es ahora la variable a explicar, y la inflacion
% interanual (en vez del TC digital) es el regresor exogeno del Modelo C.
%   A: SARIMA puro
%   B: SARIMA + dummies Bloqueo/shocks puntuales (COVID y dic-2020 se
%      excluyen: preceden al inicio real del TC digital, oct-2023)
%   C: SARIMAX (+ inflacion interanual)
%
% Muestra corta (n~32, nov-2023 a jun-2026): grilla mas acotada y un solo
% horizonte de backtest (6 meses) frente a SARIMA_M2.m.
%
% Requiere: matlab workspace ARIMA.mat (BDBruta_IPC, BDBruta_TC).

clc; clear; close all;
PROJ_ROOT   = 'C:\Users\Juande\Documents\Scripts Python\TC e π';
scriptDir   = fullfile(PROJ_ROOT, 'matlab scripts');
cd(scriptDir);

%% 0. Configuracion
FECHA_INICIO = datetime(2023,11,1);
FECHA_FIN    = datetime(2026,6,1);
H_BACKTEST  = 6;
H_FORECAST  = 7;
MAX_P = 2; MAX_Q = 2; MAX_PS = 1; MAX_QS = 1;
RUTA_OUTPUT = fullfile(PROJ_ROOT, 'Outputs', 'Proyecciones_MATLAB.xlsx');

%% 1. Carga de datos

load('matlab workspace ARIMA.mat');   % BDBruta_IPC, BDBruta_TC

tc_m = retime(BDBruta_TC, 'monthly', 'mean');
tc_m = tc_m(tc_m.Time >= FECHA_INICIO & tc_m.Time <= FECHA_FIN, :);
tc_m = sortrows(tc_m, 'Time');

fechas = tc_m.Time;
TC     = tc_m.TC_DIGITAL;
n      = length(TC);
logTC  = log(TC);

fprintf('Serie TC digital: %s a %s (n=%d)\n', ...
    datestr(fechas(1),'mmm yyyy'), datestr(fechas(end),'mmm yyyy'), n);
fprintf('  Nota: muestra corta (< 3 anios) frente a los modelos de inflacion.\n');

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
if any(isnan(inf_ia))
    warning('inf_ia tiene %d valores NaN (fechas sin match exacto en IPC).', sum(isnan(inf_ia)));
end

%% 2. Analisis preliminar

dlogTC = diff(logTC);

figure('Name','Preliminar SARIMA TC','Position',[50 50 1000 700]);
subplot(2,2,1); plot(fechas, logTC, '-', 'Color',[0.17 0.24 0.31]);
title('log(TC digital)'); grid on;
subplot(2,2,2); plot(fechas(2:end), dlogTC, '-', 'Color',[0.16 0.50 0.73]); yline(0,'r--');
title('\Delta log(TC digital)'); grid on;
subplot(2,2,3); autocorr(dlogTC, 'NumLags', min(24,n-2));
title('ACF \Delta log(TC) — rezagos estacionales');
subplot(2,2,4); parcorr(dlogTC, 'NumLags', min(24,n-2));
title('PACF \Delta log(TC)');

%% 3. Pruebas de raiz unitaria

fprintf('\n===== PRUEBAS DE RAIZ UNITARIA =====\n');
k_adf_log  = floor((n - 1)^(1/3));
k_adf_dlog = floor((length(dlogTC) - 1)^(1/3));
l_pp_log   = floor(4 * ((n - 1)/100)^0.25);
l_pp_dlog  = floor(4 * ((length(dlogTC) - 1)/100)^0.25);

[~, p_adf_log]   = adftest(logTC,  'Model','TS', 'Lags',k_adf_log);
[~, p_pp_log]    = pptest(logTC,   'Model','TS', 'Lags',l_pp_log);
[~, p_kpss_log]  = kpsstest(logTC, 'Trend', true);
[~, p_adf_dlog]  = adftest(dlogTC, 'Model','TS', 'Lags',k_adf_dlog);
[~, p_pp_dlog]   = pptest(dlogTC,  'Model','TS', 'Lags',l_pp_dlog);
[~, p_kpss_dlog] = kpsstest(dlogTC);

za_ok = false;
try
    [~, p_za_log,  ~, ~, bd_log]  = zatest(logTC,  'Model','LS');
    [~, p_za_dlog, ~, ~, bd_dlog] = zatest(dlogTC, 'Model','LS');
    za_ok = true;
catch ME
    warning('zatest: %s', ME.message);
end

fprintf('  ADF  log=%.4f Dlog=%.4f | PP  log=%.4f Dlog=%.4f | KPSS log=%.4f Dlog=%.4f\n', ...
    p_adf_log, p_adf_dlog, p_pp_log, p_pp_dlog, p_kpss_log, p_kpss_dlog);
if za_ok
    fprintf('  ZA   log~%.4f  Dlog~%.4f\n', p_za_log, p_za_dlog);
end

d_rec = 1;
if p_adf_log < 0.05 && p_pp_log < 0.05 && p_kpss_log >= 0.05
    d_rec = 0;
    fprintf('\n  => log(TC) aparentemente estacionaria — usando d=0\n');
else
    fprintf('\n  => log(TC) integrada I(1) — usando d=1\n');
end

%% 4. Estacionalidad y determinacion de D

acf_vals = autocorr(dlogTC, 'NumLags', min(24,n-2));
if numel(acf_vals) >= 13
    acf_lag12 = acf_vals(13);
else
    acf_lag12 = 0;
end
umbral = 2 / sqrt(n-1);
D_rec = double(acf_lag12 > umbral);
fprintf('\n  ACF(lag=12) de Dlog(TC) = %.4f (umbral=%.4f)\n', acf_lag12, umbral);
fprintf('  => D recomendado = %d\n', D_rec);
fprintf('  NOTA: con n=%d, D=1 deja solo %d obs efectivas para la grilla — resultado fragil.\n', n, n-12);

%% 5. Dummies de intervencion (solo las que caen dentro de la muestra)

dummy_bloqueo = double(fechas >= datetime(2025,5,1) & fechas <= datetime(2025,6,1));
dummy_jun2025 = double(fechas == datetime(2025,6,1));
dummy_jul2025 = double(fechas == datetime(2025,7,1));
dummy_jun2026 = double(fechas == datetime(2026,6,1));
X_dummies = [dummy_bloqueo, dummy_jun2025, dummy_jul2025, dummy_jun2026];
X_C = [X_dummies, inf_ia];

%% 6. Grilla de estimacion: Modelos A, B, C

fprintf('\n===== ESTIMACION (p=0:%d, q=0:%d, P=0:%d, Q=0:%d, d=%d, D=%d fijo) =====\n', ...
    MAX_P, MAX_Q, MAX_PS, MAX_QS, d_rec, D_rec);

fprintf('Modelo A (SARIMA puro)...\n');
[res_A, mdl_A] = estimar_sarima_tc(logTC, n, d_rec, D_rec, MAX_P, MAX_Q, MAX_PS, MAX_QS, [], 'A');

fprintf('Modelo B (SARIMA + dummies)...\n');
[res_B, mdl_B] = estimar_sarima_tc(logTC, n, d_rec, D_rec, MAX_P, MAX_Q, MAX_PS, MAX_QS, X_dummies, 'B');

fprintf('Modelo C (SARIMAX + inflacion)...\n');
[res_C, mdl_C] = estimar_sarima_tc(logTC, n, d_rec, D_rec, MAX_P, MAX_Q, MAX_PS, MAX_QS, X_C, 'C');

fprintf('\n--- Top 10 Modelo A ---\n');
disp(res_A(1:min(10,height(res_A)), {'Modelo','D','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));
fprintf('Ganador A: %s (Score=%.4f)\n', mdl_A.nombre, mdl_A.score);

fprintf('\n--- Top 10 Modelo B ---\n');
disp(res_B(1:min(10,height(res_B)), {'Modelo','D','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));
fprintf('Ganador B: %s (Score=%.4f)\n', mdl_B.nombre, mdl_B.score);

fprintf('\n--- Top 10 Modelo C ---\n');
disp(res_C(1:min(10,height(res_C)), {'Modelo','D','AICc','BIC','HQ','RMSE','MAPE','LB12','Score'}));
fprintf('Ganador C: %s (Score=%.4f)\n', mdl_C.nombre, mdl_C.score);

%% 7. Diagnostico de residuos

n_mod = 3;
figure('Name','Diagnosticos SARIMA TC','Position',[50 50 1000 700]);
mod_list  = {mdl_A, mdl_B, mdl_C};
col_list  = {[0.15 0.68 0.38], [0.91 0.50 0.13], [0.56 0.27 0.68]};
lbl_list  = {'A: SARIMA', 'B: +Dummies', 'C: +Inflacion'};
for mi = 1:n_mod
    mm = mod_list{mi}; cc = col_list{mi};
    res_plot = nan(n,1);
    res_plot(end-length(mm.res)+1:end) = mm.res;
    subplot(n_mod,3,(mi-1)*3+1);
    plot(fechas, res_plot, '-', 'Color',cc); yline(0,'r--'); grid on;
    xlim([fechas(1) fechas(end)]);
    title(sprintf('%s: Residuos', lbl_list{mi}));
    subplot(n_mod,3,(mi-1)*3+2); autocorr(mm.res,'NumLags',min(15,length(mm.res)-2));
    title(sprintf('%s: ACF', lbl_list{mi}));
    subplot(n_mod,3,(mi-1)*3+3); qqplot(mm.res);
    title(sprintf('%s: Q-Q', lbl_list{mi}));
end

for mi = 1:n_mod
    mm = mod_list{mi};
    nc = mm.n_coef;
    try
        [~,p_lb] = lbqtest(mm.res,'Lags',min(12,length(mm.res)-nc-1),'DoF',max(min(12,length(mm.res)-nc-1)-nc,1));
    catch
        p_lb = NaN;
    end
    [~,p_jb] = jbtest(mm.res);
    fprintf('%s: LB=%.4f [%s]  JB=%.4f [%s]\n', lbl_list{mi}, ...
        p_lb, ternary(p_lb>0.05,'OK','AUTO'), p_jb, ternary(p_jb>0.05,'Normal','No-normal'));
end

%% 8. Backtesting (un solo horizonte: H_BACKTEST meses)

fprintf('\n===== BACKTESTING — ULTIMOS %d MESES =====\n', H_BACKTEST);
nt = n - H_BACKTEST;
tc_obs_bt = TC(nt+1:end);
fechas_bt  = fechas(nt+1:end);

tc_bt = cell(n_mod,1);
X_list = {[], X_dummies, X_C};
for mi = 1:n_mod
    mm = mod_list{mi};
    tc_bt{mi} = backtest_sarima_tc(logTC, nt, n, mm, X_list{mi}, H_BACKTEST, lbl_list{mi});
end

fprintf('\n  Fecha       TC obs  ');
for mi=1:n_mod, fprintf('  %-12s', lbl_list{mi}); end
fprintf('\n');
for i=1:H_BACKTEST
    fprintf('  %-10s  %7.4f  ', datestr(fechas_bt(i),'mmm yyyy'), tc_obs_bt(i));
    for mi=1:n_mod, fprintf('  %12.4f', tc_bt{mi}(i)); end
    fprintf('\n');
end

fprintf('\n  %-16s  RMSE (TC)   MAPE (%%)\n', 'Modelo');
for mi=1:n_mod
    rmse_i = sqrt(mean((tc_bt{mi}-tc_obs_bt).^2, 'omitnan'));
    mape_i = mean(abs((tc_bt{mi}-tc_obs_bt)./tc_obs_bt), 'omitnan')*100;
    fprintf('  %-16s  %10.4f   %8.2f%%\n', lbl_list{mi}, rmse_i, mape_i);
end

figure('Name','Backtesting TC','Position',[50 50 900 500]);
plot(fechas(1:nt), TC(1:nt), '-', 'Color',[0.74 0.76 0.78], 'LineWidth',1); hold on;
plot(fechas_bt, tc_obs_bt, '-o', 'Color',[0.17 0.24 0.31], 'LineWidth',1.5);
for mi=1:n_mod
    plot(fechas_bt, tc_bt{mi}, '--o', 'Color',col_list{mi}, 'LineWidth',1.3);
end
legend([{'Entrenamiento','Observado'}, lbl_list], 'Location','northwest');
title('Backtesting TC digital — Modelos SARIMA'); grid on; hold off;

rmse_ab = [sqrt(mean((tc_bt{1}-tc_obs_bt).^2,'omitnan')), sqrt(mean((tc_bt{2}-tc_obs_bt).^2,'omitnan'))];
[~, idx_min] = min(rmse_ab);
etiquetas_ab = {'A (SARIMA puro)','B (SARIMA + intervenciones)'};
mejor_ab = etiquetas_ab{idx_min};
fprintf('\n>>> Modelo base recomendado (backtest): %s\n', mejor_ab);
fprintf('    C se reporta como escenario exploratorio / sensibilidad a la inflacion.\n');

%% 9. Pronostico final

fprintf('\n===== PRONOSTICO FINAL =====\n');
fecha_ultima = fechas(end);
fechas_pron  = dateshift(fecha_ultima + calmonths(1:H_FORECAST)', 'end','month');

X_fut_dummy = zeros(H_FORECAST, size(X_dummies,2));
X_fut_C = [X_fut_dummy, repmat(inf_ia(end), H_FORECAST, 1)];  % inflacion futura = ultimo valor observado

tc_pron = cell(n_mod,1);
meses_es = {'Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio', ...
            'Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

X_fut_list = {[], X_fut_dummy, X_fut_C};
for mi=1:n_mod
    mm = mod_list{mi};
    [pron_tc, pron_lo95, pron_hi95] = forecast_sarima_tc(TC, n, mm, X_fut_list{mi}, H_FORECAST, lbl_list{mi});
    tc_pron{mi} = pron_tc;
    if mi==1, tc_lo95_A = pron_lo95; tc_hi95_A = pron_hi95; end
end

fprintf('\n  %-15s', 'Mes');
for mi=1:n_mod, fprintf('  %-14s', lbl_list{mi}); end
fprintf('\n');
for i=1:H_FORECAST
    fprintf('  %-15s', sprintf('%s %d',meses_es{month(fechas_pron(i))},year(fechas_pron(i))));
    for mi=1:n_mod
        fprintf('  %12.4f  ', tc_pron{mi}(i));
    end
    fprintf('\n');
end

figure('Name','Pronostico SARIMA TC','Position',[50 50 1000 500]);
fecha_fan = [fecha_ultima; fechas_pron];
plot(fechas, TC, '-', 'Color',[0.17 0.24 0.31], 'LineWidth',1.3); hold on;
fill([fechas_pron; flipud(fechas_pron)],[tc_lo95_A; flipud(tc_hi95_A)], ...
    col_list{1}, 'FaceAlpha',0.15,'EdgeColor','none');
for mi=1:n_mod
    plot(fecha_fan, [TC(end);tc_pron{mi}], '-', 'Color',col_list{mi}, 'LineWidth',1.5);
end
legend([{'Historico','IC95 A'}, lbl_list], 'Location','northwest');
title('Pronostico TC digital — SARIMA / +Dummies / +Inflacion'); grid on; hold off;

%% 10. Exportar resultados

if ~exist(fullfile(scriptDir,'..','Outputs'), 'dir')
    mkdir(fullfile(scriptDir,'..','Outputs'));
end

mes_pron_str = string(meses_es(month(fechas_pron)))';
T_proy = table(year(fechas_pron), mes_pron_str, ...
    round(tc_pron{1},4), round(tc_pron{2},4), round(tc_pron{3},4), ...
    'VariableNames',{'Anio','Mes','A_TC','B_TC','C_TC'});
writetable(T_proy, RUTA_OUTPUT, 'Sheet','TC_SARIMA_Proyeccion');

T_sc = table({'A: SARIMA';'B: +Dummies';'C: +Inflacion'}, ...
    {mdl_A.nombre; mdl_B.nombre; mdl_C.nombre}, ...
    [mdl_A.score; mdl_B.score; mdl_C.score], ...
    'VariableNames',{'Modelo','Ganador','Score'});
writetable(T_sc, RUTA_OUTPUT, 'Sheet','TC_SARIMA_Scoring');

T_bt = table(string(datestr(fechas_bt,'mmm yyyy')), round(tc_obs_bt,4), ...
    round(tc_bt{1},4), round(tc_bt{2},4), round(tc_bt{3},4), ...
    'VariableNames',{'Fecha','TC_observado','A_pronostico','B_pronostico','C_pronostico'});
writetable(T_bt, RUTA_OUTPUT, 'Sheet','TC_SARIMA_Backtesting');

fprintf('\nResultados exportados: %s\n', RUTA_OUTPUT);

%% =========================================================================
%% FUNCIONES LOCALES
%% =========================================================================

function [results, mejor] = estimar_sarima_tc(logTC, n, d_val, D_val, ...
        max_p, max_q, max_P, max_Q, X, ~)
use_xreg = ~isempty(X);

Modelo=strings(0); pV=[]; dV=[]; qV=[]; PV=[]; DsV=[]; QV=[];
AICv=[]; BICv=[]; AICcv=[]; HQv=[]; RMSEv=[]; MAPEv=[];
LogLikv=[]; LB12v=[]; RaicesOk=logical([]); Cx=[];

    if D_val == 1
        y  = logTC(13:end) - logTC(1:end-12);
        ny = n - 12;
        X_y = []; if use_xreg, X_y = X(13:end,:); end
    else
        y = logTC; ny = n; X_y = X;
    end

    for p = 0:max_p
     for q = 0:max_q
      for P = 0:max_P
       for Q = 0:max_Q
        if p==0&&q==0&&P==0&&Q==0, continue; end
        try
            cst_val = ternary(d_val==0, NaN, 0);
            base_args = {'ARLags',1:p,'D',d_val,'MALags',1:q,'Constant',cst_val};
            seas_args = {};
            if P>0||Q>0
                seas_args = {'Seasonality',12};
                if P>0, seas_args=[seas_args,{'SARLags',1:P}]; end
                if Q>0, seas_args=[seas_args,{'SMALags',1:Q}]; end
            end

            if ~use_xreg
                Mdl = arima(base_args{:}, seas_args{:});
                [EstMdl,~,logL] = estimate(Mdl, y, 'Display','off');
                res = infer(EstMdl, y);
            else
                rbase = {'ARLags',1:p,'D',d_val,'MALags',1:q,'Intercept',0};
                Mdl = regARIMA(rbase{:}, seas_args{:});
                [EstMdl,~,logL] = estimate(Mdl, y, 'X',X_y, 'Display','off');
                res = infer(EstMdl, y, 'X',X_y);
            end

            R   = summarize(EstMdl);
            k   = R.NumEstimatedParameters;
            n_eff = ny - d_val;
            aicc  = R.AIC + (2*k*(k+1)) / max(n_eff-k-1, 1);
            hq    = -2*logL + 2*k*log(log(n_eff));
            rmse  = sqrt(mean(res.^2));
            mape  = mean(abs(res ./ y))*100;

            k_coef = max(k-1,0);
            lags_lb = max(min(12, n_eff-k_coef-2), k_coef+1);
            try, [~,p_lb] = lbqtest(res,'Lags',lags_lb,'DoF',max(lags_lb-k_coef,1));
            catch, p_lb=NaN; end

            rok = true;
            try
                if p>0, ar=cell2mat(EstMdl.AR); rt=roots([-fliplr(ar) 1]); rok=rok&&all(abs(1./rt)<1); end
                if q>0, ma=cell2mat(EstMdl.MA); rt=roots([fliplr(ma)  1]); rok=rok&&all(abs(1./rt)<1); end
                if P>0, sar=cell2mat(EstMdl.SAR); rt=roots([-fliplr(sar) 1]); rok=rok&&all(abs(1./rt)<1); end
                if Q>0, sma=cell2mat(EstMdl.SMA); rt=roots([fliplr(sma)  1]); rok=rok&&all(abs(1./rt)<1); end
            catch; end

            nm = sprintf('(%d,%d,%d)(%d,%d,%d)[12]',p,d_val,q,P,D_val,Q);
            Modelo(end+1,1)=nm; pV(end+1,1)=p; dV(end+1,1)=d_val;
            qV(end+1,1)=q; PV(end+1,1)=P; DsV(end+1,1)=D_val; QV(end+1,1)=Q;
            AICv(end+1,1)=R.AIC; BICv(end+1,1)=R.BIC; AICcv(end+1,1)=aicc;
            HQv(end+1,1)=hq; RMSEv(end+1,1)=rmse; MAPEv(end+1,1)=mape;
            LogLikv(end+1,1)=logL; LB12v(end+1,1)=p_lb;
            RaicesOk(end+1,1)=rok; Cx(end+1,1)=p+q+P+Q;
        catch; continue; end
       end
      end
     end
    end

if isempty(Modelo)
    results = table();
    mejor = struct('p',0,'d',d_val,'q',0,'P',0,'D',D_val,'Q',0,'nombre','(sin convergencia)','score',NaN);
    return
end

results = table(Modelo,pV,dV,qV,PV,DsV,QV,AICv,BICv,AICcv,HQv, ...
    RMSEv,MAPEv,LogLikv,LB12v,RaicesOk,Cx, ...
    'VariableNames',{'Modelo','p','d','q','P','D','Q','AIC','BIC','AICc','HQ', ...
                     'RMSE','MAPE','LogLik','LB12','Raices_ok','Complexity'});

results.RMSE_s = zscore(results.RMSE); results.MAPE_s = zscore(results.MAPE);
results.AICc_s = zscore(results.AICc); results.BIC_s  = zscore(results.BIC);
results.HQ_s   = zscore(results.HQ);   results.Cx_s   = zscore(results.Complexity);
results.LB_pen = double(isnan(results.LB12) | results.LB12 < 0.05);
results.R_pen  = double(~results.Raices_ok);
results.Score  = results.RMSE_s + results.MAPE_s + ...
    (results.AICc_s+results.BIC_s+results.HQ_s)/3 + ...
    0.3*results.Cx_s + 2*results.LB_pen + 2*results.R_pen;
results = sortrows(results,'Score');

bst = results(1,:);
mejor.p = bst.p; mejor.d = bst.d; mejor.q = bst.q;
mejor.P = bst.P; mejor.D = bst.D; mejor.Q = bst.Q;
mejor.nombre = bst.Modelo{1}; mejor.score = bst.Score;

D_w = mejor.D;
if D_w == 1
    y_w = logTC(13:end)-logTC(1:end-12);
    X_w = []; if use_xreg, X_w=X(13:end,:); end
else
    y_w = logTC; X_w = X;
end
cst_w = ternary(mejor.d==0, NaN, 0);
base_w = {'ARLags',1:mejor.p,'D',mejor.d,'MALags',1:mejor.q,'Constant',cst_w};
seas_w = {};
if mejor.P>0||mejor.Q>0
    seas_w={'Seasonality',12};
    if mejor.P>0, seas_w=[seas_w,{'SARLags',1:mejor.P}]; end
    if mejor.Q>0, seas_w=[seas_w,{'SMALags',1:mejor.Q}]; end
end
if ~use_xreg
    Mdl_w = arima(base_w{:}, seas_w{:});
    [mejor.EstMdl,~] = estimate(Mdl_w, y_w, 'Display','off');
    mejor.res = infer(mejor.EstMdl, y_w);
else
    rb_w = {'ARLags',1:mejor.p,'D',mejor.d,'MALags',1:mejor.q,'Intercept',0};
    Mdl_w = regARIMA(rb_w{:}, seas_w{:});
    [mejor.EstMdl,~] = estimate(Mdl_w, y_w, 'X',X_w, 'Display','off');
    mejor.res = infer(mejor.EstMdl, y_w, 'X',X_w);
end
mejor.n_coef = max(summarize(mejor.EstMdl).NumEstimatedParameters - 1, 0);
mejor.y_fit  = y_w;
mejor.use_xreg = use_xreg;
end

% -------------------------------------------------------------------------
function tc_bt = backtest_sarima_tc(logTC, nt, n, mm, X_mi, H_BT, etiqueta)
D_w   = mm.D; d_w = mm.d;
cst_w = ternary(d_w==0, NaN, 0);
base_w = {'ARLags',1:mm.p,'D',d_w,'MALags',1:mm.q,'Constant',cst_w};
seas_w = {};
if mm.P>0||mm.Q>0
    seas_w={'Seasonality',12};
    if mm.P>0, seas_w=[seas_w,{'SARLags',1:mm.P}]; end
    if mm.Q>0, seas_w=[seas_w,{'SMALags',1:mm.Q}]; end
end
if D_w==1
    y_tr = logTC(13:nt)-logTC(1:nt-12);
    X_tr = []; X_te = [];
    if mm.use_xreg, X_tr=X_mi(13:nt,:); X_te=X_mi(nt+1:n,:); end
else
    y_tr = logTC(1:nt);
    if mm.use_xreg
        X_tr=X_mi(1:nt,:); X_te=X_mi(nt+1:n,:);
    else
        X_tr=[]; X_te=[];
    end
end
try
    if ~mm.use_xreg
        Mdl_bt = arima(base_w{:}, seas_w{:});
        EstMdl_bt = estimate(Mdl_bt, y_tr, 'Display','off');
        Ybt = forecast(EstMdl_bt, H_BT, 'Y0', y_tr);
    else
        rb = {'ARLags',1:mm.p,'D',d_w,'MALags',1:mm.q,'Intercept',0};
        Mdl_bt = regARIMA(rb{:}, seas_w{:});
        EstMdl_bt = estimate(Mdl_bt, y_tr, 'X',X_tr, 'Display','off');
        Ybt = forecast(EstMdl_bt, H_BT, 'Y0', y_tr, 'XF', X_te);
    end
    if D_w==1
        logTC_bt = nan(H_BT,1);
        for i=1:H_BT, logTC_bt(i)=logTC(nt-12+i)+Ybt(i); end
    else
        logTC_bt = Ybt;
    end
    tc_bt = exp(logTC_bt);
catch ME
    warning('backtest_sarima_tc:fallo', 'Backtest de %s fallo (%s): %s', etiqueta, ME.identifier, ME.message);
    tc_bt = nan(H_BT,1);
end
end

% -------------------------------------------------------------------------
function [tc_p, lo95, hi95] = forecast_sarima_tc(TC, n, mm, X_fut, H, etiqueta)
D_w=mm.D; d_w=mm.d;
try
    if ~mm.use_xreg
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit);
    else
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit, 'XF', X_fut);
    end
    if D_w==1
        logTC_p=nan(H,1);
        for i=1:H, logTC_p(i)=log(TC(n-12+i))+Yf(i); end
    else
        logTC_p=Yf;
    end
    tc_p = exp(logTC_p);
    lo95 = exp(logTC_p - norminv(0.975)*sqrt(YMSEf));
    hi95 = exp(logTC_p + norminv(0.975)*sqrt(YMSEf));
catch ME
    warning('forecast_sarima_tc:fallo', 'Pronostico de %s fallo (%s): %s', etiqueta, ME.identifier, ME.message);
    tc_p=nan(H,1); lo95=nan(H,1); hi95=nan(H,1);
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out=a; else, out=b; end
end
