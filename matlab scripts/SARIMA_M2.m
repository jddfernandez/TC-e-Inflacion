%% SARIMA_M2.m — Modelos Estacionales para Inflacion en Bolivia
%
% Replica en MATLAB de SARIMA_Inflacion.Rmd:
%   A: SARIMA puro
%   B: SARIMA + dummies COVID / Bloqueos
%   C: SARIMAX (+ tipo de cambio)
%
% Todos sobre log(IPC) desde 2019.  Mismo scoring compuesto que ARIMA_M2.m.
% Requiere workspace actualizado (ejecutar actualizar_workspace.m primero).

clc; clear; close all;
PROJ_ROOT   = 'C:\Users\Juande\Documents\Scripts Python\TC e π';
scriptDir   = fullfile(PROJ_ROOT, 'matlab scripts');
cd(scriptDir);

%% 0. Configuracion
ANIO_INICIO = 2019;
H_BACKTEST  = 6;
H_FORECAST  = 7;
MAX_P = 2; MAX_Q = 2; MAX_PS = 1; MAX_QS = 1;   % grilla SARIMA
RUTA_OUTPUT = fullfile(PROJ_ROOT, 'Outputs', 'Proyecciones_MATLAB.xlsx');

%% 1. Carga de datos

load('matlab workspace ARIMA.mat');   % BDBruta_IPC, BDBruta_TC

% --- IPC desde 2019 ---
ipc_tbl = sortrows(BDBruta_IPC(BDBruta_IPC.A_o >= ANIO_INICIO, :), 'Time');
fechas  = ipc_tbl.Time;
IPC     = ipc_tbl.Valor;
logIPC  = log(IPC);
n       = length(logIPC);

fprintf('Serie SARIMA: %s a %s (n=%d)\n', ...
    datestr(fechas(1),'mmm yyyy'), datestr(fechas(end),'mmm yyyy'), n);

% --- TC mensual (media de TC diario por mes) ---
tc_ok = false;
try
    % retime agrega por mes (usa primer dia de cada mes como indice)
    tc_m = retime(BDBruta_TC, 'monthly', 'mean');
    % Recortar a periodo de la muestra IPC
    tc_m = tc_m(tc_m.Time >= fechas(1) & tc_m.Time <= fechas(end), :);
    % Asegurar misma longitud que IPC (join por mes)
    tc_dig_v = nan(n,1);  tc_ref_v = nan(n,1);
    for i = 1:n
        tgt = tc_m.Time(year(tc_m.Time)==year(fechas(i)) & month(tc_m.Time)==month(fechas(i)));
        if ~isempty(tgt)
            idx = find(tc_m.Time == tgt, 1);
            tc_dig_v(i) = tc_m.TC_DIGITAL(idx);
            tc_ref_v(i) = tc_m.TC_REFERENCIAL(idx);
        end
    end
    % Interpolacion lineal + forward-fill + fallback 6.96
    tc_dig_v = fillmissing(tc_dig_v, 'linear');
    tc_ref_v = fillmissing(tc_ref_v, 'linear');
    tc_dig_v = fillmissing(tc_dig_v, 'previous');
    tc_ref_v = fillmissing(tc_ref_v, 'previous');
    tc_dig_v(isnan(tc_dig_v)) = 6.96;
    tc_ref_v(isnan(tc_ref_v)) = 6.96;
    tc_ok = true;
catch ME
    warning('No se pudo procesar TC: %s — Modelo C omitido.', ME.message);
end

%% 2. Analisis preliminar

dlogIPC = diff(logIPC);

figure('Name','Preliminar SARIMA','Position',[50 50 1000 700]);
subplot(2,2,1); plot(fechas, logIPC, '-', 'Color',[0.17 0.24 0.31]);
title('log(IPC)'); grid on;
subplot(2,2,2); plot(fechas(2:end), dlogIPC, '-', 'Color',[0.16 0.50 0.73]); yline(0,'r--');
title('\Delta log(IPC)'); grid on;
subplot(2,2,3); autocorr(dlogIPC, 'NumLags', 36);
title('ACF \Delta log(IPC) — rezagos estacionales');
subplot(2,2,4); parcorr(dlogIPC, 'NumLags', 36);
title('PACF \Delta log(IPC)');

%% 3. Pruebas de raiz unitaria (incl. Zivot-Andrews)

fprintf('\n===== PRUEBAS DE RAIZ UNITARIA =====\n');
fprintf('  H0 (ADF/PP/ZA): raiz unitaria | H0 (KPSS): estacionaria\n\n');

[~, p_adf_log]   = adftest(logIPC);
[~, p_pp_log]    = pptest(logIPC);
[~, p_kpss_log]  = kpsstest(logIPC);
[~, p_adf_dlog]  = adftest(dlogIPC);
[~, p_pp_dlog]   = pptest(dlogIPC);
[~, p_kpss_dlog] = kpsstest(dlogIPC);

za_ok = false;
try
    [~, p_za_log,  ~, ~, bd_log]  = zatest(logIPC,  'Model','ARD');
    [~, p_za_dlog, ~, ~, bd_dlog] = zatest(dlogIPC, 'Model','ARD');
    za_ok = true;
catch ME
    warning('zatest: %s', ME.message);
    p_za_log = NaN; p_za_dlog = NaN;
end

fprintf('  %-24s  log(IPC)                     Dlog(IPC)\n','Prueba');
fprintf('  %s\n', repmat('-',1,72));
fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n','ADF', ...
    p_adf_log, ternary(p_adf_log<0.05,'[Estacionaria]  ','[Raiz unitaria]'), ...
    p_adf_dlog,ternary(p_adf_dlog<0.05,'[Estacionaria]','[Raiz unitaria]'));
fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n','Phillips-Perron', ...
    p_pp_log,  ternary(p_pp_log<0.05, '[Estacionaria]  ','[Raiz unitaria]'), ...
    p_pp_dlog, ternary(p_pp_dlog<0.05,'[Estacionaria]','[Raiz unitaria]'));
fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n','KPSS (H0: estac.)', ...
    p_kpss_log, ternary(p_kpss_log<0.05, '[Raiz unitaria] ','[Estacionaria] '), ...
    p_kpss_dlog,ternary(p_kpss_dlog<0.05,'[Raiz unitaria]','[Estacionaria]'));
if za_ok
    fprintf('  %-24s  p=%.4f  %-18s  p=%.4f  %s\n','Zivot-Andrews (*)', ...
        p_za_log, ternary(p_za_log<0.05, '[Estac.+quiebre]','[Raiz unitaria]'), ...
        p_za_dlog,ternary(p_za_dlog<0.05,'[Estac.+quiebre]','[Raiz unitaria]'));
    try
        fprintf('  Quiebre: log(IPC)=%s  |  Dlog(IPC)=%s\n', ...
            datestr(fechas(bd_log),'mmm yyyy'), ...
            datestr(fechas(min(bd_dlog+1,n)),'mmm yyyy'));
    catch; end
    fprintf('  (*) Estacionariedad condicionada al quiebre estructural detectado\n');
end

% Orden de integracion recomendado
d_rec = 1;
if p_adf_log < 0.05 && p_pp_log < 0.05 && p_kpss_log >= 0.05
    d_rec = 0;
    fprintf('\n  => log(IPC) aparentemente estacionaria — usando d=0\n');
else
    fprintf('\n  => log(IPC) integrada I(1) — usando d=1\n');
end

%% 4. Estacionalidad y determinacion de D

% ACF de diff(logIPC) en lag 12 para detectar si hace falta diferenciacion estacional
acf_vals = autocorr(dlogIPC, 'NumLags', 24);
acf_lag12 = acf_vals(13);  % indice 1=lag 0, 13=lag 12
umbral = 2 / sqrt(n-1);
D_rec = double(acf_lag12 > umbral);
fprintf('\n  ACF(lag=12) de Dlog(IPC) = %.4f (umbral=%.4f)\n', acf_lag12, umbral);
fprintf('  => D recomendado = %d  (la grilla evaluara D=0 y D=1)\n', D_rec);

%% 5. Dummies de intervencion

dummy_covid   = double(fechas >= datetime(2020,3,1) & fechas <= datetime(2020,8,1));
dummy_bloqueo = double(fechas >= datetime(2025,5,1) & fechas <= datetime(2025,6,1));
X_dummies = [dummy_covid, dummy_bloqueo];

%% 6. Grilla de estimacion: Modelos A, B, C

fprintf('\n===== ESTIMACION (p=0:%d, q=0:%d, P=0:%d, Q=0:%d, d=%d, D=0/1) =====\n', ...
    MAX_P, MAX_Q, MAX_PS, MAX_QS, d_rec);

X_C = [];
if tc_ok, X_C = [X_dummies, tc_dig_v, tc_ref_v]; end

fprintf('Modelo A (SARIMA puro)...\n');
[res_A, mdl_A] = estimar_sarima(logIPC, n, d_rec, MAX_P, MAX_Q, MAX_PS, MAX_QS, [], 'A');

fprintf('Modelo B (SARIMA + dummies)...\n');
[res_B, mdl_B] = estimar_sarima(logIPC, n, d_rec, MAX_P, MAX_Q, MAX_PS, MAX_QS, X_dummies, 'B');

if tc_ok
    fprintf('Modelo C (SARIMAX + TC)...\n');
    [res_C, mdl_C] = estimar_sarima(logIPC, n, d_rec, MAX_P, MAX_Q, MAX_PS, MAX_QS, X_C, 'C');
end

%% 7. Cuadros de scoring (Top 10 por modelo)

fprintf('\n--- Top 10 Modelo A ---\n');
disp(res_A(1:min(10,height(res_A)), {'Modelo','D','AICc','BIC','HQ','RMSE','MAPE','LB24','Score'}));
fprintf('Ganador A: %s (Score=%.4f)\n', mdl_A.nombre, mdl_A.score);

fprintf('\n--- Top 10 Modelo B ---\n');
disp(res_B(1:min(10,height(res_B)), {'Modelo','D','AICc','BIC','HQ','RMSE','MAPE','LB24','Score'}));
fprintf('Ganador B: %s (Score=%.4f)\n', mdl_B.nombre, mdl_B.score);

if tc_ok
    fprintf('\n--- Top 10 Modelo C ---\n');
    disp(res_C(1:min(10,height(res_C)), {'Modelo','D','AICc','BIC','HQ','RMSE','MAPE','LB24','Score'}));
    fprintf('Ganador C: %s (Score=%.4f)\n', mdl_C.nombre, mdl_C.score);
end

%% 8. Diagnostico de residuos (3x3 o 2x3)

n_mod = 2 + tc_ok;
figure('Name','Diagnosticos SARIMA','Position',[50 50 1000 700]);
mod_list  = {mdl_A, mdl_B}; if tc_ok, mod_list{3} = mdl_C; end
col_list  = {[0.15 0.68 0.38], [0.91 0.50 0.13], [0.56 0.27 0.68]};
lbl_list  = {'A: SARIMA', 'B: +Dummies', 'C: +TC'};
for mi = 1:n_mod
    mm = mod_list{mi}; cc = col_list{mi};
    subplot(n_mod,3,(mi-1)*3+1);
    plot(fechas, mm.res, '-', 'Color',cc); yline(0,'r--'); grid on;
    title(sprintf('%s: Residuos', lbl_list{mi}));
    subplot(n_mod,3,(mi-1)*3+2); autocorr(mm.res,'NumLags',24);
    title(sprintf('%s: ACF', lbl_list{mi}));
    subplot(n_mod,3,(mi-1)*3+3); qqplot(mm.res);
    title(sprintf('%s: Q-Q', lbl_list{mi}));
end

% Ljung-Box y Jarque-Bera para ganadores
for mi = 1:n_mod
    mm = mod_list{mi};
    nc = mm.n_coef;
    [~,p_lb] = lbqtest(mm.res,'Lags',24,'DoF',max(24-nc,1));
    [~,p_jb] = jbtest(mm.res);
    fprintf('%s: LB(24)=%.4f [%s]  JB=%.4f [%s]\n', lbl_list{mi}, ...
        p_lb, ternary(p_lb>0.05,'OK','AUTO'), p_jb, ternary(p_jb>0.05,'Normal','No-normal'));
end

%% 9. Backtesting (ultimos H_BACKTEST meses)

fprintf('\n===== BACKTESTING — ULTIMOS %d MESES =====\n', H_BACKTEST);
nt = n - H_BACKTEST;
ipc_obs_bt = IPC(nt+1:end);
fechas_bt  = fechas(nt+1:end);

ipc_bt = cell(n_mod,1);
for mi = 1:n_mod
    mm = mod_list{mi};
    X_mi = []; if mi==2, X_mi=X_dummies; elseif mi==3, X_mi=X_C; end
    ipc_bt{mi} = backtest_sarima(logIPC, nt, n, mm, X_mi, H_BACKTEST);
end

fprintf('\n  Fecha       IPC obs  ');
for mi=1:n_mod, fprintf('  %-12s', lbl_list{mi}); end
fprintf('\n  %s\n', repmat('-',1,50+13*n_mod));
for i=1:H_BACKTEST
    fprintf('  %-10s  %7.2f  ', datestr(fechas_bt(i),'mmm yyyy'), ipc_obs_bt(i));
    for mi=1:n_mod, fprintf('  %12.2f', ipc_bt{mi}(i)); end
    fprintf('\n');
end

fprintf('\n  %-16s  RMSE (IPC)   MAPE (%%)\n', 'Modelo');
fprintf('  %s\n', repmat('-',1,40));
for mi=1:n_mod
    rmse_i = sqrt(mean((ipc_bt{mi}-ipc_obs_bt).^2));
    mape_i = mean(abs((ipc_bt{mi}-ipc_obs_bt)./ipc_obs_bt))*100;
    fprintf('  %-16s  %10.4f   %8.2f%%\n', lbl_list{mi}, rmse_i, mape_i);
end

figure('Name','Backtesting','Position',[50 50 900 500]);
plot(fechas(1:nt), IPC(1:nt), '-', 'Color',[0.74 0.76 0.78], 'LineWidth',1); hold on;
plot(fechas_bt, ipc_obs_bt, '-o', 'Color',[0.17 0.24 0.31], 'LineWidth',1.5);
for mi=1:n_mod
    plot(fechas_bt, ipc_bt{mi}, '--o', 'Color',col_list{mi}, 'LineWidth',1.3);
end
legend([{'Entrenamiento','Observado'}, lbl_list(1:n_mod)], 'Location','northwest');
title('Backtesting IPC — Modelos SARIMA'); grid on; hold off;

%% 10. Pronostico final

fprintf('\n===== PRONOSTICO FINAL =====\n');
fecha_ultima = fechas(end);
fechas_pron  = dateshift(fecha_ultima + calmonths(1:H_FORECAST)', 'end','month');

% X futuros: dummies = 0, TC = ultimo valor disponible
X_fut_dummy = zeros(H_FORECAST, 2);
X_fut_C = [];
if tc_ok
    X_fut_C = [X_fut_dummy, repmat(tc_dig_v(end), H_FORECAST, 1), ...
                             repmat(tc_ref_v(end), H_FORECAST, 1)];
end

ipc_pron = cell(n_mod,1); inf_mens  = cell(n_mod,1); inf_ia = cell(n_mod,1);
meses_es = {'Enero','Febrero','Marzo','Abril','Mayo','Junio','Julio', ...
            'Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

for mi=1:n_mod
    mm = mod_list{mi};
    X_fut_mi = []; if mi==2, X_fut_mi=X_fut_dummy; elseif mi==3, X_fut_mi=X_fut_C; end
    [pron_ipc, pron_lo95, pron_hi95] = forecast_sarima(logIPC, IPC, fechas, n, mm, X_fut_mi, H_FORECAST);
    ipc_pron{mi} = pron_ipc;
    % Inflacion mensual
    inf_mens{mi} = (pron_ipc ./ [IPC(end); pron_ipc(1:end-1)] - 1) * 100;
    % Inflacion interanual
    ia = nan(H_FORECAST,1);
    for i=1:H_FORECAST
        f_base = fechas_pron(i) - calmonths(12);
        [~,ix] = min(abs(fechas - f_base));
        ia(i) = (pron_ipc(i)/IPC(ix) - 1)*100;
    end
    inf_ia{mi} = ia;
    if mi==1
        ipc_lo95_A = pron_lo95; ipc_hi95_A = pron_hi95;
    end
end

fprintf('\n  %-15s', 'Mes');
for mi=1:n_mod, fprintf('  %-22s', sprintf('%s (mens%%/ia%%)', lbl_list{mi})); end
fprintf('\n  %s\n', repmat('-',1,20+24*n_mod));
for i=1:H_FORECAST
    fprintf('  %-15s', sprintf('%s %d',meses_es{month(fechas_pron(i))},year(fechas_pron(i))));
    for mi=1:n_mod
        fprintf('  %7.4f%% / %6.2f%%  ', inf_mens{mi}(i), inf_ia{mi}(i));
    end
    fprintf('\n');
end

% Fan chart (Modelo A con IC, mas lineas B y C)
figure('Name','Pronostico SARIMA','Position',[50 50 1000 500]);
fecha_fan = [fecha_ultima; fechas_pron];
plot(fechas, IPC, '-', 'Color',[0.17 0.24 0.31], 'LineWidth',1.3); hold on;
fill([fechas_pron; flipud(fechas_pron)],[ipc_lo95_A; flipud(ipc_hi95_A)], ...
    col_list{1}, 'FaceAlpha',0.12,'EdgeColor','none');
for mi=1:n_mod
    plot(fecha_fan, [IPC(end);ipc_pron{mi}], '-', 'Color',col_list{mi}, 'LineWidth',1.5);
end
legend([{'Historico','IC95 A'}, lbl_list(1:n_mod)], 'Location','northwest');
title('Pronostico IPC — SARIMA / +Dummies / +TC'); grid on; hold off;

%% 11. Contrafactual (Modelo B vs. escenario post-bloqueos)

esc_pb = struct('fechas', [datetime(2026,7,1),datetime(2026,8,1),datetime(2026,9,1), ...
                            datetime(2026,10,1),datetime(2026,11,1),datetime(2026,12,1)], ...
                'ia',     [12.70, 13.40, 14.63, 15.38, 16.29, 16.90]);

idx_jul  = fechas_pron >= datetime(2026,7,1);
ia_B_all = inf_ia{2}(idx_jul);
n_cf     = min(length(ia_B_all), length(esc_pb.ia));  % trim al periodo del escenario
ia_B_fut = ia_B_all(1:n_cf);

fprintf('\n===== CONTRAFACTUAL (Modelo B) vs. POST-BLOQUEOS =====\n');
fprintf('  Mes         B (Contraf.)   Post-bloqueos  Dif (pp)\n');
fprintf('  %s\n', repmat('-',1,55));
for i = 1:length(esc_pb.ia)
    fprintf('  %-12s  %8.2f%%       %8.2f%%      %+.2f\n', ...
        meses_es{month(esc_pb.fechas(i))}, ia_B_fut(i), esc_pb.ia(i), ...
        esc_pb.ia(i) - ia_B_fut(i));
end
fprintf('  Brecha dic 2026: %.1f pp\n', esc_pb.ia(end) - ia_B_fut(end));

%% 12. Exportar resultados

if ~exist(fullfile(scriptDir,'..','Outputs'), 'dir')
    mkdir(fullfile(scriptDir,'..','Outputs'));
end

% Hoja proyeccion comparada
mes_pron_str = string(meses_es(month(fechas_pron)))';
T_proy = table(year(fechas_pron), mes_pron_str, ...
    round(inf_mens{1},4), round(inf_ia{1},2), ...
    round(inf_mens{2},4), round(inf_ia{2},2), ...
    'VariableNames',{'Anio','Mes','A_mens','A_ia','B_mens','B_ia'});
if tc_ok
    T_proy.C_mens = round(inf_mens{3},4);
    T_proy.C_ia   = round(inf_ia{3},2);
end
writetable(T_proy, RUTA_OUTPUT, 'Sheet','SARIMA_Proyeccion');

% Hoja contrafactual
fechas_cf  = esc_pb.fechas(1:n_cf)';   % 6×1 datetime
mes_pb_str = string(meses_es(month(fechas_cf)));
T_cf = table(year(fechas_cf), mes_pb_str, round(ia_B_fut,2), ...
    round(esc_pb.ia(1:n_cf)',2), round(esc_pb.ia(1:n_cf)' - ia_B_fut, 2), ...
    'VariableNames',{'Anio','Mes','B_contrafactual','Post_bloqueos','Dif_pp'});
writetable(T_cf, RUTA_OUTPUT, 'Sheet','SARIMA_Contrafactual');

% Hoja scoring resumida
if tc_ok
    T_sc = table({'A: SARIMA';'B: +Dummies';'C: +TC'}, ...
        {mdl_A.nombre; mdl_B.nombre; mdl_C.nombre}, ...
        [mdl_A.score; mdl_B.score; mdl_C.score], ...
        'VariableNames',{'Modelo','Ganador','Score'});
else
    T_sc = table({'A: SARIMA';'B: +Dummies'}, ...
        {mdl_A.nombre; mdl_B.nombre}, ...
        [mdl_A.score; mdl_B.score], ...
        'VariableNames',{'Modelo','Ganador','Score'});
end
writetable(T_sc, RUTA_OUTPUT, 'Sheet','SARIMA_Scoring');

fprintf('\nResultados exportados: %s\n', RUTA_OUTPUT);

%% =========================================================================
%% FUNCIONES LOCALES
%% =========================================================================

function [results, mejor] = estimar_sarima(logIPC, n, d_val, ...
        max_p, max_q, max_P, max_Q, X, ~)
% Grid search SARIMA(p,d,q)(P,D,Q)[12], D in {0,1}
% X: matrix de regresores (n x k) o [] para SARIMA puro
% Devuelve tabla de resultados ordenada por Score + struct del ganador

use_xreg = ~isempty(X);

Modelo=strings(0); pV=[]; dV=[]; qV=[]; PV=[]; DsV=[]; QV=[];
AICv=[]; BICv=[]; AICcv=[]; HQv=[]; RMSEv=[]; MAPEv=[];
LogLikv=[]; LB24v=[]; RaicesOk=logical([]); Cx=[];

for D_val = 0:1
    % Diferencia estacional (lag-12 manual, NO diff(y,12) que es de orden 12)
    if D_val == 1
        y  = logIPC(13:end) - logIPC(1:end-12);  % y(t) - y(t-12), len n-12
        ny = n - 12;
        X_y = []; if use_xreg, X_y = X(13:end,:); end
    else
        y = logIPC; ny = n; X_y = X;
    end

    for p = 0:max_p
     for q = 0:max_q
      for P = 0:max_P
       for Q = 0:max_Q
        if p==0&&q==0&&P==0&&Q==0, continue; end
        try
            % Argumentos de la parte no-estacional
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
            try, [~,p_lb] = lbqtest(res,'Lags',24,'DoF',max(24-k_coef,1));
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
            LogLikv(end+1,1)=logL; LB24v(end+1,1)=p_lb;
            RaicesOk(end+1,1)=rok; Cx(end+1,1)=p+q+P+Q;
        catch; continue; end
       end
      end
     end
    end
end

results = table(Modelo,pV,dV,qV,PV,DsV,QV,AICv,BICv,AICcv,HQv, ...
    RMSEv,MAPEv,LogLikv,LB24v,RaicesOk,Cx, ...
    'VariableNames',{'Modelo','p','d','q','P','D','Q','AIC','BIC','AICc','HQ', ...
                     'RMSE','MAPE','LogLik','LB24','Raices_ok','Complexity'});

% Scoring compuesto (identico a ARIMA_M2.m)
results.RMSE_s = zscore(results.RMSE); results.MAPE_s = zscore(results.MAPE);
results.AICc_s = zscore(results.AICc); results.BIC_s  = zscore(results.BIC);
results.HQ_s   = zscore(results.HQ);   results.Cx_s   = zscore(results.Complexity);
results.LB_pen = double(isnan(results.LB24) | results.LB24 < 0.05);
results.R_pen  = double(~results.Raices_ok);
results.Score  = results.RMSE_s + results.MAPE_s + ...
    (results.AICc_s+results.BIC_s+results.HQ_s)/3 + ...
    0.3*results.Cx_s + 2*results.LB_pen + 2*results.R_pen;
results = sortrows(results,'Score');

% Struct del ganador para uso en backtesting y pronostico
bst = results(1,:);
mejor.p = bst.p; mejor.d = bst.d; mejor.q = bst.q;
mejor.P = bst.P; mejor.D = bst.D; mejor.Q = bst.Q;
mejor.nombre = bst.Modelo{1}; mejor.score = bst.Score;

% Re-estimar ganador en muestra completa para guardar EstMdl y residuos
D_w = mejor.D;
if D_w == 1
    y_w = logIPC(13:end)-logIPC(1:end-12);
    X_w = []; if use_xreg, X_w=X(13:end,:); end
else
    y_w = logIPC; X_w = X;
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
mejor.y_fit  = y_w;  % serie usada en estimacion (para Y0 en forecast)
mejor.use_xreg = use_xreg;
end

% -------------------------------------------------------------------------
function ipc_bt = backtest_sarima(logIPC, nt, n, mm, X_mi, H_BT)
% Re-estima en muestra de entrenamiento y genera pronostico de backtest
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
    y_tr = logIPC(13:nt)-logIPC(1:nt-12);
    X_tr = []; X_te = [];
    if mm.use_xreg, X_tr=X_mi(13:nt,:); X_te=X_mi(nt+1:n,:); end
else
    y_tr = logIPC(1:nt);
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
        Ybt = forecast(EstMdl_bt, H_BT, 'Y0', y_tr, 'X0', X_te);
    end
    if D_w==1
        logIPC_bt = nan(H_BT,1);
        for i=1:H_BT, logIPC_bt(i)=logIPC(nt-12+i)+Ybt(i); end
    else
        logIPC_bt = Ybt;
    end
    ipc_bt = exp(logIPC_bt);
catch
    ipc_bt = nan(H_BT,1);
end
end

% -------------------------------------------------------------------------
function [ipc_p, lo95, hi95] = forecast_sarima(logIPC, IPC, ~, n, mm, X_fut, H)
% Genera pronostico H-pasos con IC 95%
D_w=mm.D; d_w=mm.d;
cst_w=ternary(d_w==0,NaN,0);
base_w={'ARLags',1:mm.p,'D',d_w,'MALags',1:mm.q,'Constant',cst_w};
seas_w={};
if mm.P>0||mm.Q>0
    seas_w={'Seasonality',12};
    if mm.P>0, seas_w=[seas_w,{'SARLags',1:mm.P}]; end
    if mm.Q>0, seas_w=[seas_w,{'SMALags',1:mm.Q}]; end
end
try
    if ~mm.use_xreg
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit);
    else
        [Yf, YMSEf] = forecast(mm.EstMdl, H, 'Y0', mm.y_fit, 'X0', X_fut);
    end
    if D_w==1
        logIPC_p=nan(H,1);
        for i=1:H, logIPC_p(i)=logIPC(n-12+i)+Yf(i); end
    else
        logIPC_p=Yf;
    end
    ipc_p   = exp(logIPC_p);
    lo95 = exp(logIPC_p - norminv(0.975)*sqrt(YMSEf));
    hi95 = exp(logIPC_p + norminv(0.975)*sqrt(YMSEf));
catch
    ipc_p=nan(H,1); lo95=nan(H,1); hi95=nan(H,1);
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out=a; else, out=b; end
end
