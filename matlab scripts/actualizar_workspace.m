%% actualizar_workspace.m — Re-importa BD Bruta.xlsx y actualiza el workspace
% Ejecutar cada vez que se agreguen nuevas observaciones a BD Bruta.xlsx.
% IMPORTANTE: BD Bruta.xlsx debe estar CERRADO antes de correr este script.

clc;
RUTA_BD = 'C:\Users\Juande\Documents\Scripts Python\TC e π\Inputs\BD Bruta.xlsx';
RUTA_WS = 'C:\Users\Juande\Documents\Scripts Python\TC e π\matlab scripts\matlab workspace ARIMA.mat';

meses_nom = {'Enero','Febrero','Marzo','Abril','Mayo','Junio', ...
             'Julio','Agosto','Septiembre','Octubre','Noviembre','Diciembre'};

%% 1. Importar IPC (columnas por posicion: [1]=ano, [2]=mes_nombre, [3]=valor)
fprintf('=== Importando hoja IPC ===\n');
T_ipc = readtable(RUTA_BD, 'Sheet', 'IPC');

anio_raw = T_ipc{:,1};  % double (ano)
mes_raw  = T_ipc{:,2};  % string / cell (nombre del mes)
val_raw  = T_ipc{:,3};  % double o mixed (valor IPC)

% Convertir valor a numerico de forma robusta (tolera celdas vacias y strings)
if iscell(val_raw)
    val_num = cellfun(@(x) str2double(string(x)), val_raw);
elseif isstring(val_raw) || ischar(val_raw)
    val_num = str2double(string(val_raw));
else
    val_num = double(val_raw);
end

valid    = ~isnan(val_num) & ~isnan(double(anio_raw));
n_ok     = sum(valid);
years_ok = double(anio_raw(valid));
mes_ok   = string(mes_raw(valid));
vals_ok  = val_num(valid);

% Construir vector de fechas (primer dia de cada mes)
fechas_ipc = NaT(n_ok, 1);
for i = 1:n_ok
    m_idx = find(strcmp(meses_nom, char(mes_ok(i))), 1);
    if ~isempty(m_idx)
        fechas_ipc(i) = datetime(years_ok(i), m_idx, 1);
    end
end

% Timetable con las mismas columnas que el workspace original
BDBruta_IPC = timetable(fechas_ipc, years_ok, vals_ok, 'VariableNames', {'A_o','Valor'});
BDBruta_IPC.Properties.DimensionNames{1} = 'Time';
BDBruta_IPC = sortrows(BDBruta_IPC, 'Time');
BDBruta_IPC = BDBruta_IPC(~isnat(BDBruta_IPC.Time), :);  % purgar NaT

fprintf('  IPC importado: %s → %s (n=%d)\n', ...
    datestr(BDBruta_IPC.Time(1),'mmm yyyy'), ...
    datestr(BDBruta_IPC.Time(end),'mmm yyyy'), height(BDBruta_IPC));
fprintf('  Ultimo valor: IPC = %.4f (%s %d)\n', ...
    BDBruta_IPC.Valor(end), meses_nom{month(BDBruta_IPC.Time(end))}, ...
    year(BDBruta_IPC.Time(end)));

%% 2. Importar TC (columnas: [1]=dia_habil, [2]=fecha, [3..6]=TC_PREF/DIG/OFIC/REF)
fprintf('\n=== Importando hoja TC ===\n');
T_tc = readtable(RUTA_BD, 'Sheet', 'TC');

% Columna fecha (col 2) — puede llegar como datetime, numeric (serial Excel) o cell
fecha_raw = T_tc{:,2};
if isdatetime(fecha_raw)
    fechas_tc = fecha_raw;
elseif isnumeric(fecha_raw)
    fechas_tc = datetime(fecha_raw, 'ConvertFrom', 'excel');
elseif iscell(fecha_raw)
    try
        fechas_tc = datetime(fecha_raw);
    catch
        fechas_tc = NaT(height(T_tc), 1);
        for i = 1:height(T_tc)
            try, fechas_tc(i) = datetime(fecha_raw{i}); catch, end
        end
    end
else
    fechas_tc = datetime(fecha_raw);
end

% Columnas TC — leer como double (NaN para celdas vacias)
safe_dbl = @(col) str2double(string(col));  % tolera cualquier tipo entrada
pref = safe_dbl(T_tc{:,3});
dig  = safe_dbl(T_tc{:,4});
ofic = safe_dbl(T_tc{:,5});
ref  = safe_dbl(T_tc{:,6});

% Solo se exige fecha valida: TC PREFERENCIAL puede faltar en filas nuevas
% (p.ej. actualizaciones parciales desde un CSV que solo trae oficial/digital);
% exigir pref no-NaN aqui las descartaba silenciosamente aunque tuvieran
% TC-DIGITAL/TC-OFICIAL validos.
valid_tc = ~isnat(fechas_tc);
BDBruta_TC = timetable(fechas_tc(valid_tc), pref(valid_tc), dig(valid_tc), ...
    ofic(valid_tc), ref(valid_tc), ...
    'VariableNames', {'TC_PREFERENCIAL','TC_DIGITAL','TC_OFICIAL','TC_REFERENCIAL'});
BDBruta_TC.Properties.DimensionNames{1} = 'Time';
BDBruta_TC = sortrows(BDBruta_TC, 'Time');

fprintf('  TC importado : %s → %s (n=%d dias)\n', ...
    datestr(BDBruta_TC.Time(1),'dd/mm/yyyy'), ...
    datestr(BDBruta_TC.Time(end),'dd/mm/yyyy'), height(BDBruta_TC));

%% 3. Guardar workspace (sobreescribe IPC y TC, preserva otras variables si existen)
if isfile(RUTA_WS)
    save(RUTA_WS, 'BDBruta_IPC', 'BDBruta_TC', '-append');
else
    save(RUTA_WS, 'BDBruta_IPC', 'BDBruta_TC');
end
fprintf('\nWorkspace actualizado en: %s\n', RUTA_WS);
