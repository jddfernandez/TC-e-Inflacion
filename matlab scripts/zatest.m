function [h, pValue, stat, critValue, breakDate] = zatest(y, varargin)
%ZATEST  Zivot-Andrews (1992) unit root test with unknown structural break.
%
%  [h, pValue, stat, critValue, breakDate] = zatest(y)
%  [h, pValue, stat, critValue, breakDate] = zatest(y, 'Model', m)
%
%  Tests H0: unit root (no break) vs H1: stationary with one structural break.
%  Same output interface as MATLAB's built-in zatest (Econometrics Toolbox).
%
%  Model options (MATLAB toolbox convention → ZA 1992 notation):
%    'ARD' (default) — break in intercept + trend  (Model C)
%    'LS'            — break in intercept only      (Model A)
%    'TS'            — break in trend only          (Model B)
%
%  Lag selection: AIC over 0..floor(12*(n/100)^0.25) at each break date.
%  Trimming: 10% from each end (default).

% --- Parse inputs ---
mdl_arg = 'ARD';
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i}, 'Model')
        mdl_arg = upper(varargin{i+1});
    end
end

% Map toolbox convention → ZA 1992 model letter
switch mdl_arg
    case 'ARD', za_model = 'C';
    case 'LS',  za_model = 'A';
    case 'TS',  za_model = 'B';
    otherwise,  error('zatest:badModel', 'Model must be ''ARD'', ''LS'', or ''TS''.');
end

y      = y(:);
y      = y(~isnan(y));
n      = numel(y);
maxLag = min(floor(12*(n/100)^0.25), max(0, floor((n-8)/2)));
trim   = 0.10;

firstTb = max(3, floor(trim*n));
lastTb  = min(n-2, ceil((1-trim)*n));

if firstTb >= lastTb
    error('zatest:trim', 'Series too short for 10%% trimming (n=%d).', n);
end

% --- Grid over break dates (AIC lag selection at each tb) ---
bestStat = Inf;
bestTb   = firstTb;
bestLag  = 0;

for tb = firstTb:lastTb
    bestAIC = Inf;
    for k = 0:maxLag
        reg = za_reg(y, tb, k, za_model);
        if isempty(reg), continue; end
        if reg.aic < bestAIC
            bestAIC  = reg.aic;
            best_reg = reg;
        end
    end
    if ~exist('best_reg','var') || isempty(best_reg), continue; end
    if best_reg.tAlpha < bestStat
        bestStat = best_reg.tAlpha;
        bestTb   = tb;
        bestLag  = best_reg.lag;  %#ok<NASGU>
    end
    clear best_reg
end

stat      = bestStat;
breakDate = bestTb;

% --- Critical values (Zivot & Andrews 1992, Table 2, T→∞) ---
switch za_model
    case 'C', cv = struct('cv1',-5.57, 'cv5',-5.08, 'cv10',-4.82);
    case 'A', cv = struct('cv1',-5.34, 'cv5',-4.80, 'cv10',-4.58);
    case 'B', cv = struct('cv1',-4.93, 'cv5',-4.42, 'cv10',-4.11);
end
critValue = cv.cv5;

% --- p-value (linear interpolation, capped at [0.01, 0.10]) ---
pValue = za_pval(stat, cv);

h = double(stat < critValue);

end


% =========================================================================
function reg = za_reg(y, tb, k, model)
% ADF regression with structural-break dummies at break index tb.
% Column order: [const, y_{t-1}, trend, dummies..., lag_diffs...]
% t-stat on y_{t-1} is beta(2)/se(2).

    dy  = diff(y);
    n   = numel(y);
    idx = (k+2 : n)';

    if numel(idx) < 8
        reg = []; return;
    end

    Y    = dy(idx-1);
    yLag = y(idx-1);
    tr   = idx;
    DU   = double(idx > tb);
    DT   = (idx - tb) .* DU;

    switch upper(model)
        case 'C', X = [ones(numel(idx),1), yLag, tr, DU, DT];
        case 'A', X = [ones(numel(idx),1), yLag, tr, DU];
        case 'B', X = [ones(numel(idx),1), yLag, tr, DT];
        otherwise, reg = []; return;
    end

    for j = 1:k
        X = [X, dy(idx-1-j)]; %#ok<AGROW>
    end

    nObs = size(X,1);
    nPar = size(X,2);
    if nObs <= nPar + 2 || rank(X) < nPar
        reg = []; return;
    end

    beta = X \ Y;
    e    = Y - X*beta;
    s2   = (e'*e) / (nObs - nPar);
    se   = sqrt(diag(s2 * ((X'*X) \ eye(nPar))));

    reg       = struct();
    reg.tAlpha = beta(2) / se(2);
    reg.aic    = nObs * log((e'*e)/nObs) + 2*nPar;
    reg.lag    = k;
end


% =========================================================================
function p = za_pval(stat, cv)
% Linear interpolation between 1%, 5%, 10% critical values.

    if stat <= cv.cv1
        p = 0.01;
    elseif stat <= cv.cv5
        p = 0.01 + 0.04 * (stat - cv.cv1) / (cv.cv5  - cv.cv1);
    elseif stat <= cv.cv10
        p = 0.05 + 0.05 * (stat - cv.cv5) / (cv.cv10 - cv.cv5);
    else
        p = 0.10;
    end
end
