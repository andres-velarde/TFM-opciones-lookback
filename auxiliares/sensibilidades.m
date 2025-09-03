%% ======================================
%  SENSIBILIDADES (DERIVADAS) PARA 4 TIPOS
%  - Sólo derivadas (tablas compactas)
%  - 5 nodos por parámetro
%  - CN + paridades para flotantes
% =======================================

%% Parámetros base
S0   = 120;
K0   = 100;
T    = 10/12;
Smax = 360;
r0   = 0.10;
q0   = 0.25;
sg0  = 0.25;
N    = 150; 
M    = 150;

% Helpers constantes
rfunC    = @(rv) @(t) rv + 0*t;
qfunC    = @(qv) @(t) qv + 0*t;
sigmaC   = @(sv) @(t) sv + 0*t;
disc_r   = @(rv) exp(-rv*T);
disc_q   = @(qv) exp(-qv*T);

% ---- PRICERS (CN) ----
% Fixed PUT / CALL (asumen funciones existentes)
price_fix_put  = @(S0,K,r,q,sg) lookback_fixed_put (T,N,M, rfunC(r), qfunC(q), sigmaC(sg), K, S0, Smax);
price_fix_call = @(S0,K,r,q,sg) lookback_fixed_call(T,N,M, rfunC(r), qfunC(q), sigmaC(sg), K, S0, Smax);

% Floating via PARIDADES (sin solver extra)
price_float_call = @(S0,r,q,sg) ( S0*disc_q(q) - S0*disc_r(r) ) + price_fix_put (S0, S0, r, q, sg);
price_float_put  = @(S0,r,q,sg) (-S0*disc_q(q) + S0*disc_r(r) ) + price_fix_call(S0, S0, r, q, sg);

% ---- Diferencia centrada genérica (con paso relativo robusto) ----
cent_diff = @(f, x, hx) (f(x+hx) - f(x-hx)) / (2*hx);

relStep = struct( ...
    'sigma', 1e-3, ...                 % relativo a sigma
    'r',     1e-3, ...                 % relativo a max(1, r)
    'q',     1e-3, ...                 % relativo a max(1, q)
    'K',     1e-3  ...                 % relativo a K
);

%% 1) VEGA: dV/dsigma en 5 puntos
sigGrid = linspace(0.12, 0.48, 5)';    % evita rozar 0 y bordes
dVsig_fix_put   = zeros(5,1);
dVsig_fix_call  = zeros(5,1);
dVsig_float_put = zeros(5,1);
dVsig_float_cal = zeros(5,1);

for i = 1:numel(sigGrid)
    sgi = sigGrid(i);
    hs  = max(1e-4, relStep.sigma*sgi);

    % Fixed put / call
    dVsig_fix_put(i) = cent_diff(@(x) price_fix_put (S0, K0, r0, q0, x), sgi, hs);
    dVsig_fix_call(i)= cent_diff(@(x) price_fix_call(S0, K0, r0, q0, x), sgi, hs);

    % Floating put / call por paridad
    dVsig_float_put(i) = cent_diff(@(x) price_float_put (S0, r0, q0, x), sgi, hs);
    dVsig_float_cal(i) = cent_diff(@(x) price_float_call(S0, r0, q0, x), sgi, hs);
end

Tbl_VEGA = table(sigGrid, dVsig_fix_put, dVsig_fix_call, dVsig_float_put, dVsig_float_cal, ...
    'VariableNames', {'sigma','dV_dsigma_fix_put','dV_dsigma_fix_call','dV_dsigma_float_put','dV_dsigma_float_call'});
disp('== VEGA (dV/dsigma) =='); disp(Tbl_VEGA);

%% 2) RHO: dV/dr en 5 puntos
rGrid = linspace(0.00, 0.20, 5)';      % incluye 0
dVr_fix_put   = zeros(5,1);
dVr_fix_call  = zeros(5,1);
dVr_float_put = zeros(5,1);
dVr_float_cal = zeros(5,1);

for i = 1:numel(rGrid)
    ri = rGrid(i);
    hr = max(1e-5, relStep.r*max(1,abs(ri)));

    dVr_fix_put(i)   = cent_diff(@(x) price_fix_put (S0, K0, x, q0, sg0), ri, hr);
    dVr_fix_call(i)  = cent_diff(@(x) price_fix_call(S0, K0, x, q0, sg0), ri, hr);
    dVr_float_put(i) = cent_diff(@(x) price_float_put (S0, x,  q0, sg0), ri, hr);
    dVr_float_cal(i) = cent_diff(@(x) price_float_call(S0, x,  q0, sg0), ri, hr);
end

Tbl_RHO = table(rGrid, dVr_fix_put, dVr_fix_call, dVr_float_put, dVr_float_cal, ...
    'VariableNames', {'r','dV_dr_fix_put','dV_dr_fix_call','dV_dr_float_put','dV_dr_float_call'});
disp('== RHO (dV/dr) =='); disp(Tbl_RHO);

%% 3) DIVIDEND-RHO: dV/dq en 5 puntos
qGrid = linspace(0.00, 0.10, 5)';      % incluye 0
dVq_fix_put   = zeros(5,1);
dVq_fix_call  = zeros(5,1);
dVq_float_put = zeros(5,1);
dVq_float_cal = zeros(5,1);

for i = 1:numel(qGrid)
    qi = qGrid(i);
    hq = max(1e-5, relStep.q*max(1,abs(qi)));

    dVq_fix_put(i)   = cent_diff(@(x) price_fix_put (S0, K0, r0, x,  sg0), qi, hq);
    dVq_fix_call(i)  = cent_diff(@(x) price_fix_call(S0, K0, r0, x,  sg0), qi, hq);
    dVq_float_put(i) = cent_diff(@(x) price_float_put (S0, r0, x,  sg0), qi, hq);
    dVq_float_cal(i) = cent_diff(@(x) price_float_call(S0, r0, x,  sg0), qi, hq);
end

Tbl_QRHO = table(qGrid, dVq_fix_put, dVq_fix_call, dVq_float_put, dVq_float_cal, ...
    'VariableNames', {'q','dV_dq_fix_put','dV_dq_fix_call','dV_dq_float_put','dV_dq_float_call'});
disp('== DIVIDEND-RHO (dV/dq) =='); disp(Tbl_QRHO);

%% 4) dV/dK (sólo FIXED) en 5 puntos
KGrid = linspace(60, 120, 5)';         % 5 nodos de K
dVK_fix_put  = zeros(5,1);
dVK_fix_call = zeros(5,1);

for i = 1:numel(KGrid)
    Ki = KGrid(i);
    hK = max(1e-3, relStep.K*max(1,abs(Ki)));

    dVK_fix_put(i)  = cent_diff(@(x) price_fix_put (S0, x,  r0, q0, sg0), Ki, hK);
    dVK_fix_call(i) = cent_diff(@(x) price_fix_call(S0, x,  r0, q0, sg0), Ki, hK);
end

Tbl_dK = table(KGrid, dVK_fix_put, dVK_fix_call, ...
    'VariableNames', {'K','dV_dK_fix_put','dV_dK_fix_call'});
disp('== dV/dK (sólo fixed) =='); disp(Tbl_dK);
