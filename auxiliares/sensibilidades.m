% ======================================
%  SENSIBILIDADES (DERIVADAS) PARA 3 TIPOS
% =======================================

% Parámetros base
S0   = 120;
T    = 10/12;
Smax = 360;
K0   = 100;
r0   = 0.10;
q0   = 0.25;
sg0  = 0.25;
N    = 100; 
M    = 100;

% Funciones constantes
rfunC    = @(rv) @(t) rv + 0*t;
qfunC    = @(qv) @(t) qv + 0*t;
sigmaC   = @(sv) @(t) sv + 0*t;
disc_r   = @(rv) exp(-rv*T);
disc_q   = @(qv) exp(-qv*T);

% ---- PRICERS (CN) ----
price_fix_put    = @(S0,K,r,q,sg) lookback_fixed_put (T,N,M, rfunC(r), qfunC(q), sigmaC(sg), K, S0, Smax);
price_fix_call   = @(S0,K,r,q,sg) lookback_fixed_call (T,N,M, rfunC(r), qfunC(q), sigmaC(sg), K, S0, Smax);
price_float_call = @(S0,r,q,sg) lookback_float_call (T,N,M, rfunC(r), qfunC(q), sigmaC(sg), S0, max(2, sg*10));
price_float_put  = @(S0,r,q,sg) lookback_float_put (T,N,M, rfunC(r), qfunC(q), sigmaC(sg), S0);

% ---- Diferencia centrada genérica ----
cent_diff = @(f, x, hx) (f(x+hx) - f(x-hx)) / (2*hx);

relStep = struct( ...
    'sigma', 1e-3, ...                 % relativo a sigma
    'r',     1e-3, ...                 % relativo a max(1, r)
    'q',     1e-3 ...                  % relativo a max(1, q)
);

% 1) VEGA: dV/dsigma en 5 puntos
sigGrid = linspace(0.12, 0.48, 5)';
dVsig_fix_put   = zeros(5,1);
dVsig_fix_call  = zeros(5,1);
dVsig_float_put = zeros(5,1);
dVsig_float_cal = zeros(5,1);

for i = 1:numel(sigGrid)
    
    sgi = sigGrid(i);
    hs  = max(1e-4, relStep.sigma*sgi);

    dVsig_fix_put(i) = cent_diff(@(x) price_fix_put (S0, K0, r0, q0, x), sgi, hs);
    dVsig_fix_call(i)= cent_diff(@(x) price_fix_call(S0, K0, r0, q0, x), sgi, hs);
    dVsig_float_put(i) = cent_diff(@(x) price_float_put (S0, r0, q0, x), sgi, hs);
    dVsig_float_cal(i) = cent_diff(@(x) price_float_call(S0, r0, q0, x), sgi, hs);
end

Tbl_VEGA = table(sigGrid, dVsig_fix_put, dVsig_fix_call, dVsig_float_put, dVsig_float_cal, ...
    'VariableNames', {'sigma','dV_dsigma_fix_put','dV_dsigma_fix_call','dV_dsigma_float_put','dV_dsigma_float_call'});
disp('== VEGA (dV/dsigma) =='); disp(Tbl_VEGA);

% 2) RHO: dV/dr en 5 puntos
rGrid = linspace(0.00, 0.20, 5)';
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

% 3) DIVIDEND-RHO: dV/dq en 5 puntos
qGrid = linspace(0.00, 0.10, 5)';
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
