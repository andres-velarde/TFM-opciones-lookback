%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Comparación gráfica: Opciones Lookback PUT (flotante vs fija)
%   -----------------------------------------------------------------------
%   Se muestran en una única figura (subplots):
%     • Izquierda: lookback PUT flotante
%     • Derecha:  lookback PUT fija
%
%   Para cada caso se comparan tres métodos:
%     - Fórmula cerrada (coeficientes constantes)
%     - Aproximación por Crank–Nicolson (diferencias finitas)
%     - Simulación de Monte Carlo
%
%   Parámetros de entrada (constantes en el tiempo):
%       r_val   - tipo de interés libre de riesgo
%       q_val   - tasa de dividendos
%       sig_val - volatilidad
%       T       - madurez (años)
%       K       - strike fijo (solo para el caso fijo)
%       Smax    - cota superior en malla espacial (caso fijo)
%       N, M    - nodos temporales y espaciales
%
%   Salida:
%       Una figura con dos subplots comparando los tres métodos.
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ========= Parámetros comunes =========
r_val   = 0.10;
q_val   = 0.25;
sig_val = 0.10;
T       = 10/12;
N       = 150;
M       = 150;

% Constantes (handles triviales para el solver)
rfun     = @(t) r_val+t-t;
qfun     = @(t) q_val+t-t;
sigmafun = @(t) sig_val+t-t;

% ========= Mallado de S0 =========
S_vals = linspace(0.01, 100, 80);   % evita S=0
nS = length(S_vals);

% ========= Caso FLOTTANTE =========
val_formula_float = zeros(1,nS);
val_cn_float      = zeros(1,nS);
val_mc_float      = zeros(1,nS);

for i=1:nS
    S0 = S_vals(i);
    % Fórmula cerrada (flotante put)
    val_formula_float(i) = lookback_float_put_formulae(S0, r_val, q_val, sig_val, T);
    % Diferencias finitas (CN)
    val_cn_float(i)      = lookback_float_put(T, N, M, rfun, qfun, sigmafun, S0);
    % Monte Carlo
    val_mc_float(i)      = mc_lookback_floating_put(S0, T, 1000, 10000, rfun, qfun, sigmafun);
end

% ========= Caso FIJO =========
K     = 60;
S_vals_fix = linspace(0.01, 100, 60);
nS_fix = length(S_vals_fix);

val_formula_fix = zeros(1,nS_fix);
val_cn_fix      = zeros(1,nS_fix);
val_mc_fix      = zeros(1,nS_fix);

for i=1:nS_fix
    S0 = S_vals_fix(i);
    % Fórmula cerrada (fijo put)
    val_formula_fix(i) = lookback_fixed_put_formulae(S0, K, r_val, q_val, sig_val, T);
    % Diferencias finitas (CN)
    val_cn_fix(i)      = lookback_fixed_put(T, N, M, rfun, qfun, sigmafun, K, S0, max(S0*5, 10)*max(sig_val*5, 1));
    % Monte Carlo
    val_mc_fix(i)      = mc_lookback_fixed_put(S0, K, T, 1000, 10000, rfun, qfun, sigmafun);
end

% ========= Gráfico combinado =========
figure;

% Subplot izquierda: flotante
subplot(1,2,1);
plot(S_vals, val_formula_float, 'k-',  'LineWidth', 2); hold on;
plot(S_vals, val_cn_float,      'b--', 'LineWidth', 1.5);
plot(S_vals, val_mc_float,      'r:',  'LineWidth', 1.5);
xlabel('Precio spot S_0');
ylabel('Precio opción lookback PUT flotante');
title('Lookback flotante (PUT)');
legend('Fórmula cerrada','Crank–Nicolson','Monte Carlo','Location','northwest');
grid on;

% Subplot derecha: fija
subplot(1,2,2);
plot(S_vals_fix, val_formula_fix, 'k-',  'LineWidth', 2); hold on;
plot(S_vals_fix, val_cn_fix,      'b--', 'LineWidth', 1.5);
plot(S_vals_fix, val_mc_fix,      'r:',  'LineWidth', 1.5);
xlabel('Precio spot S_0');
ylabel('Precio opción lookback PUT fija');
title('Lookback fija (PUT, K=60)');
legend('Fórmula cerrada','Crank–Nicolson','Monte Carlo','Location','northwest');
grid on;
