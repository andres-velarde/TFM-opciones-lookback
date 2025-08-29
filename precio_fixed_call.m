%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluación de precios de opciones PUT lookback con strike fijo
% ------------------------------------------------------------------------
% Este script calcula los precios de un lookback CALL de strike fijo usando:
%   1. Fórmula analítica (monitorización continua)
%   2. Método Crank-Nicolson (discreto)
%   3. Simulación Monte Carlo
%
% También calcula los errores relativos porcentuales respecto a la fórmula
% exacta para cada método y presenta los resultados en una tabla.
%
% Parámetros:
%   r_val, q_val: tasas de interés y dividendos
%   T           : madurez en años
%   N, M        : número de pasos en tiempo y espacio
%   S_eval      : precio inicial del subyacente
%   K           : strike fijo de la opción
%   sigmas      : vector de volatilidades a evaluar
%
% Salida:
%   Tbl         : tabla con columnas:
%                 sigma | exacta | Crank-Nicolson | error | Monte Carlo | error_mc
%
% Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parámetros fijos
r_val = 0.1;
q_val = 0.5;
T = 12/12;
N = 200;
M = 200;
S_eval = 100;
K = 80;

% Rango de valores para sigma
sigmas = linspace(0.2, 0.5, 9);

% Inicializar arreglos para los resultados
num_sigmas = length(sigmas);
valores_simulacion = zeros(num_sigmas, 1);
valores_exactos = zeros(num_sigmas, 1);
errores = zeros(num_sigmas, 1);
mc_valores = zeros(num_sigmas, 1); % Monte Carlo
error_mc = zeros(num_sigmas, 1);    % Error Monte Carlo

% Bucle para cada valor de sigma
for i = 1:num_sigmas
    sig_val = sigmas(i);

    % Calcular el precio exacto (fórmula analítica)
    exacto_val = lookback_fixed_call_formulae(S_eval, K, r_val, q_val, sig_val, T);
    valores_exactos(i) = exacto_val;
   
    % Calcular el precio mediante Crank-Nicolson
    simulacion_val = lookback_fixed_call(T, N, M, @(t)r_val+t-t, @(t)q_val+t-t, @(t)sig_val+t-t, K, S_eval, 200);
    valores_simulacion(i) = simulacion_val;

    % Error relativo porcentual para Crank-Nicolson
    errores(i) = abs(valores_simulacion(i) - exacto_val) / abs(exacto_val) * 100;

    % ---------------------------------------------------------------------
    % Calcular precio con Monte Carlo y error relativo
    % ---------------------------------------------------------------------
    mc_valores(i) = mc_lookback_fixed_call(S_eval,K, T,1000,1000,@(t)r_val+t-t, @(t)q_val+t-t, @(t)sig_val+t-t);
    if exacto_val ~= 0
        error_mc(i) = abs(mc_valores(i) - exacto_val) / mc_valores(i) * 100;
    else
        errores(i) = NaN; % Evitar división por cero
    end
end

% Crear la tabla de resultados
Tbl = table(sigmas', valores_exactos, valores_simulacion, errores, mc_valores, error_mc, ...
            'VariableNames', {'sigma', 'exacta', 'Crank-Nicolson', 'error', 'Monte Carlo', 'error_mc'});

% Mostrar la tabla
disp(Tbl);
