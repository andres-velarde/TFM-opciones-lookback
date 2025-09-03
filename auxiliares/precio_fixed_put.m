%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Barrido de volatilidad para opción Lookback PUT de strike fijo
%   ----------------------------------------------------------------------
%   Calcula los precios de la opción usando:
%       - Fórmula exacta (analítica)
%       - Simulación Crank-Nicolson (backward PDE)
%       - Monte Carlo
%
%   Se generan columnas adicionales con el error relativo:
%       - Error respecto a la solución exacta (Crank-Nicolson)
%       - Error respecto a la solución exacta (Monte Carlo)
%
%   Parámetros fijos:
%       r_val, q_val   - tasas de interés y dividendos
%       T              - tiempo hasta madurez
%       N, M           - mallado temporal y espacial
%       S_eval         - precio del subyacente
%       K              - strike de la opción
%
%   Salida:
%       Tbl - tabla con resultados para cada sigma
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parámetros fijos
r_val = 0.1;
q_val = 0.5;
T = 12/12;
N = 100;
M = 100;
S_eval = 100;

% Rango de valores para sigma
sigmas = linspace(0.2, 0.5, 9);

% Inicializar arreglos para los resultados
num_sigmas = length(sigmas);
valores_simulacion = zeros(num_sigmas, 1);
mc_valores = zeros(num_sigmas, 1); % Monte Carlo
error_mc = zeros(num_sigmas, 1);   % Error Monte Carlo

% Bucle para cada valor de sigma
for i = 1:num_sigmas
    sig_val = sigmas(i);
    
    % Calcular la aproximación con Crank–Nicolson
    simulacion_val = lookback_fixed_put(T, N, M, @(t)r_val+t-t, @(t)q_val+t-t, ...
                                         @(t)sig_val+t-t, K, S_eval, max(5, S_eval * 2) * (1+sig_val));
    valores_simulacion(i) = simulacion_val;

    % Simulación Monte Carlo
    mc_valores(i) = mc_lookback_fixed_put(S_eval,K,T,10000,100000, ...
                     @(t)r_val+t-t, @(t)q_val+t-t, @(t)sig_val+t-t);
    
    % Error relativo [%] de MC respecto a la fórmula exacta
    if mc_valores(i) ~= 0
        error_mc(i) = abs(mc_valores(i) - simulacion_val) / mc_valores(i) * 100;
    else
        error_mc(i) = NaN; % Evitar división por cero
    end
    
end

% Crear la tabla de resultados
Tbl = table(sigmas', valores_simulacion, ...
            mc_valores, error_mc, ...
            'VariableNames', {'sigma', 'Crank-Nicolson', ...
                              'Monte Carlo', 'error_mc'});

% Mostrar la tabla
disp(Tbl);

