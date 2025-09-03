%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Experimento de validación para Lookback Call flotante
%   ---------------------------------------------------------------
%   Este script compara:
%       • Fórmula exacta (monitorización continua, coeficientes constantes)
%       • Aproximación numérica (Crank–Nicolson en PDE)
%       • Estimación por Monte Carlo discreto
%
%   Se calcula el error relativo porcentual respecto de la fórmula exacta
%   para diferentes valores de sigma.
%
%   Columnas en la tabla de salida:
%       sigma           - Volatilidad
%       exacta          - Valor analítico (fórmula cerrada)
%       Crank-Nicolson  - Valor obtenido por PDE (CN)
%       error           - Error relativo [%] de CN vs exacta
%       Monte Carlo     - Valor obtenido por simulación MC
%       error_mc        - Error relativo [%] de MC vs exacta
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
    simulacion_val = lookback_float_call(T, N, M, @(t)r_val+t-t, @(t)q_val+t-t, ...
                                         @(t)sig_val+t-t, S_eval, max(2, sig_val*10));
    valores_simulacion(i) = simulacion_val;

    % Simulación Monte Carlo
    mc_valores(i) = mc_lookback_floating_call(S_eval,T,10000,100000, ...
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

