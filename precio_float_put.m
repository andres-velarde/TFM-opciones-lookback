%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Experimento numérico: Lookback PUT flotante
%   ---------------------------------------------------------------
%   Compara los precios obtenidos con:
%       1. Fórmula analítica cerrada
%       2. Esquema Crank–Nicolson (PDE con similaridad)
%       3. Simulación Monte Carlo
%
%   Se calcula además el error relativo (%) frente a la solución exacta.
%
%   Columnas finales de la tabla:
%       sigma           - Volatilidad
%       exacta          - Valor por fórmula analítica
%       Crank-Nicolson  - Valor por PDE
%       error           - Error relativo (%) CN vs exacta
%       Monte Carlo     - Valor por simulación MC
%       error_mc        - Error relativo (%) MC vs exacta
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parámetros fijos
r_val = 0.1;       % tipo de interés
q_val = 0.5;       % dividendos
T = 12/12;         % horizonte temporal (1 año)
N = 100;           % pasos en tiempo (CN)
M = 100;           % pasos en espacio (CN)
S_eval = 100;      % precio spot inicial

% Rango de valores para sigma
sigmas = linspace(0.2, 0.5, 9);

% Inicializar arreglos para los resultados
num_sigmas = length(sigmas);
valores_simulacion = zeros(num_sigmas, 1);  % CN
valores_exactos    = zeros(num_sigmas, 1);  % fórmula
errores            = zeros(num_sigmas, 1);  % error CN
mc_valores         = zeros(num_sigmas, 1);  % Monte Carlo
error_mc           = zeros(num_sigmas, 1);  % error MC

% Bucle para cada valor de sigma
for i = 1:num_sigmas
    sig_val = sigmas(i);

    % --- (1) Precio exacto (fórmula analítica)
    exacto_val = lookback_float_put_formulae(S_eval, r_val, q_val, sig_val, T);
    valores_exactos(i) = exacto_val;
   
    % --- (2) Precio con CN
    simulacion_val = lookback_float_put(...
        T, N, M, @(t)r_val+t-t, @(t)q_val+t-t, @(t)sig_val+t-t, S_eval);
    valores_simulacion(i) = simulacion_val;
    
    % Error relativo CN
    errores(i) = abs(valores_simulacion(i) - exacto_val) / abs(exacto_val) * 100;

    % --- (3) Precio con Monte Carlo
    mc_valores(i) = mc_lookback_floating_put(...
        S_eval, T, 1000, 1000, @(t)r_val+t-t, @(t)q_val+t-t, @(t)sig_val+t-t);

    % Error relativo MC
    if exacto_val ~= 0
        error_mc(i) = abs(mc_valores(i) - exacto_val) / abs(exacto_val) * 100;
    else
        error_mc(i) = NaN; % evitar división por cero
    end
end

% Crear la tabla de resultados
Tbl = table(sigmas', valores_exactos, valores_simulacion, errores, ...
            mc_valores, error_mc, ...
    'VariableNames', {'sigma', 'exacta', 'Crank_Nicolson', 'error', ...
                      'Monte_Carlo', 'error_mc'});

% Mostrar la tabla en consola
disp(Tbl);
