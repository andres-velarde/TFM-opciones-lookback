%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Precio Lookback CALL flotante
%   ----------------------------------------------------------------------
%   Evalúa el precio de un lookback CALL flotante en t=0 usando la
%   transformación de variables de similaridad.
%
%   Entrada:
%       T       - Tiempo hasta madurez
%       N, M    - Número de pasos temporales y espaciales
%       r, q    - Funciones de tasa libre de riesgo y dividendos (dependientes de t)
%       sigma   - Función de volatilidad (dependiente de t)
%       S_eval  - Valor del activo subyacente al inicio
%       xi_max  - Valor máximo para la transformación de variables xi
%
%   Salida:
%       precio  - Precio de la opción lookback CALL flotante
%
%   Nota: Se evalúa en t=0 y xi=1 (primer nodo tras reconstrucción)
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function precio = lookback_float_call(T, N, M, r, q, sigma, S_eval, xi_max)

    % Llamada a la función que resuelve la PDE
    [W, ~, ~] = floating_call_similaridad(T, N, M, r, q, sigma, xi_max);

    % Evaluación en t=0, xi=1 (fila 1, columna 1 tras reconstrucción)
    precio = S_eval * W(1, 1);

end
