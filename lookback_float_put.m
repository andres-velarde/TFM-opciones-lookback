%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   lookback_float_put.m
%   ---------------------------------------------------------------
%   Calcula el precio de una opción lookback PUT flotante a partir 
%   de la solución numérica W(t,xi) obtenida con el método de 
%   similaridad y esquema Crank–Nicolson.
%
%   INPUT:
%       T       - Horizonte temporal
%       N       - Número de pasos en tiempo
%       M       - Número de pasos en espacio (xi)
%       r(t)    - Función de tipo de interés
%       q(t)    - Función de dividendos
%       sigma(t)- Función de volatilidad
%       S_eval  - Precio spot evaluado (S_0)
%
%   OUTPUT:
%       precio  - Valor de la opción lookback PUT flotante en t=0
%
%   Relación con W:
%       precio = S_eval * W(0, xi=1)
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function precio = lookback_float_put(T, N, M, r, q, sigma, S_eval)
    [W, ~, ~] = floating_put_similaridad(T, N, M, r, q, sigma);
    precio = S_eval * W(1, end);  % xi=1 en t=0, esto es, J=S en tiempo inicial
end
