%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluación de CALL lookback con strike fijo en la diagonal S = J
% ------------------------------------------------------------------------
% Resuelve la opción lookback tipo CALL con strike fijo usando Crank–Nicolson
% y devuelve C(0, S0, J0=S0) interpolando sobre la diagonal S = J.
%
% Entradas:
%   T      - tiempo hasta madurez
%   N, M   - nº de pasos en tiempo y espacio
%   r,q    - tipos libre de riesgo y dividendos (funciones o constantes)
%   sigma  - volatilidad (función o constante)
%   K      - strike fijo
%   S0     - spot para evaluar
%   Smax   - cota superior del dominio (opcional; por defecto max(K, 2*S0))
%
% Salida:
%   C0     - precio interpolado en S0 sobre la diagonal S = J al tiempo 0
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C0 = lookback_fixed_call(T,N,M,r,q,sigma,K,S0,Smax)
    % Si no se pasa Smax, usa una cota razonable
    if nargin < 9 || isempty(Smax), Smax = max(K, 2*S0); end

    % --- RESUELVE el problema 3D para la CALL fija (tu solver simétrico del PUT)
    % Debes tener implementado 'fixed_call_cn' con la misma interfaz que 'fixed_put_cn':
    %   [L, t_vec, S_vec, J_vec] = fixed_call_cn(T,K,N,M,r,q,sigma,Smax)
    [L, t_vec, S, J] = fixed_call_cn(T, K, N, M, r, q, sigma, Smax);

    % --- Selecciona t = 0 (índice temporal)
    [~, j0] = min(abs(t_vec));   % robusto por si t_vec no pasa EXACTO por 0

    % --- Extrae la DIAGONAL S = J en t = 0
    Kmax  = min(numel(S), numel(J));
    diagV = nan(1, Kmax);
    for k = 1:Kmax
        diagV(k) = L(j0, k, k);
    end

    % --- Interpolación 1D SOLO sobre la diagonal (evita NaNs del triángulo)
    %     (lineal; sin extrapolar fuera del rango de S)
    C0 = interp1(S(1:Kmax), diagV, S0, 'linear');
end
