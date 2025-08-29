%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluación de PUT lookback con strike fijo en la diagonal S = m
% ------------------------------------------------------------------------
% Resuelve la opción lookback tipo PUT con strike fijo usando Crank-Nicolson
% y devuelve L(0, S0, m0=S0) interpolando sobre la diagonal S = m.
%
% Entradas:
%   T      - tiempo hasta madurez
%   N, M   - número de pasos en tiempo y espacio
%   r,q    - tasas libres de riesgo y dividendos (pueden depender del tiempo)
%   sigma  - volatilidad (constante o función de tiempo)
%   K      - strike fijo
%   S0     - precio del subyacente para evaluación
%   Smax   - precio máximo considerado (opcional, defecto max(K, 2*S0))
%
% Salida:
%   L0     - precio interpolado en S0 sobre la diagonal S = m al tiempo 0
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L0 = lookback_fixed_put(T,N,M,r,q,sigma,K,S0,Smax)
    % Resuelve y evalúa L(0,S0,m0=S0) SOBRE LA DIAGONAL S=m (interpolación 1D).

    if nargin < 9 || isempty(Smax), Smax = max(K, 2*S0); end

    [L, t_vec, S, m] = fixed_put_cn(T,K,N,M,r,q,sigma,Smax);

    % Selecciona t=0
    [~, j0] = min(abs(t_vec));  % robusto

    % Extrae la DIAGONAL S=m en t=0: k = 1..Kmax
    Kmax = min(numel(S), numel(m));
    diagV = nan(1, Kmax);
    for k = 1:Kmax
        diagV(k) = L(j0, k, k);
    end

    % Interpola SOLO en la diagonal (S=m), sin NaNs de triángulo
    L0 = interp1(S(1:Kmax), diagV, S0, 'linear');
end
