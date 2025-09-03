%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fórmula exacta: Lookback PUT con strike fijo (monitorización continua)
% ------------------------------------------------------------------------
% Calcula el precio de una opción lookback tipo PUT con strike fijo K
% bajo el modelo Black-Scholes con tasas y volatilidad constantes.
%
% Asume que el mínimo inicial es el precio inicial S0 (sin histórico previo).
%
% Entradas:
%   S0    - precio inicial del subyacente
%   K     - strike fijo
%   r     - tasa libre de riesgo
%   q     - tasa de dividendos
%   sigma - volatilidad del subyacente
%   T     - tiempo hasta madurez
%
% Salida:
%   P     - precio exacto de la opción
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = lookback_fixed_put_formulae(S0,K,r,q,sigma,T)
    b   = r - q;
    dt  = sigma*sqrt(T);
    d1  = (log(S0./K) + (b + 0.5*sigma^2)*T) / dt;

    if abs(b) < 1e-12
        % Límite b -> 0
        phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
        P = K.*exp(-r*T).*normcdf(-d1 + dt) ...
          - S0.*exp(-q*T).*normcdf(-d1) ...
          + 0.5*sigma^2*T * S0.*exp(-q*T).*phi(-d1);
    else
        lam = 2*b/(sigma^2);
        P = K .*exp(-r*T).*normcdf(-d1 + dt) ...
          - S0.*exp(-q*T).*normcdf(-d1) ...
          + (sigma^2/(2*b))*S0 .* ( ...
                exp(-r*T).*(S0./K).^(-lam) .* normcdf(-d1 + lam*dt) ...
              - exp(-q*T).*normcdf(-d1) );
    end
end
