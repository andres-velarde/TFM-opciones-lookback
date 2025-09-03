%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fórmula exacta: Lookback PUT con strike fijo (monitorización continua)
% ------------------------------------------------------------------------
% Payoff: max(K - M_T, 0), con M_T = min_{0<=u<=T} S_u.
% Se asume que el mínimo inicial es M0 = S0.
%
% Entradas:
%   S0    - spot inicial
%   K     - strike fijo
%   r     - tasa libre de riesgo
%   q     - tasa de dividendos
%   sigma - volatilidad
%   T     - madurez
%
% Salida:
%   P     - precio exacto
%
% Ref: Conze & Viswanathan (1991), Goldman–Sosin–Gatto (1979)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = lookback_fixed_put_formulae(S0,K,r,q,sigma,T)

    b  = r - q;
    dt = sigma*sqrt(T);

    d1 = (log(S0./K) + (b + 0.5*sigma^2)*T) / dt;
    d2 = d1 - dt;

    if abs(b) < 1e-12
        % Límite b -> 0
        phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
        P = K*exp(-r*T).*normcdf(-d2) ...
          - S0*exp(-q*T).*normcdf(-d1) ...
          + 0.5*sigma^2*T*S0*exp(-q*T).*phi(-d1);
    else
        lam = 2*b/sigma^2;
        P = K*exp(-r*T).*normcdf(-d2) ...
          - S0*exp(-q*T).*normcdf(-d1) ...
          + S0*exp(-q*T)*(sigma^2/(2*b)) .* ( ...
                (S0./K).^(-lam) .* normcdf(-d1 + lam*dt) ...
              - exp(b*T).*normcdf(-d1) );
    end
end
