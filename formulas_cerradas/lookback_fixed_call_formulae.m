%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lookback CALL de strike fijo (Conze–Viswanathan)
% ------------------------------------------------------------------------
% Fórmula exacta con monitorización continua.
%
% Entradas:
%   S0    : precio inicial del subyacente
%   K     : strike fijo
%   r     : tasa libre de riesgo constante
%   q     : tasa de dividendos constante
%   sigma : volatilidad constante
%   T     : tiempo hasta madurez
%
% Salida:
%   C     : precio exacto de la opción lookback CALL
%
% Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = lookback_fixed_call_formulae(S0,K,r,q,sigma,T)

b   = r - q;  
dt  = sigma*sqrt(T);
d   = (log(S0./K) + (b + 0.5*sigma^2)*T) / dt;
d2  = d - dt;
db  = d - 2*b/sigma*sqrt(T);

if abs(b) < 1e-10
    % Límite b -> 0 (expansión de primer orden)
    phi = @(x) exp(-0.5*x.^2)/sqrt(2*pi);
    extra = 0.5*sigma^2*T * S0.*exp(-q*T).*phi(d);
    C = S0.*exp(-q*T).*normcdf(d) - K.*exp(-r*T).*normcdf(d2) + extra;
else
    C = S0.*exp(-q*T).*normcdf(d) ...
      - K .*exp(-r*T).*normcdf(d2) ...
      + (sigma^2/(2*b))*S0 .* ( exp(-q*T).*normcdf(d) ...
      - exp(-r*T).*(S0./K).^-((2*b)/(sigma^2)) .* normcdf(db) );
end

% Ajuste para S0 > K considerando el histórico inicial M=S0
idx = (S0 > K);
if any(idx)
    if abs(b) < 1e-10
        extraS = 0.5*sigma^2*T * S0(idx).*exp(-q*T).*0.3989422804014337; % phi(0)
        HS = S0(idx).*exp(-q*T).*0.5 - S0(idx).*exp(-r*T).*0.5 + extraS; 
    else
        dS  = (0 + (b + 0.5*sigma^2)*T) / dt;
        d2S = dS - dt;
        dbS = dS - 2*b/sigma*sqrt(T);
        HS = S0(idx).*exp(-q*T).*normcdf(dS) ...
           - S0(idx).*exp(-r*T).*normcdf(d2S) ...
           + (sigma^2/(2*b))*S0(idx).*( exp(-q*T).*normcdf(dS) ...
           - exp(-r*T).*normcdf(dbS) );
    end
    C(idx) = exp(-r*T).*(S0(idx)-K) + HS;
end

end
