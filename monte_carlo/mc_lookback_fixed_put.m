%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo: Lookback PUT con strike fijo
% ------------------------------------------------------------------------
% Calcula el precio de una opción lookback tipo PUT con strike fijo K:
% Payoff: (K - mínimo(S_t))^+ con monitorización discreta en nSteps.
% Coeficientes de tasas y volatilidad pueden depender del tiempo.
%
% Entradas:
%   S0     - precio inicial del subyacente
%   K      - strike fijo
%   T      - tiempo hasta madurez
%   nSteps - número de pasos discretos en el tiempo
%   nPaths - número de simulaciones Monte Carlo
%   rfun   - tasa libre de riesgo en el tiempo
%   qfun   - dividendos en el tiempo
%   sigfun - volatilidad en el tiempo
%
% Salida:
%   price  - precio promedio descontado de la opción
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function price = mc_lookback_fixed_put(S0,K,T,nSteps,nPaths,rfun,qfun,sigfun)
    dt  = T/(nSteps+1);
    t   = linspace(0,T,nSteps+2); tL = t(1:end-1);
    r   = rfun(tL);  q = qfun(tL);  sg = sigfun(tL);

    mu  = (r - q - 0.5*sg.^2)*dt;
    vol = sg*sqrt(dt);
    DF  = exp(-sum(r)*dt);

    Z = randn(nPaths,nSteps);
    S = S0*ones(nPaths,1);
    m = S;                              % mínimo acumulado

    for j=1:nSteps
        S = S .* exp(mu(j) + vol(j).*Z(:,j));
        m = min(m,S);
    end

    pay = max(K - m, 0);
    price = DF * mean(pay);
end

