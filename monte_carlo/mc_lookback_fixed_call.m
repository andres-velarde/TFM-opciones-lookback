%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo: Lookback CALL de strike fijo
% ------------------------------------------------------------------------
% Simula paths del subyacente bajo GBM con r(t), q(t), sigma(t) dependientes
% del tiempo usando pasos discretos (nSteps) y nPaths simulaciones.
%
% Payoff:
%   max(M_T - K, 0), donde M_T es el máximo del subyacente a lo largo del
%   camino simulado.
%
% Entradas:
%   S0     : precio inicial del subyacente
%   K      : strike fijo
%   T      : madurez
%   nSteps : número de pasos discretos
%   nPaths : número de caminos Monte Carlo
%   rfun   : handle para la tasa libre de riesgo r(t)
%   qfun   : handle para la tasa de dividendos q(t)
%   sigfun : handle para la volatilidad sigma(t)
%
% Salida:
%   price  : precio estimado de la opción
%
% Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function price = mc_lookback_fixed_call(S0,K,T,nSteps,nPaths,rfun,qfun,sigfun)

    dt  = T / nSteps;
    tL  = (0:nSteps-1) * dt;          % instantes de evaluación por paso

    r   = rfun(tL);    q = qfun(tL);  sg = sigfun(tL);
    r   = r(:)';       q = q(:)';     sg = sg(:)';

    mu  = (r - q - 0.5*sg.^2) * dt;
    vol = sg * sqrt(dt);
    DF  = exp(-sum(r) * dt);

    Z = randn(nPaths, nSteps);          % N(0,1) i.i.d.
    S = S0 * ones(nPaths,1);
    M = S;                              % máximo acumulado

    for j = 1:nSteps
        S = S .* exp(mu(j) + vol(j) .* Z(:,j));
        M = max(M, S);
    end

    pay   = max(M - K, 0);              % payoff lookback CALL
    price = DF * mean(pay);
end

