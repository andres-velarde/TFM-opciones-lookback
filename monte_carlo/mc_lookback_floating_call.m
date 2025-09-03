%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Monte Carlo Lookback CALL flotante
%   ----------------------------------------------------------------------
%   Calcula el precio de una opción Lookback CALL flotante usando simulación
%   Monte Carlo bajo el modelo GBM con coeficientes dependientes del tiempo.
%
%   Payoff: S_T - m_T >= 0 (mínimo histórico m_T)
%
%   Entrada:
%       S0      - Precio inicial del activo subyacente
%       T       - Tiempo hasta madurez
%       nSteps  - Número de pasos discretos para monitorización
%       nPaths  - Número de simulaciones Monte Carlo
%       rfun    - Tasa libre de riesgo
%       qfun    - Tasa de dividendos
%       sigfun  - Volatilidad
%
%   Salida:
%       price   - Precio estimado de la opción Lookback CALL flotante
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function price = mc_lookback_floating_call(S0,T,nSteps,nPaths,rfun,qfun,sigfun)

    % Paso temporal
    dt  = T / (nSteps + 1);
    t   = linspace(0, T, nSteps+2); tL = t(1:end-1);
    r   = rfun(tL);  
    q = qfun(tL);  
    sg = sigfun(tL);

    % Coeficientes de movimiento geométrico
    mu  = (r - q - 0.5*sg.^2) * dt;                      % media logarítmica
    vol = sg * sqrt(dt);                                 % desviación estándar
    DF  = exp(-sum(r) * dt);                             % factor de descuento

    % Generar trayectorias normales
    Z = randn(nPaths, nSteps);                          % N(0,1) i.i.d.
    S = S0 * ones(nPaths,1);
    m = S;                                              % mínimo acumulado

    % Simulación paso a paso
    for j = 1:nSteps
        S = S .* exp(mu(j) + vol(j) .* Z(:,j));
        m = min(m, S);                                  % actualizar mínimo
    end

    % Payoff CALL flotante y descuento
    pay = S - m;                                        % siempre >=0
    price = DF * mean(pay);

end

