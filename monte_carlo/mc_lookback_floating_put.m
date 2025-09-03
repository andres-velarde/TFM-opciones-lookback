%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   mc_lookback_floating_put
%   ---------------------------------------------------------------
%   Monte Carlo para el precio de una opción lookback PUT flotante 
%   bajo un modelo GBM con coeficientes dependientes del tiempo.
%
%   Características:
%       • Opción PUT flotante europea
%       • Monitorización discreta en nSteps
%       • Generador normal estándar (antitéticos opcionales, no usados)
%       • Payoff:  M_T - S_T   (siempre >= 0)
%
%   Parámetros de entrada:
%       S0      - Precio inicial del subyacente
%       T       - Horizonte temporal (madurez)
%       nSteps  - Número de pasos de tiempo
%       nPaths  - Número de trayectorias simuladas
%       rfun    - Función r(t): tipo de interés (vectorizado)
%       qfun    - Función q(t): dividendos (vectorizado)
%       sigfun  - Función sigma(t): volatilidad (vectorizado)
%
%   Salida:
%       price   - Estimación Monte Carlo del precio de la opción
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function price = mc_lookback_floating_put(S0,T,nSteps,nPaths,rfun,qfun,sigfun)

    % Paso temporal y mallado
    tau = T / (nSteps+1);
    t   = linspace(0,T,nSteps+2);
    tL  = t(1:end-1);                 % nodos a izquierda
    
    % Coeficientes dependientes del tiempo
    r  = rfun(tL);  
    q  = qfun(tL);  
    sg = sigfun(tL);

    mu  = (r - q - 0.5*sg.^2) * tau;  % drift discreto
    vol = sg * sqrt(tau);             % volatilidad discreta
    DF  = exp(-sum(r)*tau);           % factor de descuento (rectángulo izqda)

    % Muestras normales
    Z = randn(nPaths, nSteps);

    % Inicialización
    S = S0 * ones(nPaths,1);          % subyacente
    M = S;                            % máximo alcanzado

    % Evolución
    for j = 1:nSteps
        S = S .* exp(mu(j) + vol(j)*Z(:,j));
        M = max(M, S);
    end

    % Payoff y descuento
    pay     = M - S;                  % PUT flotante
    discPay = DF * pay;

    % Valor esperado
    price = mean(discPay);
end
