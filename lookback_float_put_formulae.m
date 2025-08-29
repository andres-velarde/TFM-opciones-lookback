%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   lookback_float_put_formulae
%   ---------------------------------------------------------------
%   Fórmula exacta (monitorización continua) para la opción lookback PUT 
%   flotante bajo supuestos de coeficientes constantes.
%
%   Características:
%       • Evaluación en t=0 con S = J = S_eval
%       • Modelo GBM con r, q y sigma constantes
%       • Payoff:  M_T - S_T
%
%   Parámetros de entrada:
%       S_eval  - Precio inicial del subyacente (S0)
%       r       - Tipo de interés libre de riesgo
%       q       - Tasa de dividendos
%       sigma   - Volatilidad
%       T       - Horizonte temporal (madurez)
%
%   Salida:
%       price   - Precio exacto de la opción lookback PUT flotante
%
%   Referencia:
%       Paul Willmot, Derivatives: The Theory and Practice of Financial
%       Engineering, 1998
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function price = lookback_float_put_formulae(S_eval, r, q, sigma, T)
    % Parámetros auxiliares
    kappa = 2*(r - q) / sigma^2;
    tau   = T;

    % Variables d
    d1 = ((r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau));
    d2 = d1 - sigma * sqrt(tau);
    d3 = d1 - kappa * sigma * sqrt(tau);

    % Componentes de la fórmula
    term1 = -S_eval * exp(-q*tau) * normcdf(-d1);
    term2 =  S_eval * exp(-r*tau) * normcdf(-d2);
    term3 = (S_eval / kappa) * exp(-r*tau) * ...
            (-normcdf(d3) + exp((r - q)*tau) * normcdf(d1));

    % Precio final
    price = term1 + term2 + term3;
end
