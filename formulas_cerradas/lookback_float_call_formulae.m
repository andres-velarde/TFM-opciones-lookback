%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fórmula exacta Lookback CALL flotante
%   ----------------------------------------------------------------------
%   Calcula el precio de una opción Lookback CALL flotante usando la 
%   fórmula analítica para monitorización continua.
%
%   Entrada:
%       S_eval  - Valor del activo subyacente al inicio
%       r       - Tasa libre de riesgo (constante)
%       q       - Tasa de dividendos (constante)
%       sigma   - Volatilidad (constante)
%       T       - Tiempo hasta madurez
%
%   Salida:
%       price   - Precio exacto de la opción Lookback CALL flotante
%
%   Nota:
%       - Se evalúa en t=0 con S=J=S_eval.
%       - Fórmula para parámetros constantes.
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function price = lookback_float_call_formulae(S_eval, r, q, sigma, T)

    % Constantes auxiliares
    kappa = 2*(r - q) / sigma^2;
    tau   = T;

    % Variables d1, d2, d3 según la fórmula de lookback
    d1 = ((r - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau));
    d2 = d1 - sigma * sqrt(tau);
    d3 = d1 - kappa * sigma * sqrt(tau);

    % Evaluado en S = J (S_eval)
    term1 =  S_eval * exp(-q*tau) * normcdf(d1);
    term2 = -S_eval * exp(-r*tau) * normcdf(d2);
    term3 = (S_eval / kappa) * exp(-r*tau) * (normcdf(-d3) - exp((r - q)*tau) * normcdf(-d1));

    % Precio final
    price = term1 + term2 + term3;

end