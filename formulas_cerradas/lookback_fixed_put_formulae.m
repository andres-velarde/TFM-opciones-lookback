function P = lookback_fixed_put_formulae(S0,K,r,q,sigma,T)

    b   = r - q;
    tau = T;                  % tiempo hasta madurez
    dt  = sigma*sqrt(tau);

    % Parámetro lambda
    lam = 2*b/(sigma^2);

    % Caso 1: K < S0  (mínimo inicial M=S0 > K)
    if K < S0
        X = K;
        d1 = (log(S0./X) + (b + 0.5*sigma^2)*tau) / dt;
        d2 = d1 - dt;

        term1 = K*exp(-r*tau).*normcdf(-d2) ...
              - S0*exp(-q*tau).*normcdf(-d1);

        term2 = (sigma^2/(2*b))*S0 .* exp(-q*tau) .* ...
                ( (S0./X).^(-lam) .* normcdf(-d1 + lam*dt) ...
                 - exp(b*tau).*normcdf(-d1) );

        P = term1 + term2;

    % Caso 2: K >= S0  (mínimo inicial M=S0 <= K)
    else
        M = S0;
        X = M;
        d1 = (log(S0./X) + (b + 0.5*sigma^2)*tau) / dt;
        d2 = d1 - dt;

        term1 = (K - M)*exp(-r*tau) ...
              - S0*exp(-q*tau).*normcdf(-d1) ...
              + M*exp(-r*tau).*normcdf(-d2);

        term2 = (sigma^2/(2*b))*S0 .* exp(-q*tau) .* ...
                ( (S0./X).^(-lam) .* normcdf(-d1 + lam*dt) ...
                 - exp(b*tau).*normcdf(-d1) );

        P = term1 + term2;
    end
end
