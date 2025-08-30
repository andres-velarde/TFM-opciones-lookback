%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solución Crank-Nicolson para PUT lookback de strike fijo
% ------------------------------------------------------------------------
% Calcula L(t,S,m) para una opción PUT lookback con strike fijo, considerando
% la monitorización continua y usando Crank-Nicolson en S para cada nivel m.
%
% Dominio:
%   0 <= m <= S <= Smax
% Fronteras:
%   t = T: L(T,S,m) = max(K - m, 0)
%   m = 0: L(t,S,0) = K * exp(-∫_t^T r) (Dirichlet)
%   S = Smax: L(t,Smax,m) = max(K - m,0) * exp(-∫_t^T r) (Dirichlet)
%   S = m (diagonal): ∂_m L = 0 → diferencias hacia atrás (2º orden si k>=3)
%
% Entradas:
%   T, K, N, M - vencimiento, strike y mallas en tiempo y espacio
%   r, q, sigma - funciones de tiempo o valores constantes
%   Smax        - valor máximo de S considerado
%
% Salida:
%   L      - matriz 3D L(j,i,k) con t,S,m
%   t_vec  - vector de tiempos
%   S_vec  - vector de precios S
%   m_vec  - vector de mínimos m
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, t_vec, S_vec, m_vec] = fixed_put_cn(T, K, N, M, r, q, sigma, Smax)

    % Mallado
    dt    = T/(N+1);
    t_vec = linspace(0, T, N+2);
    S_vec = linspace(0, Smax, M+2);
    m_vec = S_vec;
    h     = S_vec(2) - S_vec(1);

    % Coeficientes temporales
    r_t = arrayfun(r,     t_vec);
    q_t = arrayfun(q,     t_vec);
    s_t = arrayfun(sigma, t_vec);

    % Descuento e^{-∫ r}
    psi = -arrayfun(@(t) integral(@(u) r(u), t, T, 'ArrayValued', true), t_vec);

    % Almacén
    L = nan(N+2, M+2, M+2);

    % Terminal t = T
    for k = 1:M+2
        L(end, k:M+2, k) = max(K - m_vec(k), 0);
    end

    % Fronteras Dirichlet en m=0 y S=Smax
    for j = 1:N+2
        L(j, :,   1) = K * exp(psi(j));
        L(j, end, :) = max(K - m_vec, 0) * exp(psi(j));
    end

    % Barrido en m = 2..M
    for k = 2:M
        iL = k+1; iR = M+1;
        Si = S_vec(iL:iR)';

        % Coefs "C" iniciales
        sigC = s_t(end); rC = r_t(end); qC = q_t(end);
        alphaC = 0.5*dt*( (sigC^2)*(Si.^2)/h^2 - ((rC - qC).*Si)/h );
        betaC  =      dt*( (sigC^2)*(Si.^2)/h^2 +  rC );
        gammaC = 0.5*dt*( (sigC^2)*(Si.^2)/h^2 + ((rC - qC).*Si)/h );

        % Retroceso temporal
        for j = N+2:-1:2
            
            % 2º orden hacia atrás en J si k<=M; 1º orden si k=M+1
            if k >=3
                L(j-1, k, k) = (4*L(j-1, k, k-1) - L(j-1, k, k-2))/3;
            else
                L(j-1, k, k) = L(j-1, k, k-1);
            end

            % RHS en j
            Vij   = reshape(L(j, iL:iR, k), [], 1);
            Vij_l = [ L(j, k, k); Vij(1:end-1) ];
            Vij_r = [ Vij(2:end); L(j, end, k) ];
            d = (2 - betaC).*Vij + alphaC.*Vij_l + gammaC.*Vij_r;

            % Coefs "P" (nivel j-1)
            sigP = s_t(j-1); rP = r_t(j-1); qP = q_t(j-1);
            alphaP = 0.5*dt*( (sigP^2)*(Si.^2)/h^2 - ((rP - qP).*Si)/h );
            betaP  =      dt*( (sigP^2)*(Si.^2)/h^2 +  rP );
            gammaP = 0.5*dt*( (sigP^2)*(Si.^2)/h^2 + ((rP - qP).*Si)/h );

            main = (2 + betaP);
            sub  = -alphaP(2:end);
            sup  = -gammaP(1:end-1);
            A = diag(main) + diag(sub,-1) + diag(sup,1);

            % Bordes nivel nuevo j-1
            d(1)   = d(1)   + alphaP(1)   * L(j-1, k, k);
            d(end) = d(end) + gammaP(end) * L(j-1, end, k);

            % Resolver sistema tridiagonal
            x = tridiagonal_matrix(A, d);
            L(j-1, iL:iR, k) = reshape(x, 1, []);

            % Reciclar coefs
            alphaC = alphaP; betaC = betaP; gammaC = gammaP;
        end
    end

    % Consistencia último plano
    L(:, M+1, M+1) = L(:, M+1, M);
end
