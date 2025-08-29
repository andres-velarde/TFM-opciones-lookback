%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solución PDE para Lookback CALL de strike fijo con Crank–Nicolson
%   ----------------------------------------------------------------------
%   Dominio: 0 <= S <= J <= Smax
%   Variables de estado:
%       • S: precio instantáneo del subyacente
%       • J: máximo alcanzado por S hasta el tiempo t
%
%   EDP por plano J fijo:
%       L_t + (r-q) S L_S + 0.5 sigma^2 S^2 L_SS - r L = 0,   0 < S < J
%
%   Condiciones de contorno:
%       • Terminal en t = T:
%            L(T,S,J) = max(J - K, 0)
%       • S = 0:
%            L(t,0,J) = max(J - K, 0) * exp(-∫_t^T r)
%       • J = Smax (cota superior):
%            L(t,S,Smax) = max(Smax - K, 0) * exp(-∫_t^T r)
%       • S = J (diagonal):
%            ∂L/∂J = 0  →  se impone con esquema unilateral de 2º orden
%
%   Entrada:
%       T       - Tiempo hasta madurez
%       K       - Strike fijo
%       N, M    - Número de pasos en tiempo y espacio
%       r, q    - Funciones de tasa libre de riesgo y dividendos (dependen de t)
%       sigma   - Función de volatilidad (dependiente de t)
%       Smax    - Cota superior del dominio para J (≈ múltiplo de max(K,S0))
%
%   Salida:
%       L       - Tensor 3D de la solución L(j,i,k) en (t,S,J)
%       t_vec   - Vector de tiempos discretizados
%       S_vec   - Vector de nodos espaciales para S
%       J_vec   - Vector de nodos espaciales para J
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, t_vec, S_vec, J_vec] = fixed_call_cn(T, K, N, M, r, q, sigma, Smax)

    % === Mallado temporal y espacial ===
    dt    = T/(N+1);
    t_vec = linspace(0, T, N+2);     % discretización temporal
    S_vec = linspace(0, Smax, M+2);  % discretización espacial en S
    J_vec = S_vec;                   % reutilizamos misma malla para J
    h     = S_vec(2) - S_vec(1);     % paso espacial

    % === Coeficientes temporales en nodos ===
    r_t = arrayfun(r,     t_vec);
    q_t = arrayfun(q,     t_vec);
    s_t = arrayfun(sigma, t_vec);

    % Descuento exacto e^{-∫ r}
    psi = -arrayfun(@(t) integral(@(u) r(u), t, T, 'ArrayValued', true), t_vec);

    % === Inicialización del tensor solución ===
    L = nan(N+2, M+2, M+2);          % (j,i,k) = (tiempo, S, J)

    % === Condición terminal en t = T ===
    for k = 1:M+2
        L(end, 1:k, k) = max(J_vec(k) - K, 0);
    end

    % === Contornos Dirichlet globales ===
    for j = 1:N+2
        L(j, 1,   :) = max(J_vec - K, 0) * exp(psi(j));      % S=0
        L(j, :, end) = max(J_vec(end) - K, 0) * exp(psi(j)); % J=Smax
    end

    % === Barrido descendente en J (desde el borde hacia la diagonal) ===
    for k = M+1:-1:3

        % Interior en S dentro del plano J_k
        iL = 2; iR = k-1;   
        if iR < iL, continue; end
        Si = S_vec(iL:iR)';               % columna de nodos interiores

        % Coeficientes iniciales (nivel C, al empezar en t=T)
        sigC = s_t(end); rC = r_t(end); qC = q_t(end);
        alphaC = 0.5*dt*( (sigC^2)*(Si.^2)/h^2 - ((rC - qC).*Si)/h );
        betaC  =      dt*( (sigC^2)*(Si.^2)/h^2 +  rC );
        gammaC = 0.5*dt*( (sigC^2)*(Si.^2)/h^2 + ((rC - qC).*Si)/h );

        % Retroceso en el tiempo
        for j = N+2:-1:2

            % Neumann en la diagonal S=J: imponemos ∂L/∂J=0 en nivel actual
            L(j-1, k, k) = L(j-1, k, k+1);

            % === Construcción del RHS (nivel C) ===
            Vij   = reshape(L(j, iL:iR, k), [], 1);
            Vij_l = [ L(j, 1, k); Vij(1:end-1) ];   % vecino izq (S=0)
            Vij_r = [ Vij(2:end); L(j, k, k)   ];   % vecino der (S=J diagonal)

            d = (2 - betaC).*Vij + alphaC.*Vij_l + gammaC.*Vij_r;

            % === Coeficientes en nivel P (j-1) ===
            sigP = s_t(j-1); rP = r_t(j-1); qP = q_t(j-1);
            alphaP = 0.5*dt*( (sigP^2)*(Si.^2)/h^2 - ((rP - qP).*Si)/h );
            betaP  =      dt*( (sigP^2)*(Si.^2)/h^2 +  rP );
            gammaP = 0.5*dt*( (sigP^2)*(Si.^2)/h^2 + ((rP - qP).*Si)/h );

            % Matriz tridiagonal del LHS
            main = (2 + betaP);
            sub  = -alphaP(2:end);
            sup  = -gammaP(1:end-1);
            A = diag(main) + diag(sub, -1) + diag(sup, 1);

            % Ajuste de contornos en RHS:
            d(1)   = d(1) + alphaP(1)   * L(j-1, 1, k); % borde S=0
            d(end) = d(end) + gammaP(end)* L(j-1, k, k); % borde S=J diagonal

            % Resolver sistema tridiagonal
            L(j-1, iL:iR, k) = tridiagonal_matrix(A, d);

            % Reimponer condición Neumann en j-1
            L(j-1, k, k) = L(j-1, k, k+1);

            % Reciclar coeficientes para siguiente paso
            alphaC = alphaP; betaC = betaP; gammaC = gammaP;
        end
    end

    % === Consistencia en el último plano (k=2, solo un interior) ===
    L(:, 2, 2) = L(:, 2, 3);

end
