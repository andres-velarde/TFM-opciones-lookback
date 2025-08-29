%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solución PDE para Lookback CALL flotante en variables de similaridad
%   ----------------------------------------------------------------------
%   Dominio: xi ∈ [1, xi_max]
%   Condiciones de frontera:
%       • Dirichlet en xi = xi_max
%       • Condición unilateral de 2º orden en xi = 1:
%         W'(1) ≈ (4 W_1 - W_2)/(2h + 3)
%
%   Entrada:
%       T       - Tiempo hasta madurez
%       N, M    - Número de pasos temporales y espaciales
%       r, q    - Funciones de tasa libre de riesgo y dividendos (dependientes de t)
%       sigma   - Función de volatilidad (dependiente de t)
%       xi_max  - Valor máximo en xi
%
%   Salida:
%       W        - Matriz de solución W(t_j, xi_i)
%       vector_t - Vector de tiempos t_j
%       vector_xi- Vector de nodos xi_i
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, vector_t, vector_xi] = floating_call_similaridad(T, N, M, r, q, sigma, xi_max)

    % Mallado temporal y espacial
    dt        = T / (N + 1);
    vector_t  = linspace(0, T, N + 2);        % j=1..N+2 (t_j)
    vector_xi = linspace(1, xi_max, M + 2);   % i=0..M+1 (xi_i)
    h         = (xi_max - 1) / (M + 1);
    xi_i      = vector_xi(2:M+1)';            % nodos interiores

    % Inicialización de la solución
    W = zeros(N + 2, M + 2);

    % Coeficientes temporales en nodos
    sig_t = arrayfun(sigma, vector_t);
    r_t   = arrayfun(r,     vector_t);
    q_t   = arrayfun(q,     vector_t);

    % Condición terminal: W(T, xi) = xi - 1
    W(end, :) = vector_xi - 1;

    % Condición Dirichlet en xi_max
    psi = -arrayfun(@(t) integral(@(s) r(s), t, T), vector_t);
    mu  = -arrayfun(@(t) integral(@(s) q(s), t, T), vector_t);
    W(:, end) = vector_xi(end)*exp(mu) - exp(psi);

    % Coeficientes iniciales para la primera iteración temporal
    sigC = sig_t(end); rc = r_t(end); qc = q_t(end);
    alphaC = 0.5*dt * ( (sigC^2)*(xi_i.^2)/h^2 - ((rc - qc).*xi_i)/h );
    betaC  =      dt * ( (sigC^2)*(xi_i.^2)/h^2 +  rc );
    gammaC = 0.5*dt * ( (sigC^2)*(xi_i.^2)/h^2 + ((rc - qc).*xi_i)/h );

    % Marcha hacia atrás en el tiempo
    for j = N+2 : -1 : 2

        % Nodos interiores en el tiempo j
        Wj = W(j, 2:M+1)';  % columna, tamaño M

        % ---------- RHS (nivel C) ----------
        % Vecino izquierdo con Neumann en xi=1
        W0_C  = (4*Wj(1) - Wj(2)) / (2*h + 3);
        Vleft = [ W0_C ; Wj(1:end-1) ];

        % Vecino derecho: último nodo Dirichlet
        Vright = [ Wj(2:end) ; W(j, end) ];

        d = (2 - betaC).*Wj + alphaC.*Vleft + gammaC.*Vright;

        % ---------- LHS (nivel P) ----------
        sigP = sig_t(j-1); rp = r_t(j-1); qp = q_t(j-1);
        alphaP = 0.5*dt * ( (sigP^2)*(xi_i.^2)/h^2 - ((rp - qp).*xi_i)/h );
        betaP  =      dt * ( (sigP^2)*(xi_i.^2)/h^2 +  rp );
        gammaP = 0.5*dt * ( (sigP^2)*(xi_i.^2)/h^2 + ((rp - qp).*xi_i)/h );

        % Matriz tridiagonal
        A = diag(2 + betaP) - diag(alphaP(2:end), -1) - diag(gammaP(1:end-1), 1);

        % Ajuste Neumann en el primer nodo
        A(1,1) = A(1,1) - (4/(2*h + 3)) * alphaP(1);
        A(1,2) = A(1,2) + (1/(2*h + 3)) * alphaP(1);

        % Dirichlet en xi_max pasa a RHS
        d(end) = d(end) + gammaP(end) * W(j-1, end);

        % Resolver el sistema
        W(j-1, 2:M+1) = tridiagonal_matrix(A, d);

        % Reconstrucción del nodo xi=1 en el nivel nuevo
        W(j-1, 1) = (4*W(j-1, 2) - W(j-1, 3)) / (2*h + 3);

        % Actualizar coeficientes para siguiente paso
        alphaC = alphaP; betaC = betaP; gammaC = gammaP;
    end
end
