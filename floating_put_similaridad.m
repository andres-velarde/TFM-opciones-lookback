%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   floating_put_similaridad.m
%   ---------------------------------------------------------------
%   Resuelve la EDP asociada a una lookback PUT flotante mediante 
%   un esquema Crank–Nicolson en el dominio reducido (t, xi).
%
%   INPUT:
%       T      - Horizonte temporal
%       N      - Número de pasos en tiempo
%       M      - Número de pasos en espacio (en variable xi)
%       r(t)   - Función de tipo de interés
%       q(t)   - Función de dividendos
%       sigma  - Función de volatilidad
%
%   OUTPUT:
%       W         - Matriz solución en nodos reducidos (tiempo x espacio)
%       vector_t  - Mallado temporal
%       vector_xi - Mallado espacial reducido
%
%   Condiciones de contorno:
%       - Dirichlet en xi = 0
%       - Neumann en xi = 1
%
%   Terminal condition: W(T,xi) = 1 - xi
%
%   Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, vector_t, vector_xi] = floating_put_similaridad(T, N, M, r, q, sigma)
% Resuelve W(t,xi) para una lookback PUT flotante
% Dominio: [0,T] x [0,1]; Neumann en xi=1, Dirichlet en xi=0.

    % Mallado
    dt        = T / (N + 1);
    vector_t  = linspace(0, T, N + 2);     % t_0,...,t_{N+1} (t_{N+1}=T)
    vector_xi = linspace(0, 1, M + 2);     % xi_0=0,...,xi_{M+1}=1
    h         = 1 / (M + 1);
    xi_i      = vector_xi(2:M+1)';         % nodos interiores (i=1..M)

    % Solución (tiempo x espacio)
    W = zeros(N + 2, M + 2);

    % Funciones temporales en nodos
    sig_t = arrayfun(sigma, vector_t);
    r_t   = arrayfun(r,     vector_t);
    q_t   = arrayfun(q,     vector_t);

    % Terminal: W(T,xi) = 1 - xi
    W(end, :) = 1 - vector_xi;

    % Dirichlet en xi=0: W(t,0) = e^{-∫_t^T r}
    int_r = -arrayfun(@(t) integral(@(s) r(s), t, T), vector_t);
    W(:, 1) = exp(int_r);

    % Coeficientes en nivel j (RHS), empezar en t=T
    sigC = sig_t(end); rc = r_t(end); qc = q_t(end);
    alphaC = 0.5 * dt * ( (sigC^2) * (xi_i.^2) / h^2 - (rc - qc) * xi_i / h );
    betaC  =       dt * ( (sigC^2) * (xi_i.^2) / h^2 + rc );
    gammaC = 0.5 * dt * ( (sigC^2) * (xi_i.^2) / h^2 + (rc - qc) * xi_i / h );

    % Marcha hacia atrás
    for j = N + 2 : -1 : 2
        
        Wj = W(j, 2:M+1)';         % interior i=1..M en t_j
        d  = (2 - betaC).*Wj ...
           + alphaC .* [W(j,1); Wj(1:end-1)] ...
           + gammaC .* [Wj(2:end); 0];
        
        d(1) = d(1) + alphaC(1) * W(j, 1);
            
        % Neumann en xi=1 (última entrada): sustituye fantasma en el RHS
        d(end) = (2 - betaC(end))*Wj(end) ...
       + alphaC(end)*Wj(end-1) ...
       + gammaC(end) * ( (-4*Wj(end) + Wj(end-1)) / (2*h - 3) );
      
        % Coeficientes en j-1 (LHS)
        sigP = sig_t(j-1); rp = r_t(j-1); qp = q_t(j-1);
        alphaP = 0.5 * dt * ( (sigP^2) * (xi_i.^2) / h^2 - (rp - qp) * xi_i / h );
        betaP  =       dt * ( (sigP^2) * (xi_i.^2) / h^2 + rp );
        gammaP = 0.5 * dt * ( (sigP^2) * (xi_i.^2) / h^2 + (rp - qp) * xi_i / h );

        % Matriz A y ajuste Neumann en la última fila
        A = diag(2 + betaP) - diag(alphaP(2:end), -1) - diag(gammaP(1:end-1), 1);
        A(end,end)   = A(end,end)   + (4/(2*h - 3)) * gammaP(end);
        A(end,end-1) = A(end,end-1) - (1/(2*h - 3)) * gammaP(end);

        % Dirichlet izq. en LHS: pasar αP(1)*W_0^{j-1} al RHS con signo -
        d(1) = d(1) + alphaP(1) * W(j-1, 1);

        % Resolver
        W(j-1, 2:M+1) = tridiagonal_matrix(A, d)';
        
        % Nodo en xi=1
        W(j-1, M+2) = (-4*W(j-1, M+1) + W(j-1, M)) / (2*h-3);
        
        % Actualizar coeficientes para siguiente paso
        alphaC = alphaP; betaC = betaP; gammaC = gammaP;

    end
end
