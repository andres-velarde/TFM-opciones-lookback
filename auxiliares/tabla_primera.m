%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errores y tiempos: opciones lookback flotantes y fijas (CN vs fórmula)
% -------------------------------------------------------------------------
% Calcula errores relativos (%) y tiempos de ejecución promedio
% para diferentes tamaños de malla (N,M) y volatilidades sigma.
%
% Autor: Andrés Velarde Náñez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% --- Parámetros comunes ---
r_val = 0.1;
q_val = 0.25;
T     = 10/12;
K     = 100;
Smax  = 300;
spots_float = 50:20:150;   % rango para flotante
spots_fixed = 100:10:150;  % rango para fija
sigmas = [0.05 0.12 0.19 0.26 0.33 0.40];  % volatilidades

% --- Configuraciones de malla ---
grids = [10 10; 100 100; 250 250];

% --- Resultados ---
for g = 1:size(grids,1)
    N = grids(g,1);
    M = grids(g,2);

    fprintf('== Grid N=%d, M=%d ==\n', N, M);

    for s = 1:length(sigmas)
        sig_val = sigmas(s);

        % --- FLOATING PUT ---
        tic;
        Vcn = zeros(size(spots_float));
        Vex = zeros(size(spots_float));
        for idx = 1:length(spots_float)
            S0 = spots_float(idx);
            % Precio CN
            Vcn(idx) = lookback_float_put(T,N,M,@(t)r_val+t-t,@(t)q_val+t-t,@(t)sig_val+t-t,S0);
            % Fórmula cerrada
            Vex(idx) = lookback_float_put_formulae(S0,r_val,q_val,sig_val,T);
        end
        t_float = toc/length(spots_float);
        err_float = max(abs(Vcn-Vex)) / max(Vex)*100;

        % --- FIXED PUT ---
        tic;
        errs_fixed = zeros(size(spots_fixed));
        Vcn = zeros(size(spots_float));
        Vex = zeros(size(spots_float));
        for idx = 1:length(spots_fixed)
            S0 = spots_fixed(idx);
            % Precio CN
            Vcn(idx) = lookback_fixed_put(T,N,M,@(t)r_val+t-t,@(t)q_val+t-t,@(t)sig_val+t-t,K,S0,300);
            % Fórmula cerrada
            Vex(idx) = lookback_fixed_put_formulae(S0,K,r_val,q_val,sig_val,T);
        end
        t_fixed = toc/length(spots_fixed);
        err_fixed = max(abs(Vcn-Vex)) / max(Vex)*100;

        % --- Mostrar fila ---
        fprintf('sigma=%.2f | Float Err=%.2f%%, T=%.2f s | Fixed Err=%.2f%%, T=%.2f s\n',...
            sig_val*100, err_float, t_float, err_fixed, t_fixed);
    end
    fprintf('\n');
end
