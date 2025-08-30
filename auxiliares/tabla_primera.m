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
Smax  = 500;
spots_float = 60:20:140;   % rango para flotante
spots_fixed = 100:10:140;  % rango para fija
sigmas = [0.10 0.14 0.18 0.22 0.26 0.30];  % volatilidades

% --- Configuraciones de malla ---
grids = [10 10; 100 100; 200 200];

% --- Resultados ---
for g = 1:size(grids,1)
    N = grids(g,1);
    M = grids(g,2);

    fprintf('== Grid N=%d, M=%d ==\n', N, M);

    for s = 1:length(sigmas)
        sig_val = sigmas(s);

        % --- FLOATING PUT ---
        tic;
        errs_float = zeros(size(spots_float));
        for idx = 1:length(spots_float)
            S0 = spots_float(idx);
            % Precio CN
            Vcn = lookback_float_put(T,N,M,@(t)r_val+t-t,@(t)q_val+t-t,@(t)sig_val+t-t,S0);
            % Fórmula cerrada
            Vex = lookback_float_put_formulae(S0,r_val,q_val,sig_val,T);
            errs_float(idx) = abs(Vcn - Vex)/Vex*100;
        end
        t_float = toc/length(spots_float);
        err_float = mean(errs_float);

        % --- FIXED PUT ---
        tic;
        errs_fixed = zeros(size(spots_fixed));
        for idx = 1:length(spots_fixed)
            S0 = spots_fixed(idx);
            % Precio CN
            Vcn = lookback_fixed_put(T,N,M,@(t)r_val+t-t,@(t)q_val+t-t,@(t)sig_val+t-t,K,S0,max(200, S0*2));
            % Fórmula cerrada
            Vex = lookback_fixed_put_formulae(S0,K,r_val,q_val,sig_val,T);
            errs_fixed(idx) = abs(Vcn - Vex)/Vex*100;
        end
        t_fixed = toc/length(spots_fixed);
        err_fixed = mean(errs_fixed);

        % --- Mostrar fila ---
        fprintf('sigma=%.2f | Float Err=%.2f%%, T=%.2f s | Fixed Err=%.2f%%, T=%.2f s\n',...
            sig_val*100, err_float, t_float, err_fixed, t_fixed);
    end
    fprintf('\n');
end
