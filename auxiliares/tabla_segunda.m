%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generar tabla de errores CN vs fórmula para opciones Lookback Fixed Put
% Errores en función de sigma y S0 ∈ [60, 96]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Parámetros comunes ---
r_val   = 0.10;
q_val   = 0.25;
T       = 10/12;
K  = 100;

rfun = @(t) r_val + t-t;
qfun = @(t) q_val + t-t;

N = 100;   % pasos temporales
M = 100;   % pasos espaciales

% --- Rejilla de evaluación ---
sigmas  = [5, 12, 19, 26, 33, 40]/100;   % % → en decimales
S0_vals = linspace(10, 90, 9);
nS_fix = length(S0_vals);

% --- Almacén de resultados ---
CN = zeros(length(sigmas), length(S0_vals));
MC = zeros(length(sigmas), length(S0_vals));

% --- Bucle principal ---
for i = 1:length(sigmas)
    sigma = sigmas(i);
    sigmafun = @(t) sigma + t-t;
    
    for j = 1:length(S0_vals)
        S0 = S0_vals(j);

        % Valor CN
        CN(i, j) = lookback_fixed_put(T,N,M,@(t)r_val+t-t,@(t)q_val+t-t,@(t)sigma+t-t,K,S0,max(5, 2*S0) * (1+sigma));

        % Valor exacto
        MC(i, j) = mc_lookback_fixed_put(S0,K,T,2e4,2e4, ...
        %             @(t)r_val+t-t, @(t)q_val+t-t, @(t)sigma+t-t);
        
    end
end

errores = (CN-MC) ./ MC*100;
% --- Mostrar resultados ---
disp('=== Errores relativos (%) ===');
disp(array2table(errores, 'VariableNames', ...
    strcat("S0_", string(S0_vals)), 'RowNames', strcat("sigma_", string(sigmas*100))));
