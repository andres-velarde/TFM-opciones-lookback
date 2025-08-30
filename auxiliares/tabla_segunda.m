%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generar tabla de errores CN vs fórmula para opciones Lookback Fixed Put
% Errores en función de sigma y S0 ∈ [60, 96]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Parámetros comunes ---
r  = 0.1;
q  = 0.25;
T  = 10/12;
K  = 80;

rfun = @(t) r + 0*t;
qfun = @(t) q + 0*t;

N = 200;   % pasos temporales
M = 200;   % pasos espaciales

% --- Rejilla de evaluación ---
sigmas  = [10, 14, 18, 22, 26, 30]/100;   % % → en decimales
S0_vals = [60 66 72 78 84 90 96];

% --- Almacén de resultados ---
Errores = zeros(length(sigmas), length(S0_vals));

% --- Bucle principal ---
for i = 1:length(sigmas)
    sigma = sigmas(i);
    sigmafun = @(t) sigma + 0*t;
    
    for j = 1:length(S0_vals)
        S0 = S0_vals(j);

        % Valor CN
        CN = lookback_fixed_put(T,N,M,rfun,qfun,sigmafun,K,S0,S0*2.75);

        % Valor exacto
        Exacto = lookback_fixed_put_formulae(S0,K,r,q,sigma,T);

        % Error relativo (%)
        Errores(i,j) = 100 * (CN - Exacto) / Exacto;
    end
end

% --- Mostrar resultados ---
disp('=== Errores relativos (%) ===');
disp(array2table(Errores, 'VariableNames', ...
    strcat("S0_", string(S0_vals)), 'RowNames', strcat("sigma_", string(sigmas*100))));
