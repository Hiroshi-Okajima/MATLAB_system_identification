%% subspace_sysid_demo.m  (v3 — blog figure version)
%  Subspace Identification: N4SID, MOESP, and CVA
%
%  Generates 4 figures for the blog article:
%    Figure 1: Model Order Selection (bar chart of loss function)
%    Figure 2: Bode Plot — True vs N4SID vs ssest (N4SID+PEM)
%    Figure 3: Residual Analysis (autocorrelation + cross-correlation)
%    Figure 4: Pole-Zero Map — True vs Identified
%
%  Three-method comparison (N4SID/MOESP/CVA) is printed as a table
%  in the console — suitable for inclusion in the blog as text.
%
%  Blog: https://blog.control-theory.com/entry/subspace-identification
%  Hub:  https://blog.control-theory.com/entry/system-identification
%
%  Author: Hiroshi Okajima, Kumamoto University
%  Repository: https://github.com/Hiroshi-Okajima/MATLAB_system_identification

clear; close all; clc;

%% 1. Define the true plant (3rd-order, 1 input, 2 outputs)
Ts = 0.1;

A = diag([0.9, 0.7, 0.4]);
B = [1; 1; 1];
C = [1 1 0.5;
     0.5 -1 1];
D = [0; 0];

sys_true = ss(A, B, C, D, Ts);
n_true = size(A, 1);

fprintf('=== True Plant ===\n');
fprintf('Order: %d,  Inputs: %d,  Outputs: %d\n', ...
    n_true, size(B,2), size(C,1));
fprintf('Poles: %.1f, %.1f, %.1f\n\n', eig(A));

%% 2. Generate input-output data
N = 3000;
rng(123);
u = randn(N, 1);

y_clean = lsim(sys_true, u);
noise_std = 0.15;
e = noise_std * randn(N, size(C,1));
y = y_clean + e;

data = iddata(y, u, Ts);
data_est = data(1:2000);
data_val = data(2001:end);

fprintf('Data: %d samples (est: %d, val: %d)\n', N, 2000, 1000);
fprintf('Noise std: %.3f\n\n', noise_std);

%% 3. Figure 1: Model Order Selection
fprintf('=== Model Order Selection ===\n');

max_order = 10;
loss = zeros(max_order, 1);
for n_test = 1:max_order
    sys_test = n4sid(data_est, n_test);
    loss(n_test) = sys_test.Report.Fit.LossFcn;
end

% Find the elbow
rel_drop = -diff(loss) ./ abs(loss(1:end-1));
[~, n_recommended] = max(rel_drop);

figure('Name', 'Figure 1: Model Order Selection', ...
       'Position', [100, 100, 560, 380]);
colors = repmat([0.3 0.5 0.8], max_order, 1);
colors(n_recommended, :) = [0.9 0.3 0.2];
b = bar(1:max_order, loss);
b.FaceColor = 'flat';
b.CData = colors;
xlabel('Model Order n', 'FontSize', 13);
ylabel('Loss Function', 'FontSize', 13);
title('Model Order Selection via N4SID', 'FontSize', 14);
grid on;
set(gca, 'XTick', 1:max_order, 'FontSize', 11);
text(n_recommended, loss(n_recommended)*1.08, ...
    sprintf('  n = %d (recommended)', n_recommended), ...
    'FontSize', 12, 'Color', [0.9 0.3 0.2], 'FontWeight', 'bold');

fprintf('Recommended order: %d (true order = %d)\n\n', n_recommended, n_true);

%% 4. Identify with three methods + ssest (console table only)
n = n_true;

opt_n4sid = n4sidOptions('N4Weight', 'auto');
sys_n4sid = n4sid(data_est, n, opt_n4sid);

opt_moesp = n4sidOptions('N4Weight', 'MOESP');
sys_moesp = n4sid(data_est, n, opt_moesp);

opt_cva = n4sidOptions('N4Weight', 'CVA');
sys_cva = n4sid(data_est, n, opt_cva);

sys_ssest = ssest(data_est, n, 'Ts', Ts);

% Compute validation fit
u_val = data_val.u;
y_val = data_val.y;

calc_fit = @(y_hat, y_ref) max(0, ...
    (1 - norm(y_ref - y_hat) / norm(y_ref - mean(y_ref))) * 100);

methods = {'N4SID', 'MOESP', 'CVA', 'ssest'};
systems = {sys_n4sid, sys_moesp, sys_cva, sys_ssest};
n_methods = length(methods);
n_outputs = size(C, 1);
fit_table = zeros(n_outputs, n_methods);

for j = 1:n_methods
    y_sim = lsim(systems{j}, u_val);
    for i = 1:n_outputs
        fit_table(i, j) = calc_fit(y_sim(:,i), y_val(:,i));
    end
end

% Print comparison table (for blog article text)
fprintf('=== Validation Fit (%%) — for blog table ===\n');
fprintf('| Output | N4SID | MOESP | CVA   | ssest (N4SID+PEM) |\n');
fprintf('|--------|-------|-------|-------|-------------------|\n');
for i = 1:n_outputs
    fprintf('| %d      | %.1f | %.1f | %.1f | %.1f              |\n', ...
        i, fit_table(i,1), fit_table(i,2), fit_table(i,3), fit_table(i,4));
end
fprintf('\n');

%% 5. Figure 2: Bode Plot — True vs N4SID vs ssest
figure('Name', 'Figure 2: Bode Plot', ...
       'Position', [100, 100, 600, 500]);

w = logspace(-1, log10(pi/Ts), 300);
[mag_true, ph_true]     = bode(sys_true, w);
[mag_n4sid_b, ph_n4sid_b] = bode(sys_n4sid, w);
[mag_ssest_b, ph_ssest_b] = bode(sys_ssest, w);

for out_idx = 1:2
    % Magnitude
    subplot(2, 2, 2*(out_idx-1)+1);
    semilogx(w, 20*log10(squeeze(mag_true(out_idx,1,:))), ...
        'b-', 'LineWidth', 2); hold on;
    semilogx(w, 20*log10(squeeze(mag_n4sid_b(out_idx,1,:))), ...
        'r--', 'LineWidth', 1.2);
    semilogx(w, 20*log10(squeeze(mag_ssest_b(out_idx,1,:))), ...
        'g-.', 'LineWidth', 1.2);
    hold off; grid on;
    ylabel('Magnitude (dB)', 'FontSize', 11);
    title(sprintf('Output %d: Magnitude', out_idx), 'FontSize', 12);
    if out_idx == 1
        legend('True', 'N4SID', 'ssest', ...
            'FontSize', 9, 'Location', 'southwest');
    end
    
    % Phase
    subplot(2, 2, 2*(out_idx-1)+2);
    semilogx(w, squeeze(ph_true(out_idx,1,:)), ...
        'b-', 'LineWidth', 2); hold on;
    semilogx(w, squeeze(ph_n4sid_b(out_idx,1,:)), ...
        'r--', 'LineWidth', 1.2);
    semilogx(w, squeeze(ph_ssest_b(out_idx,1,:)), ...
        'g-.', 'LineWidth', 1.2);
    hold off; grid on;
    ylabel('Phase (deg)', 'FontSize', 11);
    title(sprintf('Output %d: Phase', out_idx), 'FontSize', 12);
    if out_idx == 2
        xlabel('Frequency (rad/s)', 'FontSize', 11);
    end
end

sgtitle('Bode Plot: True vs N4SID vs ssest (N4SID+PEM)', 'FontSize', 13);

%% 6. Figure 3: Residual Analysis (ssest model)
fprintf('=== Residual Analysis ===\n');

y_pred = predict(sys_ssest, data_val, 1);
residuals = y_val - y_pred.y;

max_lag = 25;
N_val = size(residuals, 1);
conf_99 = 2.576 / sqrt(N_val);

figure('Name', 'Figure 3: Residual Analysis', ...
       'Position', [100, 100, 700, 450]);

for i = 1:n_outputs
    % Autocorrelation
    [acf, lags] = xcorr(residuals(:,i), max_lag, 'coeff');
    
    subplot(n_outputs, 2, 2*(i-1)+1);
    stem(lags, acf, 'b', 'MarkerSize', 3, 'LineWidth', 0.8);
    hold on;
    yline(conf_99, 'r--', 'LineWidth', 1);
    yline(-conf_99, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('Lag', 'FontSize', 11);
    ylabel('Autocorrelation', 'FontSize', 11);
    title(sprintf('Output %d: Residual Autocorrelation', i), 'FontSize', 12);
    grid on;
    ylim([-0.15, 1.1]);
    
    % Cross-correlation (residual vs input)
    [ccf, lags_c] = xcorr(residuals(:,i), u_val, max_lag, 'coeff');
    
    subplot(n_outputs, 2, 2*(i-1)+2);
    stem(lags_c, ccf, 'b', 'MarkerSize', 3, 'LineWidth', 0.8);
    hold on;
    yline(conf_99, 'r--', 'LineWidth', 1);
    yline(-conf_99, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('Lag', 'FontSize', 11);
    ylabel('Cross-correlation', 'FontSize', 11);
    title(sprintf('Output %d: Residual-Input Cross-corr.', i), 'FontSize', 12);
    grid on;
end

sgtitle('Residual Analysis: ssest model (red = 99% confidence)', 'FontSize', 13);

%% 7. Figure 4: Pole-Zero Map — True vs Identified
figure('Name', 'Figure 4: Pole-Zero Map', ...
       'Position', [100, 100, 480, 420]);

% Unit circle
theta = linspace(0, 2*pi, 200);
plot(cos(theta), sin(theta), 'k-', 'LineWidth', 0.5); hold on;

% True poles
p_true = eig(A);
plot(real(p_true), imag(p_true), 'bo', ...
    'MarkerSize', 12, 'LineWidth', 2);

% N4SID poles
p_n4sid = eig(sys_n4sid.A);
plot(real(p_n4sid), imag(p_n4sid), 'rx', ...
    'MarkerSize', 10, 'LineWidth', 2);

% ssest poles
p_ssest = eig(sys_ssest.A);
plot(real(p_ssest), imag(p_ssest), 'g^', ...
    'MarkerSize', 9, 'LineWidth', 1.5);

hold off; grid on; axis equal;
xlabel('Real', 'FontSize', 12);
ylabel('Imaginary', 'FontSize', 12);
title('Pole Locations: True vs Identified', 'FontSize', 13);
legend('Unit circle', 'True', 'N4SID', 'ssest (N4SID+PEM)', ...
    'FontSize', 10, 'Location', 'northwest');
xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);

%% 8. Eigenvalue comparison (console)
fprintf('\n=== Eigenvalue Comparison ===\n');
fprintf('%-12s  %-12s  %-12s\n', 'True', 'N4SID', 'ssest');
eig_true_s  = sort(real(eig(A)));
eig_n4sid_s = sort(real(eig(sys_n4sid.A)));
eig_ssest_s = sort(real(eig(sys_ssest.A)));
for i = 1:n_true
    fprintf('%12.4f  %12.4f  %12.4f\n', ...
        eig_true_s(i), eig_n4sid_s(i), eig_ssest_s(i));
end

fprintf('\n=== All done. 4 figures generated for the blog. ===\n');
fprintf('Figure 1: Model Order Selection  → "Model Order Selection" section\n');
fprintf('Figure 2: Bode Plot              → "MATLAB Implementation" section\n');
fprintf('Figure 3: Residual Analysis      → "Model Validation" section\n');
fprintf('Figure 4: Pole-Zero Map          → "MATLAB Implementation" section\n');
