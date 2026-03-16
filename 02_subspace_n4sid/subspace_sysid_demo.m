%% subspace_sysid_demo.m  (v2 — close poles + colored noise)
%  Subspace Identification: N4SID, MOESP, and CVA Comparison
%
%  This script demonstrates the three major subspace identification
%  algorithms on the same dataset and compares their results.
%  The plant has close poles and the noise is colored, creating
%  a challenging scenario where the three methods show distinct behavior.
%
%  Blog: https://blog.control-theory.com/entry/subspace-identification
%  Hub:  https://blog.control-theory.com/entry/system-identification
%
%  Author: Hiroshi Okajima, Kumamoto University
%  Repository: https://github.com/Hiroshi-Okajima/MATLAB_system_identification

clear; close all; clc;

%% 1. Define the true plant (3rd-order, 1 input, 2 outputs)
%  Close poles make identification harder and reveal algorithmic differences.
Ts = 0.1;

% 3rd-order state-space: three close real poles at 0.85, 0.80, 0.75
A = diag([0.85, 0.80, 0.75]);
B = [1; 0.8; 0.5];
C = [1 1 0.5;
     0.5 -1 1];
D = [0; 0];

sys_true = ss(A, B, C, D, Ts);
n_true = size(A, 1);  % True system order = 3

fprintf('=== True Plant ===\n');
fprintf('Order: %d,  Inputs: %d,  Outputs: %d\n', ...
    n_true, size(B,2), size(C,1));
fprintf('Poles: %.2f, %.2f, %.2f (close together)\n\n', eig(A));

%% 2. Generate input-output data with colored noise
N = 3000;
rng(42);
u = randn(N, 1);

y_clean = lsim(sys_true, u);

% Colored noise: white noise filtered through a 1st-order system
%   H(z) = 1 / (1 - 0.8 z^{-1})  — strongly colored
noise_std = 0.3;
e_white = noise_std * randn(N, size(C,1));
noise_filter = tf(1, [1, -0.8], Ts);
e_colored = zeros(N, size(C,1));
for i = 1:size(C,1)
    e_colored(:,i) = lsim(noise_filter, e_white(:,i));
end
y = y_clean + e_colored;

data = iddata(y, u, Ts);
data_est = data(1:2000);
data_val = data(2001:end);

fprintf('Data: %d samples (est: %d, val: %d)\n', N, 2000, 1000);
fprintf('Noise std: %.2f (colored, pole at 0.8)\n\n', noise_std);

%% 3. Model order selection (no GUI)
fprintf('=== Model Order Selection ===\n');

max_order = 8;
loss = zeros(max_order, 1);
for n_test = 1:max_order
    sys_test = n4sid(data_est, n_test);
    loss(n_test) = sys_test.Report.Fit.LossFcn;
end

% Find the elbow: largest relative drop in loss
rel_drop = -diff(loss) ./ abs(loss(1:end-1));
[~, n_recommended] = max(rel_drop);

% Plot
figure('Name', 'Model Order Selection');
colors = repmat([0.3 0.5 0.8], max_order, 1);
colors(n_recommended, :) = [0.9 0.3 0.2];  % Highlight recommended
b = bar(1:max_order, loss);
b.FaceColor = 'flat';
b.CData = colors;
xlabel('Model Order n', 'FontSize', 12);
ylabel('Loss Function (Prediction Error Variance)', 'FontSize', 12);
title('Model Order Selection via Loss Function', 'FontSize', 13);
grid on;
set(gca, 'XTick', 1:max_order);
text(n_recommended, loss(n_recommended)*1.05, ...
    sprintf('\\leftarrow n = %d (recommended)', n_recommended), ...
    'FontSize', 11, 'Color', [0.9 0.3 0.2]);

fprintf('Recommended order: %d (true order = %d)\n\n', n_recommended, n_true);

%% 4. Identify using the three subspace methods
n = n_true;  % Use true order for fair comparison

opt_n4sid = n4sidOptions('N4Weight', 'auto');
sys_n4sid = n4sid(data_est, n, opt_n4sid);

opt_moesp = n4sidOptions('N4Weight', 'MOESP');
sys_moesp = n4sid(data_est, n, opt_moesp);

opt_cva = n4sidOptions('N4Weight', 'CVA');
sys_cva = n4sid(data_est, n, opt_cva);

fprintf('=== Identification Complete (order = %d) ===\n', n);

%% 5. Compare on validation data
u_val = data_val.u;
y_val = data_val.y;
t_val = (0:length(u_val)-1)' * Ts;

y_n4sid = lsim(sys_n4sid, u_val);
y_moesp = lsim(sys_moesp, u_val);
y_cva   = lsim(sys_cva, u_val);

% Compute NRMSE fit
calc_fit = @(y_hat, y_ref) max(0, ...
    (1 - norm(y_ref - y_hat) / norm(y_ref - mean(y_ref))) * 100);

fit_n4sid = zeros(size(C,1), 1);
fit_moesp = zeros(size(C,1), 1);
fit_cva   = zeros(size(C,1), 1);
for i = 1:size(C,1)
    fit_n4sid(i) = calc_fit(y_n4sid(:,i), y_val(:,i));
    fit_moesp(i) = calc_fit(y_moesp(:,i), y_val(:,i));
    fit_cva(i)   = calc_fit(y_cva(:,i), y_val(:,i));
end

fprintf('\n=== Validation Fit (%%) ===\n');
for i = 1:size(C, 1)
    fprintf('Output %d:  N4SID = %.1f%%,  MOESP = %.1f%%,  CVA = %.1f%%\n', ...
        i, fit_n4sid(i), fit_moesp(i), fit_cva(i));
end

% Plot validation comparison
figure('Name', 'Validation: N4SID vs MOESP vs CVA');
for i = 1:size(C,1)
    subplot(size(C,1), 1, i);
    plot(t_val, y_val(:,i), 'k-', 'LineWidth', 1); hold on;
    plot(t_val, y_n4sid(:,i), 'r--', 'LineWidth', 0.8);
    plot(t_val, y_moesp(:,i), 'b-.', 'LineWidth', 0.8);
    plot(t_val, y_cva(:,i), 'g:', 'LineWidth', 1);
    hold off; grid on;
    ylabel(sprintf('Output %d', i), 'FontSize', 12);
    title(sprintf('Output %d:  N4SID=%.1f%%,  MOESP=%.1f%%,  CVA=%.1f%%', ...
        i, fit_n4sid(i), fit_moesp(i), fit_cva(i)), 'FontSize', 12);
    if i == 1
        legend('Measured', 'N4SID', 'MOESP', 'CVA', ...
            'FontSize', 10, 'Location', 'northeast');
    end
    xlim([0 30]);  % Show first 30 seconds for clarity
end
xlabel('Time (s)', 'FontSize', 12);
sgtitle('Validation: N4SID vs MOESP vs CVA (Simulation)', 'FontSize', 13);

%% 6. PEM refinement using ssest (N4SID + PEM)
fprintf('\n=== ssest: N4SID initialization + PEM refinement ===\n');
sys_ssest = ssest(data_est, n, 'Ts', Ts);

y_ssest = lsim(sys_ssest, u_val);
fit_ssest = zeros(size(C,1), 1);
for i = 1:size(C,1)
    fit_ssest(i) = calc_fit(y_ssest(:,i), y_val(:,i));
end

for i = 1:size(C, 1)
    fprintf('Output %d:  ssest = %.1f%%\n', i, fit_ssest(i));
end

%% 7. Bode plot comparison
figure('Name', 'Bode Plot: True vs Identified');

w = logspace(-1, log10(pi/Ts), 200);
[mag_true, ph_true] = bode(sys_true, w);
[mag_n4sid_b, ph_n4sid_b] = bode(sys_n4sid, w);
[mag_moesp_b, ph_moesp_b] = bode(sys_moesp, w);
[mag_cva_b, ph_cva_b] = bode(sys_cva, w);
[mag_ssest_b, ph_ssest_b] = bode(sys_ssest, w);

% Plot magnitude for output 1, input 1
subplot(2,1,1);
semilogx(w, 20*log10(squeeze(mag_true(1,1,:))), 'k-', 'LineWidth', 2); hold on;
semilogx(w, 20*log10(squeeze(mag_n4sid_b(1,1,:))), 'r--', 'LineWidth', 1);
semilogx(w, 20*log10(squeeze(mag_moesp_b(1,1,:))), 'b-.', 'LineWidth', 1);
semilogx(w, 20*log10(squeeze(mag_cva_b(1,1,:))), 'g:', 'LineWidth', 1.2);
semilogx(w, 20*log10(squeeze(mag_ssest_b(1,1,:))), 'm-', 'LineWidth', 1);
hold off; grid on;
ylabel('Magnitude (dB)', 'FontSize', 12);
title('Bode Plot: Output 1 / Input 1', 'FontSize', 13);
legend('True', 'N4SID', 'MOESP', 'CVA', 'ssest (N4SID+PEM)', ...
    'FontSize', 9, 'Location', 'southwest');

% Plot phase for output 1, input 1
subplot(2,1,2);
semilogx(w, squeeze(ph_true(1,1,:)), 'k-', 'LineWidth', 2); hold on;
semilogx(w, squeeze(ph_n4sid_b(1,1,:)), 'r--', 'LineWidth', 1);
semilogx(w, squeeze(ph_moesp_b(1,1,:)), 'b-.', 'LineWidth', 1);
semilogx(w, squeeze(ph_cva_b(1,1,:)), 'g:', 'LineWidth', 1.2);
semilogx(w, squeeze(ph_ssest_b(1,1,:)), 'm-', 'LineWidth', 1);
hold off; grid on;
xlabel('Frequency (rad/s)', 'FontSize', 12);
ylabel('Phase (deg)', 'FontSize', 12);
legend('True', 'N4SID', 'MOESP', 'CVA', 'ssest (N4SID+PEM)', ...
    'FontSize', 9, 'Location', 'southwest');

%% 8. Residual analysis
fprintf('\n=== Residual Analysis ===\n');

y_pred = predict(sys_ssest, data_val, 1);
residuals = y_val - y_pred.y;

max_lag = 25;
figure('Name', 'Residual Analysis (ssest)');

for i = 1:size(C,1)
    [acf, lags] = xcorr(residuals(:,i), max_lag, 'coeff');
    
    subplot(size(C,1), 2, 2*(i-1)+1);
    stem(lags, acf, 'b', 'MarkerSize', 3);
    xlabel('Lag', 'FontSize', 11);
    ylabel('Autocorrelation', 'FontSize', 11);
    title(sprintf('Output %d: Residual Autocorrelation', i), 'FontSize', 12);
    grid on;
    N_val = length(residuals);
    conf = 2.576 / sqrt(N_val);
    hold on;
    plot(xlim, [conf conf], 'r--', 'LineWidth', 0.8);
    plot(xlim, [-conf -conf], 'r--', 'LineWidth', 0.8);
    hold off;
    
    [ccf, lags_c] = xcorr(residuals(:,i), u_val, max_lag, 'coeff');
    
    subplot(size(C,1), 2, 2*(i-1)+2);
    stem(lags_c, ccf, 'b', 'MarkerSize', 3);
    xlabel('Lag', 'FontSize', 11);
    ylabel('Cross-correlation', 'FontSize', 11);
    title(sprintf('Output %d: Residual-Input Cross-corr.', i), 'FontSize', 12);
    grid on;
    hold on;
    plot(xlim, [conf conf], 'r--', 'LineWidth', 0.8);
    plot(xlim, [-conf -conf], 'r--', 'LineWidth', 0.8);
    hold off;
end

sgtitle('Residual Analysis: ssest model (99% confidence bounds in red)', 'FontSize', 13);

%% 9. Eigenvalue comparison
fprintf('\n=== Eigenvalue Comparison ===\n');
eig_true  = sort(real(eig(A)));
eig_n4sid = sort(real(eig(sys_n4sid.A)));
eig_moesp = sort(real(eig(sys_moesp.A)));
eig_cva   = sort(real(eig(sys_cva.A)));
eig_ssest = sort(real(eig(sys_ssest.A)));

fprintf('%-12s  %-12s  %-12s  %-12s  %-12s\n', ...
    'True', 'N4SID', 'MOESP', 'CVA', 'ssest');
for i = 1:n_true
    fprintf('%12.4f  %12.4f  %12.4f  %12.4f  %12.4f\n', ...
        eig_true(i), eig_n4sid(i), eig_moesp(i), eig_cva(i), eig_ssest(i));
end

%% 10. Summary table
fprintf('\n=== Summary: Validation Fit (%%) ===\n');
fprintf('%-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Output', 'N4SID', 'MOESP', 'CVA', 'ssest');
for i = 1:size(C,1)
    fprintf('%-10d  %8.1f%%  %8.1f%%  %8.1f%%  %8.1f%%\n', ...
        i, fit_n4sid(i), fit_moesp(i), fit_cva(i), fit_ssest(i));
end

fprintf('\nDone. See figures for visual results.\n');
