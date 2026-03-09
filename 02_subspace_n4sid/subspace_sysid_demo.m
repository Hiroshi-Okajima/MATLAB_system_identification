%% subspace_sysid_demo.m
%  Subspace Identification: N4SID, MOESP, and CVA Comparison
%
%  This script demonstrates the three major subspace identification
%  algorithms on the same dataset and compares their results.
%
%  Blog: https://blog.control-theory.com/entry/subspace-identification
%  Hub:  https://blog.control-theory.com/entry/system-identification
%
%  Author: Hiroshi Okajima, Kumamoto University
%  Repository: https://github.com/Hiroshi-Okajima/MATLAB_system_identification

clear; close all; clc;

%% 1. Define the true plant (3rd-order, 1 input, 2 outputs)
%  All three modes are well-separated and contribute significantly
%  to the output, so SVD clearly shows 3 significant singular values.
Ts = 0.1;

% 3rd-order state-space: three distinct real poles at 0.9, 0.5, -0.3
%  (well-separated => clear 3rd-order signature in SVD)
A = diag([0.9, 0.5, -0.3]);
B = [1; 1; 1];
C = [1 1 0.5;
     0.5 -1 1];
D = [0; 0];

sys_true = ss(A, B, C, D, Ts);
n_true = size(A, 1);  % True system order = 3

fprintf('=== True Plant ===\n');
fprintf('Order: %d,  Inputs: %d,  Outputs: %d\n', ...
    n_true, size(B,2), size(C,1));
fprintf('Poles: %.1f, %.1f, %.1f (well-separated)\n\n', eig(A));

%% 2. Generate input-output data
N = 3000;
rng(123);
u = randn(N, 1);

y_clean = lsim(sys_true, u);
noise_std = 0.05;
e = noise_std * randn(N, size(C,1));
y = y_clean + e;

data = iddata(y, u, Ts);
data_est = data(1:2000);
data_val = data(2001:end);

fprintf('Data: %d samples (est: %d, val: %d)\n', N, 2000, 1000);
fprintf('Noise std: %.3f\n\n', noise_std);

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
title('Model Order Selection via Singular Value Analysis', 'FontSize', 13);
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
figure('Name', 'Validation: N4SID vs MOESP vs CVA');
compare(data_val, sys_n4sid, sys_moesp, sys_cva);
% Set English legend explicitly
h = findobj(gcf, 'Type', 'axes');
for k = 1:length(h)
    lgd = findobj(h(k), 'Type', 'legend');
    if ~isempty(lgd)
        lgd.String = {'Measured', 'N4SID', 'MOESP', 'CVA'};
        lgd.FontSize = 10;
    end
end
sgtitle('Validation: N4SID vs MOESP vs CVA', 'FontSize', 13);

% Compute fit percentages
[~, fit_n4sid] = compare(data_val, sys_n4sid);
[~, fit_moesp] = compare(data_val, sys_moesp);
[~, fit_cva]   = compare(data_val, sys_cva);

fprintf('\n=== Validation Fit (%%) ===\n');
for i = 1:size(C, 1)
    fprintf('Output %d:  N4SID = %.1f%%,  MOESP = %.1f%%,  CVA = %.1f%%\n', ...
        i, fit_n4sid(i), fit_moesp(i), fit_cva(i));
end

%% 6. PEM refinement using ssest (N4SID + PEM)
fprintf('\n=== ssest: N4SID initialization + PEM refinement ===\n');
sys_ssest = ssest(data_est, n);
[~, fit_ssest] = compare(data_val, sys_ssest);

for i = 1:size(C, 1)
    fprintf('Output %d:  ssest = %.1f%%\n', i, fit_ssest(i));
end

%% 7. Bode plot comparison (using bode + manual plot to avoid legend issue)
figure('Name', 'Bode Plot: True vs Identified');

% Compute frequency response manually
w = logspace(-1, log10(pi/Ts), 200);
[mag_true, ph_true] = bode(sys_true, w);
[mag_n4sid, ph_n4sid] = bode(sys_n4sid, w);
[mag_ssest, ph_ssest] = bode(sys_ssest, w);

% Plot magnitude for output 1, input 1
subplot(2,1,1);
semilogx(w, 20*log10(squeeze(mag_true(1,1,:))), 'b-', 'LineWidth', 1.5); hold on;
semilogx(w, 20*log10(squeeze(mag_n4sid(1,1,:))), 'r--', 'LineWidth', 1);
semilogx(w, 20*log10(squeeze(mag_ssest(1,1,:))), 'g-.', 'LineWidth', 1);
hold off; grid on;
ylabel('Magnitude (dB)', 'FontSize', 12);
title('Bode Plot: Output 1 / Input 1', 'FontSize', 13);
legend('True', 'N4SID', 'ssest (N4SID+PEM)', 'FontSize', 10, 'Location', 'southwest');

% Plot phase for output 1, input 1
subplot(2,1,2);
semilogx(w, squeeze(ph_true(1,1,:)), 'b-', 'LineWidth', 1.5); hold on;
semilogx(w, squeeze(ph_n4sid(1,1,:)), 'r--', 'LineWidth', 1);
semilogx(w, squeeze(ph_ssest(1,1,:)), 'g-.', 'LineWidth', 1);
hold off; grid on;
xlabel('Frequency (rad/s)', 'FontSize', 12);
ylabel('Phase (deg)', 'FontSize', 12);
legend('True', 'N4SID', 'ssest (N4SID+PEM)', 'FontSize', 10, 'Location', 'southwest');

%% 8. Residual analysis
figure('Name', 'Residual Analysis (ssest)');
resid(data_val, sys_ssest);
sgtitle('Residual Analysis: ssest model', 'FontSize', 13);

%% 9. Eigenvalue comparison
fprintf('\n=== Eigenvalue Comparison ===\n');
eig_true  = sort(real(eig(A)));
eig_n4sid = sort(real(eig(sys_n4sid.A)));
eig_ssest = sort(real(eig(sys_ssest.A)));

fprintf('%-12s  %-12s  %-12s\n', 'True', 'N4SID', 'ssest');
for i = 1:n_true
    fprintf('%12.4f  %12.4f  %12.4f\n', eig_true(i), eig_n4sid(i), eig_ssest(i));
end

fprintf('\nDone. See figures for visual results.\n');
