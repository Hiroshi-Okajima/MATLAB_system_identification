%% parametric_pem_demo.m
%  Classical Parametric System Identification: ARX, ARMAX, OE, BJ, and PEM
%
%  This script demonstrates the four classical polynomial model structures
%  and the prediction error method on a system with strongly colored noise.
%
%  Key design: The noise power is comparable to the signal power (low SNR),
%  so the noise model significantly affects the plant estimate in ARX/ARMAX
%  (where plant and noise share the denominator). This makes the difference
%  between the model structures clearly visible.
%
%  Blog: https://blog.control-theory.com/entry/parametric-identification
%  Hub:  https://blog.control-theory.com/entry/system-identification
%
%  Author: Hiroshi Okajima, Kumamoto University
%  Repository: https://github.com/Hiroshi-Okajima/MATLAB_system_identification

clear; close all; clc;

%% 1. Define the true system (Box-Jenkins structure)
%  Plant: G(q) = B(q)/F(q)
%  Noise: H(q) = C(q)/D(q)   (independent from plant)
%
%  The plant has poles at z = 0.6 +/- 0.15j (inside unit circle, stable).
%  The noise has a pole at z = 0.85 (slowly decaying colored noise).
%  Key: Plant denominator F(q) != Noise denominator D(q).
%       ARX and ARMAX force them to be the same => biased plant estimate.

Ts = 1;

% Plant: G(q) = B(q) / F(q)
B_true = [0 0.8 -0.3];           % B(q) = 0.8q^{-1} - 0.3q^{-2}
F_true = [1 -1.2 0.45];          % F(q) = 1 - 1.2q^{-1} + 0.45q^{-2}

% Noise: H(q) = C(q) / D(q)  (strongly colored, slow decay)
C_true = [1 0.7];                % C(q) = 1 + 0.7q^{-1}
D_true = [1 -0.85];              % D(q) = 1 - 0.85q^{-1}

% Create true BJ model
sys_true = idpoly(1, B_true, C_true, D_true, F_true, 0, Ts);

% Plant transfer function for comparison
G_true_tf = tf(B_true, F_true, Ts);

fprintf('=== True System (Box-Jenkins) ===\n');
fprintf('Plant: G(q) = B(q)/F(q)\n');
fprintf('  B(q) = 0.8q^{-1} - 0.3q^{-2}\n');
fprintf('  F(q) = 1 - 1.2q^{-1} + 0.45q^{-2}\n');
fprintf('Noise: H(q) = C(q)/D(q)\n');
fprintf('  C(q) = 1 + 0.7q^{-1}\n');
fprintf('  D(q) = 1 - 0.85q^{-1}   [pole: 0.85 — slow decay]\n');
fprintf('Key: F(q) != D(q) => ARX/ARMAX will be biased.\n\n');

%% 2. Generate input-output data
N = 2000;
rng(42);
u = randn(N, 1);

% Control the SNR
y_clean = lsim(G_true_tf, u);
signal_power = var(y_clean);

H_true_tf = tf(C_true, D_true, Ts);
H_dc = abs(dcgain(H_true_tf));
e_std = sqrt(0.5 * signal_power) / max(H_dc, 1);
e_std = max(e_std, 0.3);

fprintf('Signal power (plant output): %.2f\n', signal_power);
fprintf('Innovation std: %.2f\n', e_std);

e = e_std * randn(N, 1);
y = sim(sys_true, [u e]);

data = iddata(y, u, Ts);
data_est = data(1:1200);
data_val = data(1201:end);

y_val_clean = lsim(G_true_tf, u(1201:end));
noise_val = y(1201:end) - y_val_clean;
SNR_dB = 10*log10(var(y_val_clean) / var(noise_val));
fprintf('Approximate SNR: %.1f dB\n\n', SNR_dB);

%% 3. Estimate ARX model
fprintf('--- ARX Estimation ---\n');
na = 2; nb = 2; nk = 1;
sys_arx = arx(data_est, [na nb nk]);
fprintf('  ARX(%d,%d,%d) estimated.\n', na, nb, nk);

%% 4. Estimate ARMAX model
fprintf('--- ARMAX Estimation ---\n');
nc = 1;
sys_armax = armax(data_est, [na nb nc nk]);
fprintf('  ARMAX(%d,%d,%d,%d) estimated.\n', na, nb, nc, nk);

%% 5. Estimate Output-Error model
fprintf('--- OE Estimation ---\n');
nf = 2;
sys_oe = oe(data_est, [nb nf nk]);
fprintf('  OE(%d,%d,%d) estimated.\n', nb, nf, nk);

%% 6. Estimate Box-Jenkins model
fprintf('--- BJ Estimation ---\n');
nd = 1;
sys_bj = bj(data_est, [nb nc nd nf nk]);
fprintf('  BJ(%d,%d,%d,%d,%d) estimated.\n', nb, nc, nd, nf, nk);

%% 7. Simulation comparison (plant model accuracy)
fprintf('\n=== Simulation Fit (%%) — tests plant model accuracy ===\n');

rng(999);
u_sim = randn(800, 1);
y_sim_true = lsim(G_true_tf, u_sim);

G_arx   = tf(sys_arx.B, sys_arx.A, Ts);
G_armax = tf(sys_armax.B, sys_armax.A, Ts);
G_oe    = tf(sys_oe.B, sys_oe.F, Ts);
G_bj    = tf(sys_bj.B, sys_bj.F, Ts);

y_sim_arx   = lsim(G_arx, u_sim);
y_sim_armax = lsim(G_armax, u_sim);
y_sim_oe    = lsim(G_oe, u_sim);
y_sim_bj    = lsim(G_bj, u_sim);

calc_fit = @(y_hat, y_ref) max(0, ...
    (1 - norm(y_ref - y_hat) / norm(y_ref - mean(y_ref))) * 100);

fit_arx   = calc_fit(y_sim_arx, y_sim_true);
fit_armax = calc_fit(y_sim_armax, y_sim_true);
fit_oe    = calc_fit(y_sim_oe, y_sim_true);
fit_bj    = calc_fit(y_sim_bj, y_sim_true);

fprintf('ARX:   %.1f%%\n', fit_arx);
fprintf('ARMAX: %.1f%%\n', fit_armax);
fprintf('OE:    %.1f%%\n', fit_oe);
fprintf('BJ:    %.1f%%\n', fit_bj);

% Plot simulation comparison
figure('Name', 'Simulation Comparison (Plant Model Accuracy)');
t_sim = (0:length(u_sim)-1)' * Ts;
plot(t_sim, y_sim_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t_sim, y_sim_arx, 'b--');
plot(t_sim, y_sim_oe, 'm:');
plot(t_sim, y_sim_bj, 'g-.');
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Output', 'FontSize', 12);
title('Simulation: True Plant Output vs Identified Models', 'FontSize', 13);
legend(sprintf('True G (reference)'), ...
       sprintf('ARX (%.1f%%)', fit_arx), ...
       sprintf('OE (%.1f%%)', fit_oe), ...
       sprintf('BJ (%.1f%%)', fit_bj), ...
       'FontSize', 11, 'Location', 'best');
xlim([0 100]);
grid on;

%% 8. 1-step prediction comparison
fprintf('\n=== 1-Step Prediction Fit (%%) — tests plant + noise model ===\n');

opt_pred = compareOptions('InitialCondition', 'z');
[~, pfit_arx]   = compare(data_val, sys_arx, 1, opt_pred);
[~, pfit_armax] = compare(data_val, sys_armax, 1, opt_pred);
[~, pfit_oe]    = compare(data_val, sys_oe, 1, opt_pred);
[~, pfit_bj]    = compare(data_val, sys_bj, 1, opt_pred);

fprintf('ARX:   %.1f%%\n', pfit_arx);
fprintf('ARMAX: %.1f%%\n', pfit_armax);
fprintf('OE:    %.1f%%\n', pfit_oe);
fprintf('BJ:    %.1f%%\n', pfit_bj);

%% 9. Bode plot comparison (manual plot — avoids bodeplot/legend issue)
figure('Name', 'Bode Plot: True Plant vs Identified');

w = logspace(-2, log10(pi/Ts), 300);

[mag_true, ~]   = bode(G_true_tf, w);
[mag_arx, ~]    = bode(G_arx, w);
[mag_armax, ~]  = bode(G_armax, w);
[mag_oe, ~]     = bode(G_oe, w);
[mag_bj, ~]     = bode(G_bj, w);

mag_true  = squeeze(mag_true);
mag_arx   = squeeze(mag_arx);
mag_armax = squeeze(mag_armax);
mag_oe    = squeeze(mag_oe);
mag_bj    = squeeze(mag_bj);

semilogx(w, 20*log10(mag_true), 'k-', 'LineWidth', 2); hold on;
semilogx(w, 20*log10(mag_arx),   'b--', 'LineWidth', 1);
semilogx(w, 20*log10(mag_armax), 'r-.', 'LineWidth', 1);
semilogx(w, 20*log10(mag_oe),    'm:',  'LineWidth', 1.5);
semilogx(w, 20*log10(mag_bj),    'g-',  'LineWidth', 1);
hold off;
grid on;
xlabel('Frequency (rad/s)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title('Plant Transfer Function G: True vs Identified', 'FontSize', 13);
legend('True G', 'ARX (biased)', 'ARMAX (biased)', 'OE', 'BJ (best)', ...
    'FontSize', 10, 'Location', 'southwest');

%% 10. Residual analysis
figure('Name', 'Residuals: ARX vs BJ');
subplot(2,1,1);
resid(data_val, sys_arx);
title('Residuals: ARX (shared denominator)', 'FontSize', 12);
subplot(2,1,2);
resid(data_val, sys_bj);
title('Residuals: BJ (independent plant/noise)', 'FontSize', 12);

%% 11. ssest
fprintf('\n=== ssest (N4SID + PEM) ===\n');
sys_ssest = ssest(data_est, 2, 'Ts', Ts);  % Force discrete-time
% Simulate using the identified state-space model (noise-free)
sys_ssest_ss = ss(sys_ssest.A, sys_ssest.B, sys_ssest.C, sys_ssest.D, Ts);
y_sim_ssest = lsim(sys_ssest_ss, u_sim);
fit_ssest = calc_fit(y_sim_ssest, y_sim_true);
fprintf('ssest simulation fit: %.1f%%\n', fit_ssest);

%% 12. Summary
fprintf('\n============================================================\n');
fprintf('  Summary: Model Comparison\n');
fprintf('============================================================\n');
fprintf('%-8s  %8s  %8s  %s\n', 'Model', 'Sim(%)', 'Pred(%)', 'Notes');
fprintf('%-8s  %7.1f%%  %7.1f%%  %s\n', 'ARX',   fit_arx,   pfit_arx,   'G=B/A, H=1/A (shared denom)');
fprintf('%-8s  %7.1f%%  %7.1f%%  %s\n', 'ARMAX', fit_armax, pfit_armax, 'G=B/A, H=C/A (shared denom)');
fprintf('%-8s  %7.1f%%  %7.1f%%  %s\n', 'OE',    fit_oe,    pfit_oe,    'G=B/F, H=1   (no noise model)');
fprintf('%-8s  %7.1f%%  %7.1f%%  %s\n', 'BJ',    fit_bj,    pfit_bj,    'G=B/F, H=C/D (fully independent)');
fprintf('%-8s  %7.1f%%  %8s  %s\n',     'ssest', fit_ssest, '---',      'N4SID + PEM (state-space)');
fprintf('============================================================\n');
fprintf('\nExpected pattern:\n');
fprintf('  Simulation fit:  OE, BJ >> ARX, ARMAX\n');
fprintf('    (OE and BJ have independent F(q), so plant model is unbiased)\n');
fprintf('  Prediction fit:  BJ > ARMAX > ARX > OE\n');
fprintf('    (BJ has the best noise model; OE has none)\n');
fprintf('  Bode plot: OE and BJ match true G closely;\n');
fprintf('    ARX and ARMAX show bias (distorted by noise model)\n');

fprintf('\nDone. See figures for visual results.\n');
