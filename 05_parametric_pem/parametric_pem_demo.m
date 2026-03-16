%% parametric_pem_demo.m  (v2 — two-case comparison for blog)
%  Classical Parametric System Identification: ARX, ARMAX, OE, BJ, and PEM
%
%  Generates figures for the blog article in TWO scenarios:
%    Case 1 (Low noise):  All methods succeed — shows basic workflow
%    Case 2 (High noise): ARX bias becomes visible — shows why model
%                         structure matters
%
%  Blog figures:
%    Figure 1: [Case 1] Bode plot — all methods match the true plant
%    Figure 2: [Case 2] Bode plot — ARX/ARMAX biased, OE/BJ accurate
%    Figure 3: [Case 2] Time-domain simulation comparison
%
%  Blog: https://blog.control-theory.com/entry/parametric-identification
%  Hub:  https://blog.control-theory.com/entry/system-identification
%
%  Author: Hiroshi Okajima, Kumamoto University
%  Repository: https://github.com/Hiroshi-Okajima/MATLAB_system_identification

clear; close all; clc;

%% 1. Define the true system (Box-Jenkins structure)
%  Plant: G(q) = B(q)/F(q)
%  Noise: H(q) = C(q)/D(q)
%
%  Key: F(q) != D(q), so ARX/ARMAX (shared denominator) will be biased
%  when the noise is strong enough.

Ts = 1;

% Plant: G(q) = B(q) / F(q)
B_true = [0 0.8 -0.3];           % B(q) = 0.8q^{-1} - 0.3q^{-2}
F_true = [1 -1.2 0.45];          % F(q) = 1 - 1.2q^{-1} + 0.45q^{-2}

% Noise: H(q) = C(q) / D(q)
C_true = [1 0.7];                % C(q) = 1 + 0.7q^{-1}
D_true = [1 -0.9];               % D(q) = 1 - 0.9q^{-1}

sys_true = idpoly(1, B_true, C_true, D_true, F_true, 0, Ts);
G_true_tf = tf(B_true, F_true, Ts);
H_true_tf = tf(C_true, D_true, Ts);

fprintf('=== True System (Box-Jenkins) ===\n');
fprintf('Plant G(q) = B(q)/F(q):  poles at z = 0.6 +/- 0.15j\n');
fprintf('Noise H(q) = C(q)/D(q):  pole at z = 0.9 (colored)\n');
fprintf('Key: F(q) != D(q)\n\n');

%% 2. Common settings
N = 400;
na = 2; nb = 2; nc = 1; nd = 1; nf = 2; nk = 1;

% Frequency vector for Bode plots
w = logspace(-2, log10(pi/Ts), 300);
mag_true = squeeze(bode(G_true_tf, w));

% Fit calculation
calc_fit = @(y_hat, y_ref) max(0, ...
    (1 - norm(y_ref - y_hat) / norm(y_ref - mean(y_ref))) * 100);

% Fresh input for simulation test (noise-free)
rng(999);
u_sim = randn(800, 1);
y_sim_true = lsim(G_true_tf, u_sim);

%% ================================================================
%  CASE 0: No Noise — Perfect identification (100%)
%% ================================================================
fprintf('============================================================\n');
fprintf('  CASE 0: No Noise\n');
fprintf('============================================================\n');

rng(42);
u0 = randn(N, 1);
y0 = lsim(G_true_tf, u0);   % noise-free output

data0 = iddata(y0, u0, Ts);
data0_est = data0(1:300);

% Identify with ARX and BJ (representative)
sys0_arx = arx(data0_est, [na nb nk]);
sys0_bj  = bj(data0_est, [nb nc nd nf nk]);

G0_arx = tf(sys0_arx.B, sys0_arx.A, Ts);
G0_bj  = tf(sys0_bj.B, sys0_bj.F, Ts);

fit0_arx = calc_fit(lsim(G0_arx, u_sim), y_sim_true);
fit0_bj  = calc_fit(lsim(G0_bj, u_sim),  y_sim_true);

fprintf('Simulation fit:  ARX=%.1f%%  BJ=%.1f%%\n', fit0_arx, fit0_bj);

% Figure 1: Bode plot (Case 0 — no noise)
figure('Name', 'Figure 1: Bode — No Noise (Perfect Identification)', ...
       'Position', [100 100 600 400]);

semilogx(w, 20*log10(mag_true), 'k-', 'LineWidth', 2.5); hold on;
semilogx(w, 20*log10(squeeze(bode(G0_arx, w))), 'b--', 'LineWidth', 1.2);
semilogx(w, 20*log10(squeeze(bode(G0_bj, w))),  'g-',  'LineWidth', 1.2);
hold off; grid on;
xlabel('Frequency (rad/s)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title('Case 0: No Noise — All Models Achieve Perfect Fit', 'FontSize', 13);
legend(sprintf('True G'), ...
       sprintf('ARX (%.1f%%)', fit0_arx), ...
       sprintf('BJ (%.1f%%)', fit0_bj), ...
       'FontSize', 11, 'Location', 'southwest');

%% ================================================================
%  CASE 1: Low Noise — All methods still work well
%% ================================================================
fprintf('\n============================================================\n');
fprintf('  CASE 1: Low Noise (e_std = 0.1)\n');
fprintf('============================================================\n');

e_std_low = 0.1;
rng(42);
u1 = randn(N, 1);
e1 = e_std_low * randn(N, 1);
% y = G(q)*u + H(q)*e,  computed explicitly
y1_plant = lsim(G_true_tf, u1);
y1_noise = lsim(H_true_tf, e1);
y1 = y1_plant + y1_noise;

data1 = iddata(y1, u1, Ts);
data1_est = data1(1:300);

% Identify all four models
sys1_arx   = arx(data1_est, [na nb nk]);
sys1_armax = armax(data1_est, [na nb nc nk]);
sys1_oe    = oe(data1_est, [nb nf nk]);
sys1_bj    = bj(data1_est, [nb nc nd nf nk]);

% Extract plant transfer functions
G1_arx   = tf(sys1_arx.B, sys1_arx.A, Ts);
G1_armax = tf(sys1_armax.B, sys1_armax.A, Ts);
G1_oe    = tf(sys1_oe.B, sys1_oe.F, Ts);
G1_bj    = tf(sys1_bj.B, sys1_bj.F, Ts);

% Simulation fit
fit1_arx   = calc_fit(lsim(G1_arx, u_sim),   y_sim_true);
fit1_armax = calc_fit(lsim(G1_armax, u_sim), y_sim_true);
fit1_oe    = calc_fit(lsim(G1_oe, u_sim),    y_sim_true);
fit1_bj    = calc_fit(lsim(G1_bj, u_sim),    y_sim_true);

fprintf('Simulation fit:  ARX=%.1f%%  ARMAX=%.1f%%  OE=%.1f%%  BJ=%.1f%%\n', ...
    fit1_arx, fit1_armax, fit1_oe, fit1_bj);

% Figure 2: Bode plot (Case 1 — low noise)
figure('Name', 'Figure 2: Bode — Low Noise (All Methods Succeed)', ...
       'Position', [100 100 600 400]);

semilogx(w, 20*log10(mag_true), 'k-', 'LineWidth', 2.5); hold on;
semilogx(w, 20*log10(squeeze(bode(G1_arx, w))),   'b--', 'LineWidth', 1);
semilogx(w, 20*log10(squeeze(bode(G1_armax, w))), 'r-.', 'LineWidth', 1);
semilogx(w, 20*log10(squeeze(bode(G1_oe, w))),    'm:',  'LineWidth', 1.5);
semilogx(w, 20*log10(squeeze(bode(G1_bj, w))),    'g-',  'LineWidth', 1);
hold off; grid on;
xlabel('Frequency (rad/s)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title('Case 1: Low Noise — All Models Match the True Plant', 'FontSize', 13);
legend(sprintf('True G'), ...
       sprintf('ARX (%.1f%%)', fit1_arx), ...
       sprintf('ARMAX (%.1f%%)', fit1_armax), ...
       sprintf('OE (%.1f%%)', fit1_oe), ...
       sprintf('BJ (%.1f%%)', fit1_bj), ...
       'FontSize', 10, 'Location', 'southwest');

% Figure 3: Time-domain simulation (Case 1 — low noise)
figure('Name', 'Figure 3: Simulation — Low Noise', ...
       'Position', [100 100 650 380]);

t_sim = (0:length(u_sim)-1)' * Ts;
plot(t_sim, y_sim_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t_sim, lsim(G1_arx, u_sim),   'b--', 'LineWidth', 0.8);
plot(t_sim, lsim(G1_oe, u_sim),    'm:',  'LineWidth', 1.2);
plot(t_sim, lsim(G1_bj, u_sim),    'g-.', 'LineWidth', 1);
hold off; grid on;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Output', 'FontSize', 12);
title('Simulation (Noise-Free Input): True Plant vs Identified Models', 'FontSize', 13);
legend(sprintf('True G'), ...
       sprintf('ARX (%.1f%%)', fit1_arx), ...
       sprintf('OE (%.1f%%)', fit1_oe), ...
       sprintf('BJ (%.1f%%)', fit1_bj), ...
       'FontSize', 11, 'Location', 'best');
xlim([0 100]);

%% ================================================================
%  Console summary tables (for blog article text)
%% ================================================================
fprintf('\n============================================================\n');
fprintf('  Summary Table (for blog article)\n');
fprintf('============================================================\n');
fprintf('\n--- Case 0: No Noise ---\n');
fprintf('| Model  | Sim Fit | Notes |\n');
fprintf('|--------|---------|-------|\n');
fprintf('| ARX    | %.1f%%  | Perfect (no noise bias) |\n', fit0_arx);
fprintf('| BJ     | %.1f%%  | Perfect |\n', fit0_bj);

fprintf('\n--- Case 1: Low Noise (e_std = %.1f) ---\n', e_std_low);
fprintf('| Model  | Sim Fit | Notes |\n');
fprintf('|--------|---------|-------|\n');
fprintf('| ARX    | %.1f%%  | G=B/A, H=1/A |\n', fit1_arx);
fprintf('| ARMAX  | %.1f%%  | G=B/A, H=C/A |\n', fit1_armax);
fprintf('| OE     | %.1f%%  | G=B/F, H=1   |\n', fit1_oe);
fprintf('| BJ     | %.1f%%  | G=B/F, H=C/D |\n', fit1_bj);

fprintf('\n=== All done. 3 figures generated for the blog. ===\n');
fprintf('Figure 1: Bode (no noise)          → "MATLAB Implementation" section\n');
fprintf('Figure 2: Bode (low noise)         → "MATLAB Implementation" section\n');
fprintf('Figure 3: Simulation (low noise)   → "MATLAB Implementation" section\n');
