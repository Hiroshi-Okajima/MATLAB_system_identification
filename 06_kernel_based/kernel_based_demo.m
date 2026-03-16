%% Kernel-Based System Identification: Regularized Impulse Response Estimation
%  =========================================================================
%  Blog article:
%    https://blog.control-theory.com/entry/kernel-based-identification
%
%  Repository: MATLAB_system_identification / 06_kernel_based
%
%  This script demonstrates kernel-based regularized system identification
%  using MATLAB's impulseest function. It generates the following figures:
%
%    Fig 1: Impulse response comparison
%           (True vs TC kernel vs SS kernel vs No Reg. vs ARX)
%    Fig 2: Kernel matrices visualization (TC, SS, DC)
%    Fig 3: Effect of data length — short data (N=80) vs long data (N=200)
%    Fig 4: Step response comparison
%           (True vs Regularized FIR vs ARX)
%
%  Requirements:
%    - MATLAB R2020b or later
%    - System Identification Toolbox
%    - Control System Toolbox
%
%  Author: Hiroshi Okajima, Kumamoto University
%  =========================================================================

clear; close all; clc;

%% ========================================================================
%  Plot Style (self-contained, per project instructions)
%  ========================================================================
% Colors
COL_TRUE   = [0.00 0.00 0.00];  % #000000 — true / reference (black)
COL_METH1  = [0.12 0.47 0.71];  % #1f77b4 — method 1 (blue)
COL_METH2  = [0.84 0.15 0.16];  % #d62728 — method 2 (red)
COL_METH3  = [0.17 0.63 0.17];  % #2ca02c — method 3 (green)
COL_METH4  = [1.00 0.50 0.05];  % #ff7f0e — method 4 (orange)
COL_CTRL   = [0.58 0.40 0.74];  % #9467bd — control input (purple)
COL_NOISE  = [0.50 0.50 0.50];  % #7f7f7f — noise / baseline (gray)

% Line widths
LW_TRUE = 2.0;
LW_EST  = 1.5;
LW_BASE = 1.5;

% Figure sizes (inches) — optimized for Hatena Blog (~700 px width)
FIGSIZE_STD   = [8, 5];    % standard time-series
FIGSIZE_WIDE  = [8, 6];    % Bode plot etc.
FIGSIZE_SQ    = [6, 6];    % pole map etc.

% Font
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 13);

% DPI for PNG export
DPI = 150;

% Create fig directory
fig_dir = 'fig';
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ========================================================================
%  True System
%  ========================================================================
%  Discrete-time transfer function: G(z) = (z + 0.5) / (z^2 - 1.5z + 0.7)
%  (Same as blog article "Complete Workflow Example")
Ts = 1;
G_true = tf([1 0.5], [1 -1.5 0.7], Ts);

% True impulse response (for comparison)
n_imp = 50;  % truncation length for impulse response
[g_true, t_imp] = impulse(G_true, 0:Ts:(n_imp-1)*Ts);
g_true = g_true(:);

fprintf('True system poles: %.4f +/- %.4fj\n', ...
    real(pole(G_true)), abs(imag(pole(G_true))));
fprintf('True system: G(z) = (z + 0.5) / (z^2 - 1.5z + 0.7)\n\n');

%% ========================================================================
%  Data Generation
%  ========================================================================
rng(42);  % for reproducibility

N = 200;                          % number of data points
sigma_e = 0.5;                    % noise standard deviation (large enough to show differences)
u = randn(N, 1);                  % random input (white noise)
y_clean = lsim(G_true, u);       % clean output
e = sigma_e * randn(N, 1);       % measurement noise
y = y_clean + e;                  % noisy output

data = iddata(y, u, Ts);

fprintf('Data: N = %d, sigma_e = %.2f, SNR ≈ %.1f dB\n', ...
    N, sigma_e, 10*log10(var(y_clean)/var(e)));

%% ========================================================================
%  Figure 1: Impulse Response Comparison
%  ========================================================================
fprintf('\n--- Figure 1: Impulse Response Comparison ---\n');

% (a) Regularized impulse response — TC kernel (default)
opt_tc = impulseestOptions('RegularizationKernel', 'TC');
sys_tc = impulseest(data, n_imp, opt_tc);
[g_tc, ~] = impulse(sys_tc, 0:Ts:(n_imp-1)*Ts);
g_tc = g_tc(:);
fprintf('  TC kernel: done\n');

% (b) Regularized impulse response — SS (stable spline) kernel
opt_ss = impulseestOptions('RegularizationKernel', 'SS');
sys_ss = impulseest(data, n_imp, opt_ss);
[g_ss, ~] = impulse(sys_ss, 0:Ts:(n_imp-1)*Ts);
g_ss = g_ss(:);
fprintf('  SS kernel: done\n');

% (c) No regularization (ordinary least squares FIR)
opt_none = impulseestOptions('RegularizationKernel', 'none');
sys_none = impulseest(data, n_imp, opt_none);
[g_none, ~] = impulse(sys_none, 0:Ts:(n_imp-1)*Ts);
g_none = g_none(:);
fprintf('  No regularization: done\n');

% (d) Parametric ARX for comparison
sys_arx = arx(data, [2 2 1]);
G_arx = tf(sys_arx.B, sys_arx.A, Ts);
[g_arx, ~] = impulse(G_arx, 0:Ts:(n_imp-1)*Ts);
g_arx = g_arx(:);
fprintf('  ARX(2,2,1): done\n');

% Plot
fig1 = figure('Units', 'inches', 'Position', [1 1 FIGSIZE_STD]);
hold on; grid on;
stem(t_imp, g_true, 'filled', 'Color', COL_TRUE, ...
    'LineWidth', LW_TRUE, 'MarkerSize', 5, 'DisplayName', 'True');
plot(t_imp, g_tc, '-', 'Color', COL_METH1, ...
    'LineWidth', LW_EST, 'DisplayName', 'TC Kernel');
plot(t_imp, g_ss, '-', 'Color', COL_METH2, ...
    'LineWidth', LW_EST, 'DisplayName', 'SS Kernel');
plot(t_imp, g_none, '--', 'Color', COL_NOISE, ...
    'LineWidth', LW_BASE, 'DisplayName', 'No Reg. (LS)');
plot(t_imp, g_arx, '-', 'Color', COL_METH3, ...
    'LineWidth', LW_EST, 'DisplayName', 'ARX(2,2,1)');
hold off;

xlabel('Time step $k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Impulse response $g_k$', 'Interpreter', 'latex', 'FontSize', 14);
title('Kernel-Based vs. Unregularized Impulse Response', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
set(gca, 'GridAlpha', 0.3);
xlim([0 (n_imp-1)*Ts]);

% Export
exportgraphics(fig1, fullfile(fig_dir, 'kernel_impulse_comparison.png'), ...
    'Resolution', DPI, 'BackgroundColor', 'white');
fprintf('  Saved: fig/kernel_impulse_comparison.png\n');

%% ========================================================================
%  Figure 2: Kernel Matrix Visualization
%  ========================================================================
fprintf('\n--- Figure 2: Kernel Matrix Visualization ---\n');

n_kern = 30;  % kernel size for visualization

% Hyperparameters (representative values)
lambda = 1.0;
alpha  = 0.85;
rho    = 0.5;   % for DC kernel

% TC kernel: K_TC(i,j) = lambda * alpha^max(i,j)
K_TC = zeros(n_kern);
for i = 1:n_kern
    for j = 1:n_kern
        K_TC(i,j) = lambda * alpha^max(i,j);
    end
end

% SS kernel: K_SS(i,j) = lambda * alpha^max(i,j)/2 * (alpha^max(i,j) - alpha^|i-j|/3)
K_SS = zeros(n_kern);
for i = 1:n_kern
    for j = 1:n_kern
        mx = max(i,j);
        K_SS(i,j) = lambda * (alpha^mx / 2) * (alpha^mx - alpha^abs(i-j) / 3);
    end
end

% DC kernel: K_DC(i,j) = lambda * alpha^((i+j)/2) * rho^|i-j|
K_DC = zeros(n_kern);
for i = 1:n_kern
    for j = 1:n_kern
        K_DC(i,j) = lambda * alpha^((i+j)/2) * rho^abs(i-j);
    end
end

fig2 = figure('Units', 'inches', 'Position', [1 1 FIGSIZE_WIDE(1) FIGSIZE_WIDE(2)+1]);

subplot(1,3,1);
imagesc(K_TC); axis equal tight; colorbar;
title('TC Kernel', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('$j$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$i$', 'Interpreter', 'latex', 'FontSize', 14);

subplot(1,3,2);
imagesc(K_SS); axis equal tight; colorbar;
title('SS Kernel', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('$j$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$i$', 'Interpreter', 'latex', 'FontSize', 14);

subplot(1,3,3);
imagesc(K_DC); axis equal tight; colorbar;
title(['DC Kernel ($\rho=' num2str(rho) '$)'], ...
    'Interpreter', 'latex', 'FontSize', 13);
xlabel('$j$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$i$', 'Interpreter', 'latex', 'FontSize', 14);

colormap('parula');

% Export
exportgraphics(fig2, fullfile(fig_dir, 'kernel_matrices.png'), ...
    'Resolution', DPI, 'BackgroundColor', 'white');
fprintf('  Saved: fig/kernel_matrices.png\n');

%% ========================================================================
%  Figure 3: Effect of Data Length (Short Data vs Long Data)
%  ========================================================================
fprintf('\n--- Figure 3: Effect of Data Length ---\n');

N_short = 80;
N_long  = 200;
n_imp_short = 30;   % reduced FIR order for short data (must be < N_short)
n_imp_long  = n_imp; % = 50

% Short data
rng(42);
u_s = randn(N_short, 1);
y_s = lsim(G_true, u_s) + sigma_e * randn(N_short, 1);
data_short = iddata(y_s, u_s, Ts);

% Long data (reuse existing)
data_long = data;

% True impulse response truncated for each case
[g_true_s, t_imp_s] = impulse(G_true, 0:Ts:(n_imp_short-1)*Ts);
g_true_s = g_true_s(:);
[g_true_l, t_imp_l] = impulse(G_true, 0:Ts:(n_imp_long-1)*Ts);
g_true_l = g_true_l(:);

% Estimation with TC kernel
opt_tc = impulseestOptions('RegularizationKernel', 'TC');
sys_tc_short = impulseest(data_short, n_imp_short, opt_tc);
sys_tc_long  = impulseest(data_long,  n_imp_long,  opt_tc);
[g_tc_s, ~] = impulse(sys_tc_short, 0:Ts:(n_imp_short-1)*Ts);
[g_tc_l, ~] = impulse(sys_tc_long,  0:Ts:(n_imp_long-1)*Ts);
g_tc_s = g_tc_s(:); g_tc_l = g_tc_l(:);

% No regularization
opt_none = impulseestOptions('RegularizationKernel', 'none');
sys_ls_short = impulseest(data_short, n_imp_short, opt_none);
sys_ls_long  = impulseest(data_long,  n_imp_long,  opt_none);
[g_ls_s, ~] = impulse(sys_ls_short, 0:Ts:(n_imp_short-1)*Ts);
[g_ls_l, ~] = impulse(sys_ls_long,  0:Ts:(n_imp_long-1)*Ts);
g_ls_s = g_ls_s(:); g_ls_l = g_ls_l(:);

fprintf('  Short data (N=%d, n=%d): done\n', N_short, n_imp_short);
fprintf('  Long  data (N=%d, n=%d): done\n', N_long, n_imp_long);

% Compute fit metrics (NRMSE)
fit_tc_s = 100 * (1 - norm(g_tc_s - g_true_s) / norm(g_true_s - mean(g_true_s)));
fit_tc_l = 100 * (1 - norm(g_tc_l - g_true_l) / norm(g_true_l - mean(g_true_l)));
fit_ls_s = 100 * (1 - norm(g_ls_s - g_true_s) / norm(g_true_s - mean(g_true_s)));
fit_ls_l = 100 * (1 - norm(g_ls_l - g_true_l) / norm(g_true_l - mean(g_true_l)));

fprintf('  Fit (NRMSE%%): TC short=%.1f, TC long=%.1f, LS short=%.1f, LS long=%.1f\n', ...
    fit_tc_s, fit_tc_l, fit_ls_s, fit_ls_l);

fig3 = figure('Units', 'inches', 'Position', [1 1 FIGSIZE_STD(1) FIGSIZE_STD(2)+2]);

% (a) Short data
subplot(2,1,1);
hold on; grid on;
stem(t_imp_s, g_true_s, 'filled', 'Color', COL_TRUE, ...
    'LineWidth', LW_TRUE, 'MarkerSize', 4, 'DisplayName', 'True');
plot(t_imp_s, g_tc_s, '-', 'Color', COL_METH1, ...
    'LineWidth', LW_EST, 'DisplayName', ...
    sprintf('TC Kernel (fit: %.1f%%)', fit_tc_s));
plot(t_imp_s, g_ls_s, '--', 'Color', COL_NOISE, ...
    'LineWidth', LW_BASE, 'DisplayName', ...
    sprintf('No Reg. (fit: %.1f%%)', fit_ls_s));
hold off;
xlabel('Time step $k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$g_k$', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('(a) Short Data: $N = %d$, FIR order $n = %d$', N_short, n_imp_short), ...
    'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'GridAlpha', 0.3);
xlim([0 (n_imp_short-1)*Ts]);

% (b) Long data
subplot(2,1,2);
hold on; grid on;
stem(t_imp_l, g_true_l, 'filled', 'Color', COL_TRUE, ...
    'LineWidth', LW_TRUE, 'MarkerSize', 4, 'DisplayName', 'True');
plot(t_imp_l, g_tc_l, '-', 'Color', COL_METH1, ...
    'LineWidth', LW_EST, 'DisplayName', ...
    sprintf('TC Kernel (fit: %.1f%%)', fit_tc_l));
plot(t_imp_l, g_ls_l, '--', 'Color', COL_NOISE, ...
    'LineWidth', LW_BASE, 'DisplayName', ...
    sprintf('No Reg. (fit: %.1f%%)', fit_ls_l));
hold off;
xlabel('Time step $k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$g_k$', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('(b) Long Data: $N = %d$, FIR order $n = %d$', N_long, n_imp_long), ...
    'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'GridAlpha', 0.3);
xlim([0 (n_imp_long-1)*Ts]);

% Export
exportgraphics(fig3, fullfile(fig_dir, 'kernel_data_length_effect.png'), ...
    'Resolution', DPI, 'BackgroundColor', 'white');
fprintf('  Saved: fig/kernel_data_length_effect.png\n');

%% ========================================================================
%  Figure 4: Step Response — Regularized FIR → State-Space Conversion
%  ========================================================================
fprintf('\n--- Figure 4: Step Response Comparison ---\n');

% Convert regularized FIR to state-space and reduce order
sys_fir_tc = impulseest(data, 70, opt_tc);
sys_ss_red = balred(idss(sys_fir_tc), 2);  % reduce to order 2

% Step responses
t_step = 0:Ts:30;
[y_step_true, ~] = step(G_true, t_step);
[y_step_tc,   ~] = step(sys_fir_tc, t_step);
[y_step_ss,   ~] = step(sys_ss_red, t_step);
[y_step_arx,  ~] = step(G_arx, t_step);

fig4 = figure('Units', 'inches', 'Position', [1 1 FIGSIZE_STD]);
hold on; grid on;
stairs(t_step, y_step_true, '-',  'Color', COL_TRUE,  ...
    'LineWidth', LW_TRUE, 'DisplayName', 'True System');
stairs(t_step, y_step_tc,   '-',  'Color', COL_METH1, ...
    'LineWidth', LW_EST,  'DisplayName', 'TC Kernel (FIR, n=70)');
stairs(t_step, y_step_ss,   '-',  'Color', COL_METH2, ...
    'LineWidth', LW_EST,  'DisplayName', 'Reduced SS (order 2)');
stairs(t_step, y_step_arx,  '-',  'Color', COL_METH3, ...
    'LineWidth', LW_EST,  'DisplayName', 'ARX(2,2,1)');
hold off;

xlabel('Time step $k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Step response', 'FontSize', 14);
title('Step Response: Regularized FIR vs. Parametric ARX', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
set(gca, 'GridAlpha', 0.3);
xlim([0 30]);

% Export
exportgraphics(fig4, fullfile(fig_dir, 'kernel_step_response.png'), ...
    'Resolution', DPI, 'BackgroundColor', 'white');
fprintf('  Saved: fig/kernel_step_response.png\n');

%% ========================================================================
%  Summary
%  ========================================================================
fprintf('\n============================\n');
fprintf('  All figures saved to fig/\n');
fprintf('============================\n');
fprintf('  1. kernel_impulse_comparison.png\n');
fprintf('     — True vs TC vs SS vs No Reg. vs ARX impulse response\n');
fprintf('  2. kernel_matrices.png\n');
fprintf('     — Visualization of TC, SS, DC kernel matrices\n');
fprintf('  3. kernel_data_length_effect.png\n');
fprintf('     — Short data (N=80) vs Long data (N=200)\n');
fprintf('  4. kernel_step_response.png\n');
fprintf('     — Step response: Regularized FIR vs Parametric ARX\n');
fprintf('============================\n');
