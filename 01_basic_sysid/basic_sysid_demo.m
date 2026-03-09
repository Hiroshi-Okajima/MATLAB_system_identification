%% basic_sysid_demo.m
%  Introduction to System Identification in MATLAB
%
%  This script demonstrates the basic workflow of system identification:
%    1. Data generation (true plant + noise)
%    2. Model estimation (ARX)
%    3. Model validation (compare, resid, bode)
%
%  Blog: https://blog.control-theory.com/entry/2024/10/03/151451
%  Hub:  https://blog.control-theory.com/entry/system-identification
%
%  Author: Hiroshi Okajima, Kumamoto University
%  Repository: https://github.com/Hiroshi-Okajima/MATLAB_system_identification

clear; close all; clc;

%% 1. Define the true plant
%  2nd-order discrete-time transfer function:
%    G(z) = (0.1z + 0.05) / (z^2 - 1.5z + 0.7)
%  Sampling period: Ts = 0.1 s

Ts = 0.1;
num = [0.1 0.05];
den = [1 -1.5 0.7];
G_true = tf(num, den, Ts);

fprintf('=== True Plant ===\n');
G_true

%% 2. Generate input-output data
N = 1000;                        % Number of data points
rng(42);                         % Fix random seed for reproducibility
u = randn(N, 1);                 % White noise input (persistent excitation)

% Simulate plant output with measurement noise
y_clean = lsim(G_true, u);
noise_std = 0.05;
e = noise_std * randn(N, 1);     % Measurement noise
y = y_clean + e;

% Create iddata object
data = iddata(y, u, Ts);

%% 3. Split into estimation and validation sets
data_est = data(1:600);          % First 600 points for estimation
data_val = data(601:end);        % Remaining 400 points for validation

fprintf('Estimation data: %d samples\n', length(data_est.y));
fprintf('Validation data: %d samples\n', length(data_val.y));

%% 4. Estimate ARX model
%  ARX model: A(q)y(k) = B(q)u(k) + e(k)
%  Try model orders: na = 2, nb = 2, nk = 1

na = 2;  nb = 2;  nk = 1;
sys_arx = arx(data_est, [na nb nk]);

fprintf('\n=== Identified ARX Model ===\n');
sys_arx

%% 5. Model validation

% 5a. Compare predicted output with validation data
figure('Name', 'Model Validation: compare');
compare(data_val, sys_arx);
title('ARX Model vs Validation Data');

% 5b. Residual analysis
figure('Name', 'Residual Analysis');
resid(data_val, sys_arx);

% 5c. Bode plot comparison
figure('Name', 'Bode Plot Comparison');
bode(G_true, 'b-', sys_arx, 'r--');
legend('True Plant', 'Identified ARX');
title('Bode Plot: True vs Identified');

%% 6. Compare identified parameters with true plant
%  Convert ARX to transfer function for comparison
[num_id, den_id] = tfdata(sys_arx, 'v');

fprintf('\n=== Parameter Comparison ===\n');
fprintf('True numerator:   [%.4f  %.4f]\n', num);
fprintf('Identified num:   [%.4f  %.4f]\n', num_id(2:end));
fprintf('True denominator: [%.4f  %.4f  %.4f]\n', den);
fprintf('Identified den:   [%.4f  %.4f  %.4f]\n', den_id);

fprintf('\nDone. See figures for visual results.\n');
