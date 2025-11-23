function full_wave_rectifier_CT()
% RUN_THRYRISTOR_MODEL  Convert of provided Python code into a single MATLAB file.
% Save this file as run_thyristor_model.m and call run_thyristor_model().

%% Main parameters (same as Python __main__)
Vrms = 60;
E = 12;
esr = 4.256;
capacity = 100;
v_drop = 1.2;
tr = 10*10^(-6);  % rise time (s)
tf = 5*10^(-6);   % fall time (s)
freq = 50;        % Hz
ileak = 0.02;     % A
Vm = sqrt(2)*Vrms;

% find alpha1 and beta considering v_drop
[alpha1, beta] = find_range_of_alphas_and_beta_with_v_drop(E, v_drop, Vm);

% make alpha array
alphas = make_alpha_array(alpha1, beta, 0.01);

% x axis for phase plots
x = linspace(0, 2*pi, 10000);

% Plots and computations
plot_time_of_charging(alphas, Vm, E, esr, capacity, beta, true, 'log', Vrms);
plot_final_soc(alphas, Vm, E, esr, 30.0, capacity, 0.5, beta, true);
plot_i_thyristor_with_t(x, alpha1, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
plot_Vak_with_t(x, alpha1, beta, Vm, E, v_drop, tr, tf, freq);
plot_power_loss_with_t(x, Vm, E, esr, ileak, v_drop, tr, tf, freq, alpha1, beta);
plot_power_loss_with_alpha(alphas, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);

end

%% ---------------------- Helper functions ----------------------

function t = norm_angle(theta)
% Normalize angle to [0,2*pi)
t = mod(theta, 2*pi);
end

function [t1, t2] = find_range_of_alphas_and_beta(E, Vm)
% Return two angles t1,t2 in [0,2*pi) satisfying sin(theta)=E/Vm
frac = E / Vm;
if abs(frac) > 1.0
    error('|E/Vm| must be <= 1. Got %g', frac);
end
a = asin(frac);
t1 = norm_angle(a);
t2 = norm_angle(pi - a);
end

function [t1, t2] = find_range_of_alphas_and_beta_with_v_drop(E, v_drop, Vm)
% Return two angles t1,t2 in [0,2*pi) satisfying sin(theta)=(E+v_drop)/Vm
frac = (E + v_drop) / Vm;
if abs(frac) > 1.0
    error('|E+v_drop|/Vm must be <= 1. Got %g', frac);
end
a = asin(frac);
t1 = norm_angle(a);
t2 = norm_angle(pi - a);
end

function arr = make_alpha_array(a1, a2, step)
% Create array from a1 to a2 (exclusive) with step, handle wrap-around
a1 = double(a1); a2 = double(a2); step = double(step);
if step <= 0
    error('step must be positive');
end
two_pi = 2*pi;
if abs(a1 - a2) < eps
    arr = a1;
    return;
end
if a1 < a2
    arr = a1:step:(a2-step);
    if isempty(arr)
        arr = a1; % fallback
    end
    return;
end
% wrap-around
part1 = a1:step:(two_pi-step);
part2 = 0:step:(a2-step);
if ~isempty(part1) && ~isempty(part2)
    arr = [part1, part2];
elseif ~isempty(part1)
    arr = part1;
else
    arr = part2;
end
end

function y = vs_t(x, Vm)
y = Vm.*sin(x);
end

function val = i_t(x, alpha, beta, Vm, E, esr)
% current without leakage and voltage drop going to battery
if alpha <= x && x <= beta
    val = (Vm*sin(x) - E) / esr;
elseif pi + alpha <= x && x <= pi + beta
    val = (-Vm*sin(x) - E) / esr;
else
    val = 0.0;
end
end

function val = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop)
if alpha <= x && x <= beta
    val = (Vm*sin(x) - E - v_drop) / esr;
elseif pi + alpha <= x && x <= pi + beta
    val = (-Vm*sin(x) - v_drop - E) / esr;
else
    val = -ileak;
end
end

function val = i_battery_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
% Convert rise/fall times to angular ranges
omega = 2*pi*freq;
d_theta_r = omega * tr;
d_theta_f = omega * tf;

natural_I = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop);

% Region 1
if x < alpha
    val = -ileak; return;
end

% Region 2
if alpha <= x && x < alpha + d_theta_r
    [m,b] = line_from_points(alpha, -ileak, alpha + d_theta_r, ...
        i_t_including_leakage_simplified(alpha + d_theta_r, alpha, beta, Vm, E, esr, ileak, v_drop));
    val = m*x + b; return;
end

% Region 3
if alpha + d_theta_r <= x && x < beta
    val = natural_I; return;
end

% Region 4
if beta <= x && x < beta + d_theta_f
    [m,b] = line_from_points(beta, i_t_including_leakage_simplified(beta, alpha, beta, Vm, E, esr, ileak, v_drop), ...
        beta + d_theta_f, -ileak);
    val = m*x + b; return;
end

% Region 5
if beta + d_theta_f <= x && x < pi + alpha
    val = -ileak; return;
end

% Region 6
if pi + alpha <= x && x < pi + alpha + d_theta_r
    [m,b] = line_from_points(pi + alpha, -ileak, pi + alpha + d_theta_r, ...
        i_t_including_leakage_simplified(pi + alpha + d_theta_r, alpha, beta, Vm, E, esr, ileak, v_drop));
    val = m*x + b; return;
end

% Region 7
if pi + alpha + d_theta_r <= x && x < pi + beta
    val = natural_I; return;
end

% Region 8
if pi + beta <= x && x < pi + beta + d_theta_f
    [m,b] = line_from_points(pi + beta, i_t_including_leakage_simplified(pi + beta, alpha, beta, Vm, E, esr, ileak, v_drop), ...
        pi + beta + d_theta_f, -ileak);
    val = m*x + b; return;
end

% Region 9
val = -ileak;
end

function val = i_thyristor_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
% same shape as battery current but simpler (thyristor 1)
omega = 2*pi*freq;
d_theta_r = omega * tr;
d_theta_f = omega * tf;

natural_I = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop);

if x < alpha
    val = -ileak; return;
end
if alpha <= x && x < alpha + d_theta_r
    [m,b] = line_from_points(alpha, -ileak, alpha + d_theta_r, ...
        i_t_including_leakage_simplified(alpha + d_theta_r, alpha, beta, Vm, E, esr, ileak, v_drop));
    val = m*x + b; return;
end
if alpha + d_theta_r <= x && x < beta
    val = natural_I; return;
end
if beta <= x && x < beta + d_theta_f
    [m,b] = line_from_points(beta, i_t_including_leakage_simplified(beta, alpha, beta, Vm, E, esr, ileak, v_drop), ...
        beta + d_theta_f, -ileak);
    val = m*x + b; return;
end
val = -ileak;
end

function val = Vak_including_drop_simplified(x, alpha, beta, Vm, E, v_drop)
if alpha <= x && x <= beta
    val = v_drop;
else
    val = Vm*sin(x) - E;
end
end

function plot_Vak_with_t_simplified(x, alpha, beta, Vm, E, v_drop)
y = arrayfun(@(xi) Vak_including_drop_simplified(xi, alpha, beta, Vm, E, v_drop), x);
figure; plot(x, y);
xlabel('Angle (rad)'); ylabel('Vak (V)');
title('Vak including drop simplified vs Angle (\omega t)');
grid on;
end

function [m,b] = line_from_points(x1, y1, x2, y2)
m = (y2 - y1) / (x2 - x1);
b = y1 - m * x1;
end

function val = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, tr, tf, freq)
omega = 2*pi*freq;
d_theta_r = omega * tr;
d_theta_f = omega * tf;

% Region 1
if x < alpha
    val = vs_t(x, Vm) - E; return;
end

% Region 2: alpha -> alpha + d_theta_f (falling)
if alpha <= x && x < alpha + d_theta_f
    [m,b] = line_from_points(alpha, vs_t(alpha, Vm) - E, alpha + d_theta_f, v_drop);
    val = m*x + b; return;
end

% Region 3: constant v_drop
if alpha + d_theta_f <= x && x < beta
    val = v_drop; return;
end

% Region 4: beta -> beta + d_theta_r (rising)
if beta <= x && x < beta + d_theta_r
    [m,b] = line_from_points(beta, v_drop, beta + d_theta_r, vs_t(beta + d_theta_r, Vm) - E);
    val = m*x + b; return;
end

% Region 5:
if beta + d_theta_r <= x && x < pi + alpha
    val = vs_t(x, Vm) - E; return;
end

% Region 6:
if pi + alpha <= x && x < pi + alpha + d_theta_f
    % note: original uses 2*vs_t(...)+v_drop as target
    [m,b] = line_from_points(pi + alpha, vs_t(pi + alpha, Vm) - E, pi + alpha + d_theta_f, 2*vs_t(pi + alpha + d_theta_f, Vm) + v_drop);
    val = m*x + b; return;
end

% Region 7:
if pi + alpha + d_theta_f <= x && x < pi + beta
    val = 2*vs_t(x, Vm) + v_drop; return;
end

% Region 8:
if pi + beta <= x && x < pi + beta + d_theta_r
    [m,b] = line_from_points(pi + beta, 2*vs_t(pi + beta, Vm) + v_drop, pi + beta + d_theta_r, vs_t(pi + beta + d_theta_r, Vm) - E);
    val = m*x + b; return;
end

% Region 9:
val = vs_t(x, Vm) - E;
end

function plot_Vak_with_t(x, alpha, beta, Vm, E, v_drop, tr, tf, freq)
y = arrayfun(@(xi) Vak_including_drop(xi, alpha, beta, Vm, E, v_drop, tr, tf, freq), x);
figure; plot(x, y);
xlabel('Angle (\omega t) rad'); ylabel('Vak (V)');
title('Vak including loss vs Angle (\omega t)');
grid on;
end

function plot_i_without_leakage(x, alpha, beta, Vm, E, esr)
y = arrayfun(@(xi) i_t(xi, alpha, beta, Vm, E, esr), x);
figure; plot(x, y);
xlabel('Angle (\omega t) rad'); ylabel('Current i (A)');
title('Current i without loss vs Angle (\omega t)');
grid on;
end

function plot_i_thyristor_with_t(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
y = arrayfun(@(xi) i_thyristor_t_including_leakage(xi, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq), x);
figure; plot(x, y);
xlabel('Angle (\omega t) rad'); ylabel('Current i (A)');
title('Current i thyristor including leakage vs Angle (\omega t)');
grid on;
end

function Iavg = compute_Iavg(Vm, E, esr, alpha, beta)
% integrate i_t over 0..2pi
f = @(xx) arrayfun(@(xval) i_t(xval, alpha, beta, Vm, E, esr), xx);
result = integral(f, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6, 'AbsTol',1e-12);
Iavg = (1.0 / (2.0 * pi)) * result;
end

function Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
f = @(xx) arrayfun(@(xval) i_battery_t_including_leakage(xval, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq).^2, xx);
integral_val = integral(f, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6, 'AbsTol',1e-12);
Irms = sqrt(integral_val / (2.0 * pi));
end

function pl = compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
Vak_thyrestor_1 = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, tr, tf, freq);
i_thyristor_1 = i_thyristor_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
resistance_loss = esr * (i_battery_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq).^2);
pl = 2*(Vak_thyrestor_1 .* i_thyristor_1) + resistance_loss;
end

function P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
f = @(xx) arrayfun(@(xval) compute_power_loss(xval, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq), xx);
result = integral(f, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6, 'AbsTol',1e-12);
P_loss_avg = (1.0 / (2.0 * pi)) * result;
end

function plot_power_loss_with_alpha(alpha_array, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
y = zeros(size(alpha_array));
for k = 1:numel(alpha_array)
    alpha = alpha_array(k);
    % Irms is computed but not used in original except maybe for reference
    Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq); %#ok<NASGU>
    P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
    y(k) = P_loss_avg;
end
figure; plot(alpha_array, y);
xlabel('Firing Angle (rad)'); ylabel('Average Power Loss (W)');
title('Average Power Loss vs Firing Angle'); grid on;
end

function plot_power_loss_with_t(x, Vm, E, esr, ileak, v_drop, tr, tf, freq, alpha, beta)
y = arrayfun(@(xi) compute_power_loss(xi, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq), x);
figure; plot(x, y);
xlabel('Angle (\omega t) rad'); ylabel('Power Loss (W)');
title('Power Loss vs Angle (\omega t)'); grid on;
end

function charging_time = compute_charging_time(Iavg, C)
% returns NaN for nonpositive Iavg
charging_time = nan(size(Iavg));
if ~isempty(C)
    pos = Iavg > 0;
    charging_time(pos) = C ./ Iavg(pos);
end
end

function plot_time_of_charging(alphas, Vm, E, esr, C, beta, degrees, charging_time_scale, Vrms)
% Plot Iavg (left axis) and charging_time (right axis).
if nargin < 11
    Vrms = [];
end
if degrees
    x = rad2deg(alphas);
    xlabel_text = 'Firing angle (deg)';
else
    x = alphas;
    xlabel_text = 'Firing angle (rad)';
end

[Iavg, charging_time] = time_of_charging(Vm, E, esr, C, alphas, beta, 0.01);

figure('Position',[100 100 1000 600]);
yyaxis left;
plot(x, Iavg, '-o', 'MarkerSize',4);
ylabel('Iavg (A)');
xlabel(xlabel_text);
grid on;

if ~all(isnan(charging_time))
    yyaxis right;
    plot(x, charging_time, '-s', 'MarkerSize',4);
    ylabel_text = 'Charging time (h)';
    if ~isempty(charging_time_scale) && ~strcmp(charging_time_scale,'linear')
        set(gca, 'YScale', charging_time_scale);
        ylabel_text = sprintf('%s (%s)', ylabel_text, charging_time_scale);
    end
    ylabel(ylabel_text);
end

title(sprintf('Iavg and charging time vs firing angle (Vrms=%g V, E=%g V, esr=%g ohms, C=%g Ah)', Vrms, E, esr, C), 'Interpreter','none');
end

function new_soc = compute_final_soc(Iavg, initial_soc, C, required_charging_time)
% initial_soc in percent
new_soc = (Iavg .* required_charging_time  + (initial_soc/100.0)*C) ./ C;
new_soc = new_soc * 100.0;
new_soc = min(max(new_soc, 0.0), 100.0);
end

function plot_final_soc(alphas, Vm, E, esr, initial_soc, C, required_charging_time, beta, degrees)
if degrees
    x = rad2deg(alphas);
    xlabel_text = 'Firing angle (deg)';
else
    x = alphas;
    xlabel_text = 'Firing angle (rad)';
end

[Iavg, charging_time] = time_of_charging(Vm, E, esr, C, alphas, beta, 0.01);

new_soc = compute_final_soc(Iavg, initial_soc, C, required_charging_time);

figure;
plot(x, new_soc, '-o', 'MarkerSize',4);
xlabel(xlabel_text); ylabel('Final State of Charge (%)');
grid on;
title(sprintf('Final SOC vs Firing Angle\nInit SOC=%g%%, Time=%gh, C=%gAh', initial_soc, required_charging_time, C), 'Interpreter','none');
end

function [Iavg_array, charging_time] = time_of_charging(Vm, E, esr, C, alphas, beta, ~)
% Compute Iavg for each alpha in alphas
Iavg_array = zeros(size(alphas));
for k = 1:numel(alphas)
    Iavg_array(k) = compute_Iavg(Vm, E, esr, alphas(k), beta);
end
charging_time = compute_charging_time(Iavg_array, C);
end
