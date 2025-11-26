function full_wave_rectifier_CT(Vrms, f, E, R, capacity, opts)
    arguments
        Vrms double = 0.0
        f double = 0.0
        E double = 0.0
        R double = 0.0
        capacity double = 0.0
        opts.Vforward double = 0.0
        opts.Ileak double = 0.0
        opts.Itr double = 0.0
        opts.Itf double = 0.0
        opts.Vtr double = 0.0
        opts.Vtf double = 0.0
        opts.inSOC double = 0.0
        opts.chTime double = 0.0 % hours
        opts.alphaStep double = 0.01
        opts.verbose logical = true
    end

    %% Derived quantities
    Vm = sqrt(2) * Vrms;

    %% Log inputs (optional)
    if opts.verbose
        fprintf("Running FWR charging with:\n");
        fprintf(" Vrms = %.2f V\n f = %.1f Hz\n Battery = %.2f V\n", Vrms, f, E);
        fprintf(" R = %.3f ohm,  Capacity = %.2f Ah\n", R, capacity);
    end

    %% Find firing interval
    [alpha1, beta] = find_range_of_alphas_and_beta_with_v_drop(E, opts.Vforward, Vm);

    %% Alpha sweep
    alphas = make_alpha_array(alpha1, beta, opts.alphaStep);

    %% Find firing interval
    [alpha1_bat_charging, beta_bat_charging] = find_range_of_alphas_and_beta(E, Vm);

    %% Alpha sweep
    alphas_bat_charging = make_alpha_array(alpha1_bat_charging, beta_bat_charging, opts.alphaStep);

    %% Time axis
    x = linspace(0, 2*pi, 10000);

    %% PLOTS / SIMULATION
    if E > 0 && R > 0 && capacity > 0 && Vm > 0
        plot_time_of_charging(alphas_bat_charging, Vm, E, R, capacity, beta_bat_charging, true, 'log');
    end
    
    if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.inSOC >= 0 && opts.chTime > 0
        plot_final_soc(alphas_bat_charging, Vm, E, R, opts.inSOC, capacity, opts.chTime, beta_bat_charging, true);
    end

    % if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.Ileak && opts.Vforward && opts.Itr && opts.Itf
    %     plot_i_thyristor_with_t(x, alpha1, beta, Vm, E, R, opts.Ileak, opts.Vforward, opts.Itr, opts.Itf, f);
    % end
    % 
    % if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.Ileak && opts.Vforward && opts.Vtr && opts.Vtf
    %     plot_Vak_with_t(x, alpha1, beta, Vm, E, opts.Vforward, opts.Vtr, opts.Vtf, f);
    % end
    % 
    % if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.Ileak && opts.Vforward && opts.Itr && opts.Itf && opts.Vtr && opts.Vtf
    %     plot_power_loss_with_t(x, Vm, E, R, opts.Ileak, opts.Vforward, opts.Itr, opts.Itf, opts.Vtr, opts.Vtf, f, alpha1, beta);
    % end

    if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.Ileak && opts.Vforward && opts.Itr && opts.Itf && opts.Vtr && opts.Vtf
        plot_power_loss_with_alpha(alphas, beta, Vm, E, R, opts.Ileak, opts.Vforward, opts.Itr, opts.Itf, opts.Vtr, opts.Vtf, f);
    end

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

function val = i_battery_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
% Convert rise/fall times to angular ranges
omega = 2*pi*freq;
d_theta_r = omega * Itr;
d_theta_f = omega * Itf;

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

function val = i_thyristor_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
% same shape as battery current but simpler (thyristor 1)
omega = 2*pi*freq;
d_theta_r = omega * Itr;
d_theta_f = omega * Itf;
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

function val = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq)
omega = 2*pi*freq;
d_theta_r = omega * Vtr;
d_theta_f = omega * Vtf;

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

function plot_Vak_with_t(x, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq)
y = arrayfun(@(xi) Vak_including_drop(xi, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq), x);
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

function plot_i_thyristor_with_t(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
y = arrayfun(@(xi) i_thyristor_t_including_leakage(xi, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq), x);
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

function Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
f = @(xx) arrayfun(@(xval) i_battery_t_including_leakage(xval, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq).^2, xx);
integral_val = integral(f, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6, 'AbsTol',1e-12);
Irms = sqrt(integral_val / (2.0 * pi));
end

function pl = compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
Vak_thyrestor_1 = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq);
i_thyristor_1 = i_thyristor_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq);
resistance_loss = esr * (i_battery_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq).^2);
pl = 2*(Vak_thyrestor_1 .* i_thyristor_1) + resistance_loss;
end

function P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
f = @(xx) arrayfun(@(xval) compute_power_loss(xval, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq), xx);
result = integral(f, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6, 'AbsTol',1e-12);
P_loss_avg = (1.0 / (2.0 * pi)) * result;
end

function plot_power_loss_with_alpha(alpha_array, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
y = zeros(size(alpha_array));
for k = 1:numel(alpha_array)
    alpha = alpha_array(k);
    % Irms is computed but not used in original except maybe for reference
    Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq); %#ok<NASGU>
    P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq);
    y(k) = P_loss_avg;
end
figure; plot(alpha_array, y);
xlabel('Firing Angle (rad)'); ylabel('Average Power Loss (W)');
title('Average Power Loss vs Firing Angle'); grid on;
end

function plot_power_loss_with_t(x, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq, alpha, beta)
y = arrayfun(@(xi) compute_power_loss(xi, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq), x);
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

function plot_time_of_charging(alphas, Vm, E, esr, C, beta, degrees, charging_time_scale)
% Plot Iavg (left axis) and charging_time (right axis).
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

title(sprintf('Iavg and charging time vs firing angle (Vpeak=%g V, E=%g V, esr=%g ohms, C=%g Ah)', Vm, E, esr, C), 'Interpreter','none');
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
