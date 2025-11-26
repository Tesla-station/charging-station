function half_wave_rectifier(Vrms, f, E, R, capacity, opts)
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
        plot_time_of_charging(alphas_bat_charging, beta_bat_charging, Vm, E, R, capacity, true, 'log');
    end
    
    if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.inSOC >= 0 && opts.chTime > 0
        plot_final_soc(alphas_bat_charging, beta_bat_charging, Vm, E, R, opts.inSOC, capacity, opts.chTime, true);
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
    %     plot_power_loss_with_t(alpha1, beta, x, Vm, E, R, opts.Ileak, opts.Vforward, opts.Itr, opts.Itf, opts.Vtr, opts.Vtf, f);
    % end

    if E > 0 && R > 0 && capacity > 0 && Vm > 0 && opts.Ileak && opts.Vforward && opts.Itr && opts.Itf && opts.Vtr && opts.Vtf
        plot_power_loss_with_alpha(alphas, beta, Vm, E, R, opts.Ileak, opts.Vforward, opts.Itr, opts.Itf, opts.Vtr, opts.Vtf, f);
    end
end

%% --- Utility functions ---------------------------------------------------

function out = norm_angle(theta)
    % Normalize angle to [0, 2*pi)
    out = mod(theta, 2*pi);
end

function [t1, t2] = find_range_of_alphas_and_beta(E, Vm)
    frac = E / Vm;
    if abs(frac) > 1
        error('|E/Vm| must be <= 1. Got %g', frac);
    end
    a = asin(frac);
    t1 = norm_angle(a);
    t2 = norm_angle(pi - a);
end

function [t1, t2] = find_range_of_alphas_and_beta_with_v_drop(E, v_drop, Vm)
    frac = (E + v_drop) / Vm;
    if abs(frac) > 1
        error('|E/Vm| must be <= 1. Got %g', frac);
    end
    a = asin(frac);
    t1 = norm_angle(a);
    t2 = norm_angle(pi - a);
end

function arr = make_alpha_array(a1, a2, step)
    % Create array from a1 to a2 (exclusive) with step, handling wrap-around
    if nargin < 3
        step = 0.01;
    end
    if step <= 0
        error('step must be positive');
    end
    if abs(a1 - a2) < eps
        arr = a1;
        return;
    end
    two_pi = 2*pi;
    if a1 < a2
        % mimic numpy.arange(a1, a2, step)
        arr = a1:step:(a2-step);
        if isempty(arr)
            arr = a1; % fallback minimal array
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
    y = Vm .* sin(x);
end

function i = i_t(x, alpha, beta, Vm, E, esr)
    % scalar-style function (use arrayfun for vector inputs)
    if (alpha <= x) && (x <= beta)
        i = (vs_t(x, Vm) - E) / esr;
    else
        i = 0.0;
    end
end

function i = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop)
    if (alpha <= x) && (x <= beta)
        i = (vs_t(x, Vm) - E - v_drop) / esr;
    else
        i = -ileak;
    end
end

function i = i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
    % Scalar x expected. For vector x, call with arrayfun.
    omega = 2*pi*freq;
    d_theta_r = omega * Itr;
    d_theta_f = omega * Itf;
    
    naturalI = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop);
    
    if x < alpha
        i = -ileak;
        return;
    end
    if (alpha <= x) && (x < alpha + d_theta_r)
        [m,b] = line_from_points(alpha, -ileak, alpha + d_theta_r, i_t_including_leakage_simplified(alpha + d_theta_r, alpha, beta, Vm, E, esr, ileak, v_drop));
        i = m*x + b;
        return;
    end
    if (alpha + d_theta_r <= x) && (x < beta)
        i = naturalI;
        return;
    end
    if (beta <= x) && (x < beta + d_theta_f)
        [m,b] = line_from_points(beta, i_t_including_leakage_simplified(beta, alpha, beta, Vm, E, esr, ileak, v_drop), beta + d_theta_f, -ileak);
        i = m*x + b;
        return;
    end
    i = -ileak;
end

function v = Vak_including_drop_simplified(x, alpha, beta, Vm, E, v_drop)
    if (alpha <= x) && (x <= beta)
        v = v_drop;
    else
        v = Vm * sin(x) - E;
    end
end

function plot_Vak_with_t_simplified(x, alpha, beta, Vm, E, v_drop)
    y = arrayfun(@(xi) Vak_including_drop_simplified(xi, alpha, beta, Vm, E, v_drop), x);
    plot(x, y);
    xlabel('Angle (ωt) rad');
    ylabel('Vak (V)');
    title('Voltage across the load Vak vs Angle (ωt) rad');
    grid on;
    figure;
end

function [m, b] = line_from_points(x1, y1, x2, y2)
    m = (y2 - y1) / (x2 - x1);
    b = y1 - m * x1;
end

function v = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq)
    omega = 2*pi*freq;
    d_theta_r = omega * Vtr;
    d_theta_f = omega * Vtf;
    natural_v = vs_t(x, Vm) - E;
    
    if x < alpha
        v = natural_v;
        return;
    end
    if (alpha <= x) && (x < alpha + d_theta_f)
        [m,b] = line_from_points(alpha, vs_t(alpha, Vm)-E, alpha + d_theta_f, v_drop);
        v = m*x + b;
        return;
    end
    if (alpha + d_theta_f <= x) && (x < beta)
        v = v_drop;
        return;
    end
    if (beta <= x) && (x < beta + d_theta_r)
        [m,b] = line_from_points(beta, v_drop, beta + d_theta_r, vs_t(beta + d_theta_r, Vm)-E);
        v = m*x + b;
        return; 
    end
    v = natural_v;
end

function plot_Vak_with_t(x, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq)
    y = arrayfun(@(xi) Vak_including_drop(xi, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq), x);
    plot(x, y);
    xlabel('Angle (ωt) rad');
    ylabel('Vak (V)');
    title('Voltage across the load Vak vs Angle (ωt) rad');
    grid on;
    figure; 
end

function plot_i_thyristor_with_t(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
    y = arrayfun(@(xi) i_t_including_leakage(xi, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq), x);
    plot(x, y);
    xlabel('Angle (ωt) rad');
    ylabel('Current i (A)');
    title('Current i vs Angle (ωt) rad');
    grid on;
    figure;
end

function Iavg = compute_Iavg(Vm, E, esr, alpha, beta)
    integrand = @(xx) arrayfun(@(xval) i_t(xval, alpha, beta, Vm, E, esr), xx);
    result = integral(integrand, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-9);
    Iavg = (1.0 / (2.0 * pi)) * result;
end

function Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq)
    integrand = @(xx) arrayfun(@(xval) i_t_including_leakage(xval, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq).^2, xx);
    integral_val = integral(integrand, 0, 2.0*pi, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-9);
    Irms = sqrt(integral_val / (2.0 * pi));
end

function P = compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
    Vak = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, Vtr, Vtf, freq);
    iVal = i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq);
    P = Vak .* iVal + esr .* (iVal.^2);
end

function P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
    power_loss_func = @(xx) arrayfun(@(xval) compute_power_loss(xval, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq), xx);
    result = integral(power_loss_func, 0, 2*pi, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-9);
    P_loss_avg = (1.0 / (2.0 * pi)) * result;
end

function plot_power_loss_with_alpha(alpha_array, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
    y = zeros(size(alpha_array));
    for k = 1:numel(alpha_array)
        alpha = alpha_array(k);
        % Irms computed but unused (kept like original)
        Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, freq);
        P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq);
        y(k) = P_loss_avg;
    end
    plot(alpha_array, y);
    xlabel('Firing Angle (rad)');
    ylabel('Average Power Loss (W)');
    title('Average Power Loss vs Firing Angle');
    grid on;
    figure;
end

function plot_power_loss_with_t(alpha, beta, x, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq)
    y = arrayfun(@(xi) compute_power_loss(xi, alpha, beta, Vm, E, esr, ileak, v_drop, Itr, Itf, Vtr, Vtf, freq), x);
    plot(x, y);
    xlabel('Angle (ωt) rad');
    ylabel('Power Loss (W)');
    title('Power Loss vs Angle');
    grid on;
    figure;
end

function charging_time = compute_charging_time(Iavg, C)
    % Compute charging time (hours) from capacity C (Ah) and Iavg (A)
    charging_time = nan(size(Iavg));
    if ~isempty(Iavg)
        positive = Iavg > 0;
        if ~isempty(C)
            charging_time(positive) = C ./ Iavg(positive);
        end
    end
end

function plot_time_of_charging(alphas, beta, Vm, E, esr, C, degreesFlag, charging_time_scale)
    if nargin < 7
        degreesFlag = true;
    end
    if nargin < 8
        charging_time_scale = 'linear';
    end
    if degreesFlag
        x = rad2deg(alphas);
        xlabelTxt = 'Firing angle (deg)';
    else
        x = alphas;
        xlabelTxt = 'Firing angle (rad)';
    end
    [Iavg, charging_time] = time_of_charging(Vm, E, esr, C, alphas, beta);
    
    figure('Units','normalized','Position',[0.1 0.1 0.6 0.6]);
    yyaxis left
    plot(x, Iavg, '-o', 'MarkerSize', 4);
    ylabel('Iavg (A)');
    xlabel(xlabelTxt, 'FontSize', 11);
    grid on;
    set(gca,'FontSize',11);
    
    if ~all(isnan(charging_time))
        yyaxis right
        plot(x, charging_time, '-s', 'MarkerSize', 4);
        ylabelText = 'Charging time (h)';
        if exist('charging_time_scale','var') && ~strcmp(charging_time_scale,'linear')
            set(gca, 'YScale', charging_time_scale);
            ylabelText = sprintf('%s (%s)', ylabelText, charging_time_scale);
        end
        ylabel(ylabelText);
    end
    
    title(sprintf('Iavg and charging time vs firing angle (Vrms=%g V, E=%g V, esr=%g ohms, C=%g Ah)', Vm/sqrt(2), E, esr, C));
    figure;
end

function new_soc = compute_final_soc(Iavg, initial_soc, C, required_charging_time)
    % initial_soc in percent (0-100)
    new_soc = (Iavg .* required_charging_time + (initial_soc/100.0) * C) ./ C;
    new_soc = new_soc * 100.0;
    new_soc = min(max(new_soc, 0.0), 100.0);
end

function plot_final_soc(alphas, beta, Vm, E, esr, initial_soc, C, required_charging_time, degreesFlag)
    if nargin < 9
        degreesFlag = true;
    end
    if degreesFlag
        x = rad2deg(alphas);
        xlabelTxt = 'Firing angle (deg)';
    else
        x = alphas;
        xlabelTxt = 'Firing angle (rad)';
    end
    
    [Iavg, ~] = time_of_charging(Vm, E, esr, C, alphas, beta);
    new_soc = compute_final_soc(Iavg, initial_soc, C, required_charging_time);
    
    figure;
    plot(x, new_soc, '-o', 'MarkerSize', 4);
    xlabel(xlabelTxt, 'FontSize', 11);
    ylabel('Final State of Charge (%)', 'FontSize', 11);
    grid on;
    title(sprintf('Final SOC vs Firing Angle\nInit SOC=%g%%, Time=%gh, C=%gAh', initial_soc, required_charging_time, C));
    ylim([0 100]);
    figure;
end

function [Iavg, charging_time] = time_of_charging(Vm, E, esr, C, alphas, beta)
    Iavg = zeros(size(alphas));
    for k = 1:numel(alphas)
        Iavg(k) = compute_Iavg(Vm, E, esr, alphas(k), beta);
    end
    charging_time = compute_charging_time(Iavg, C);
end
