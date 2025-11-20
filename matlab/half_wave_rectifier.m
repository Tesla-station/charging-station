function half_wave_rectifier()
    % HALF_WAVE_RECTIFIER_V2 MATLAB port of your Python script
    % Save this as half_wave_rectifier_v2.m and run `main`
    %
    % Author: converted from user's Python script

    % === USER PARAMETERS ===
    Vrms = 60;
    E = 12;
    esr = 4.256;
    capacity = 100;
    v_drop = 1.2;
    tr = 10e-6;        % rise time (s)
    tf = 5e-6;         % fall time (s)
    freq = 50;         % Hz
    ileak = 0.02;      % leakage current (A)

    Vm = sqrt(2) * Vrms;

    % compute alpha/beta using v_drop aware function
    [alpha1, beta] = find_range_of_alphas_and_beta_with_v_drop(E, v_drop, Vm);

    alphas = make_alpha_array(alpha1, beta, 0.01);

    x = linspace(0, 2*pi, 10000);

    % plots & calculations (uncomment what you want)
    plot_Vak_with_t(x, alpha1, beta, Vm, E, v_drop, tr, tf, freq);
    plot_i_with_t(x, alpha1, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
    plot_power_loss_with_alpha(alphas, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);

    % example: compute Iavg array
    %[alphas_out, Iavg, Tchg] = time_of_charging(Vrms, E, esr, capacity, 0.01, true, true, 'linear', 0, 0.0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility / basic helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t = norm_angle(theta)
    % Normalize angle to [0, 2*pi)
    t = mod(theta, 2*pi);
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
    % uses E + v_drop in numerator as in your python version
    frac = (E + v_drop) / Vm;
    if abs(frac) > 1
        error('|(E+v_drop)/Vm| must be <= 1. Got %g', frac);
    end
    a = asin(frac);
    t1 = norm_angle(a);
    t2 = norm_angle(pi - a);
end

function arr = make_alpha_array(a1, a2, step)
    if step <= 0
        error('step must be positive');
    end
    if abs(a1 - a2) < 1e-12
        arr = a1;
        return;
    end
    two_pi = 2*pi;
    if a1 < a2
        arr = a1:step:(a2-step);
    else
        part1 = a1:step:(two_pi-step);
        part2 = 0:step:(a2-step);
        arr = [part1, part2];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = vs_t(x, Vm)
    v = Vm .* sin(x);
end

function I = i_t(x, alpha, beta, Vm, E, esr)
    if (x >= alpha) && (x <= beta)
        I = (Vm * sin(x) - E) / esr;
    else
        I = 0.0;
    end
end

function I = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop)
    % Here python used Vm*sin - E - v_drop in simplified version
    if (x >= alpha) && (x <= beta)
        I = (Vm * sin(x) - E - v_drop) / esr;
    else
        I = ileak;
    end
end

function p = line_from_points(x1, y1, x2, y2)
    m = (y2 - y1) / (x2 - x1);
    b = y1 - m * x1;
    p = [m, b];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i(t) piecewise with leakage and switching slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Iout = i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    % Convert rise/fall times into angular widths
    omega = 2*pi*freq;
    d_theta_r = omega * tr;   % angle width for rise (after alpha)
    d_theta_f = omega * tf;   % angle width for fall (after beta)

    % natural conduction expression (mirrors simplified version)
    natural_I = i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak, v_drop);

    % Region 1: before alpha -> leakage
    if x < alpha
        Iout = ileak;
        return;
    end

    % Region 2: alpha -> alpha + d_theta_r : ramp from ileak to conduction
    if (x >= alpha) && (x < alpha + d_theta_r)
        target_angle = min(alpha + d_theta_r, beta); % avoid exceeding beta
        target_I = i_t_including_leakage_simplified(target_angle, alpha, beta, Vm, E, esr, ileak, v_drop);
        p = line_from_points(alpha, ileak, alpha + d_theta_r, target_I);
        Iout = p(1)*x + p(2);
        return;
    end

    % Region 3: conduction (natural current)
    if (x >= alpha + d_theta_r) && (x < beta)
        Iout = natural_I;
        return;
    end

    % Region 4: beta -> beta + d_theta_f : ramp from conduction back to leakage
    if (x >= beta) && (x < beta + d_theta_f)
        start_I = i_t_including_leakage_simplified(beta, alpha, beta, Vm, E, esr, ileak, v_drop);
        p = line_from_points(beta, start_I, beta + d_theta_f, ileak);
        Iout = p(1)*x + p(2);
        return;
    end

    % Region 5: after beta + d_theta_f -> leakage
    Iout = ileak;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V_AK piecewise with switching slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = Vak_including_drop_simplified(x, alpha, beta, Vm, E, v_drop)
    if (x >= alpha) && (x <= beta)
        v = v_drop;
    else
        v = Vm * sin(x) - E;
    end
end

function v = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, tr, tf, freq)
    omega = 2*pi*freq;
    d_theta_r = omega * tr;   % angle for rise (after beta)
    d_theta_f = omega * tf;   % angle for fall (after alpha)

    % Note: your last python version used (Vs - E - v_drop) as natural_v.
    % We'll replicate that to stay faithful:
    natural_v = vs_t(x, Vm) - E;

    % Region 1: before alpha -> natural waveform
    if x < alpha
        v = natural_v;
        return;
    end

    % Region 2: alpha -> alpha + d_theta_f : linear fall to v_drop
    if (x >= alpha) && (x < alpha + d_theta_f)
        v_start = vs_t(alpha, Vm) - E;
        p = line_from_points(alpha, v_start, alpha + d_theta_f, v_drop);
        v = p(1)*x + p(2);
        return;
    end

    % Region 3: alpha + d_theta_f <= x < beta : constant v_drop
    if (x >= alpha + d_theta_f) && (x < beta)
        v = v_drop;
        return;
    end

    % Region 4: beta -> beta + d_theta_r : linear rise back to natural
    if (x >= beta) && (x < beta + d_theta_r)
        v_end = vs_t(beta + d_theta_r, Vm) - E;
        p = line_from_points(beta, v_drop, beta + d_theta_r, v_end);
        v = p(1)*x + p(2);
        return;
    end

    % Region 5: after beta + d_theta_r -> natural waveform
    v = natural_v;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Vak_with_t_simplified(x, alpha, beta, Vm, E, v_drop)
    y = arrayfun(@(xi) Vak_including_drop_simplified(xi, alpha, beta, Vm, E, v_drop), x);
    plot(x, y, 'LineWidth', 1.4);
    xlabel('Angle (rad)'), ylabel('Vak (V)');
    title('Vak simplified');
    grid on;
end

function plot_Vak_with_t(x, alpha, beta, Vm, E, v_drop, tr, tf, freq)
    y = arrayfun(@(xi) Vak_including_drop(xi, alpha, beta, Vm, E, v_drop, tr, tf, freq), x);
    plot(x, y, 'LineWidth', 1.4);
    xlabel('Angle (rad)'), ylabel('Vak (V)');
    title('Vak with switching slopes');
    grid on;
end

function plot_i_with_t(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    y = arrayfun(@(xi) i_t_including_leakage(xi, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq), x);
    plot(x, y, 'LineWidth', 1.4);
    xlabel('Angle (rad)'), ylabel('i(t) (A)');
    title('Current i vs Angle');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration / RMS / power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Iavg = compute_Iavg(Vm, E, esr, alpha, beta)
    f = @(x) i_t(x, alpha, beta, Vm, E, esr);
    % integrate over full cycle (function returns 0 outside conduction)
    res = integral(@(t) f(t), 0, 2*pi, 'ArrayValued', true);
    Iavg = res / (2*pi);
end

function Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    f2 = @(x) (i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)).^2;
    res = integral(@(t) f2(t), 0, 2*pi, 'ArrayValued', true);
    Irms = sqrt(res / (2*pi));
end

function p = compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    Vak = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, tr, tf, freq);
    i   = i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
    p = Vak .* i + esr .* (i.^2);
end

function Pavg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    fun = @(x) compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
    res = integral(@(t) fun(t), 0, 2*pi, 'ArrayValued', true);
    Pavg = res / (2*pi);
end

function plot_power_loss_with_alpha(alpha_array, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    y = zeros(size(alpha_array));
    for k = 1:numel(alpha_array)
        alpha = alpha_array(k);
        % compute Irms (not used inside compute_average_power_loss, but kept for parity)
        Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
        y(k) = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq);
    end
    plot(alpha_array, y, 'LineWidth', 1.4);
    xlabel('Firing Angle (rad)'), ylabel('Average Power Loss (W)');
    title('Average Power Loss vs Firing Angle');
    grid on;
end

function plot_power_loss_with_t(x, Vm, E, esr, ileak, tr, tf, freq, alpha, beta)
    y = arrayfun(@(xi) compute_power_loss(xi, alpha, beta, Vm, E, esr, ileak, 0, tr, tf, freq), x);
    plot(x, y, 'LineWidth', 1.4);
    xlabel('Angle (rad)'), ylabel('Power Loss (W)');
    title('Power Loss vs Angle');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time of charging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alphas, Iavg_arr, charging_time] = time_of_charging(Vrms, E, esr, C, step, degrees, showPlot, charging_time_scale, initial_soc, required_charging_time)
    if nargin < 5, step = 0.01; end
    if nargin < 6, degrees = true; end
    if nargin < 7, showPlot = true; end
    if nargin < 8, charging_time_scale = 'linear'; end
    if nargin < 9, initial_soc = 0.0; end
    if nargin < 10, required_charging_time = 0.0; end

    Vm = sqrt(2)*Vrms;
    [alpha1, beta] = find_range_of_alphas_and_beta(E, Vm);
    alphas = make_alpha_array(alpha1, beta, step);

    Iavg_arr = zeros(size(alphas));
    for k = 1:numel(alphas)
        Iavg_arr(k) = compute_Iavg(Vm, E, esr, alphas(k), beta);
    end

    % charging time C (Ah) / Iavg (A) -> hours? In your python you left as seconds, here we keep same units: C in Ah, I in A -> hours
    charging_time = nan(size(Iavg_arr));
    pos = Iavg_arr > 0;
    charging_time(pos) = C ./ Iavg_arr(pos);

    % Plot
    if degrees
        x = rad2deg(alphas);
        xlabelStr = 'Firing angle (deg)';
    else
        x = alphas;
        xlabelStr = 'Firing angle (rad)';
    end

    figure;
    yyaxis left;
    plot(x, Iavg_arr, '-o'); grid on;
    ylabel('Iavg (A)');

    yyaxis right;
    plot(x, charging_time, '-s');
    ylabel('Charging time (h)'); % interpreted as hours because C (Ah) / I (A) -> hours
    if ~strcmp(charging_time_scale, 'linear')
        set(gca, 'YScale', charging_time_scale);
    end

    xlabel(xlabelStr);
    title(sprintf('Iavg and charging time vs firing angle (Vrms=%g, E=%g, esr=%g, C=%gAh)', Vrms, E, esr, C));

    if required_charging_time > 0 && initial_soc > 0
        new_soc = (Iavg_arr * required_charging_time + (initial_soc/100.0).*C) ./ C;
        new_soc = new_soc * 100.0;
        new_soc = min(max(new_soc, 0), 100);

        figure;
        plot(x, new_soc, '-o'); grid on;
        xlabel(xlabelStr);
        ylabel('Final SOC (%)');
        title(sprintf('Final SOC vs Firing Angle (Init SOC=%g%%, Time=%gh)', initial_soc, required_charging_time));
    end
end
