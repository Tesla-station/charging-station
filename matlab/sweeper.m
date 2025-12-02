alphas = 0:10:180;     % Sweep from 0° to 180° in steps of 5°
Vmean_half_wave_values = zeros(size(alphas));
Vmean_full_wave_values = zeros(size(alphas));
Vmean_full_wave_ct_values = zeros(size(alphas));

for i = 1:length(alphas)
    alpha_deg = alphas(i);   % Update firing angle
    % Load the parameter file
    run('simulink_params.m')  

    % Run simulation
    simOut_half_wave = sim('half_wave');   % half wave
    simOut_full_wave = sim('full_wave');   % full wave
    simOut_full_wave_ct = sim('full_wave_center_tapped');
    % Read mean voltage
    Vmean_half_wave_values(i) = simOut_half_wave.Vmean(end);
    Vmean_full_wave_values(i) = simOut_full_wave.Vmean(end);
    Vmean_full_wave_ct_values(i) = simOut_full_wave_ct.Vmean(end);
end

% Plot result
figure;
plot(alphas, Vmean_half_wave_values, 'LineWidth', 2);
xlabel('Firing Angle α (degrees)');
ylabel('Mean Output Voltage (V)');
title('Mean Voltage vs. Firing Angle — Half wave Rectifier');
grid on;

figure;
plot(alphas, Vmean_full_wave_values, 'LineWidth', 2);
xlabel('Firing Angle α (degrees)');
ylabel('Mean Output Voltage (V)');
title('Mean Voltage vs. Firing Angle — Full wave Rectifier');
grid on;

figure;
plot(alphas, Vmean_full_wave_ct_values, 'LineWidth', 2);
xlabel('Firing Angle α (degrees)');
ylabel('Mean Output Voltage (V)');
title('Mean Voltage vs. Firing Angle — Center tapped Rectifier');
grid on;