f = 50; % Hz
alpha_deg; % To change the firing angle you only need to change this value, leave it empty for the sweeper.
alpha_2_deg = 180 + alpha_deg;
Vin_rms = 230; % Volt
Vin_peak = sqrt(2)*Vin_rms;

alpha_rad = (alpha_deg*pi)/(180.0);
alpha_2_rad = (alpha_2_deg*pi)/(180.0);

alpha_1_t = alpha_rad/(2*pi*f);
alpha_2_t = alpha_2_rad/(2*pi*f);

R = 20; % Ohm
L = 0.1; % Henry
