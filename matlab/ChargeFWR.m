function ChargeFWR(Vrms,f,E,R,Ecap,EcapUnit,opts)
% CHARGEFWR  Simulates a fully-controlled full-wave rectifier (FWR)
%           charging a battery and evaluates its electrical performance.
%
%   CHARGEFWR(Vrms, f, E, R, Ecap, EcapUnit, opts) computes charging
%   behavior, average charging current, charge time, final state of charge,
%   and thyristor switching losses for a full-wave rectifier charging a
%   battery through a series resistance.
%
%   INPUT ARGUMENTS
%   ----------------
%   Vrms       - Source RMS voltage (V)
%   f          - Source frequency (Hz)
%   E          - Battery voltage (V)
%   R          - Series resistance (ohms)
%   Ecap       - Battery capacity (Wh or Ah depending on EcapUnit)
%   EcapUnit   - "Wh" or "Ah"
%
%   opts.chargeTime   - Charging time interval (hours)
%   opts.SoC          - Initial state of charge (%) before charging
%   opts.Vforward     - Thyristor forward voltage drop (V)
%   opts.Ileak        - Thyristor leakage current (A)
%   opts.triseV       - Voltage rise time for thyristor turn-on (s)
%   opts.tfallV       - Voltage fall time for thyristor turn-off (s)
%   opts.triseI       - Current rise time for thyristor turn-on (s)
%   opts.tfallI       - Current fall time for thyristor turn-off (s)
%
%
%   FUNCTION DESCRIPTION
%   --------------------
%   The function performs the following steps:
%
%   1. Computes the firing interval where the rectifier conducts by solving
%      Vs(t) = E + 2*Vforward.
%
%   2. Computes average charging current for all possible firing angles.
%
%   3. Computes charging time for each firing angle and plots charge time
%      vs. firing angle.
%
%   4. Computes final battery state of charge (SoC) after a specified
%      chargeTime, and plots SoC vs. firing angle.
%
%   5. Calculates detailed thyristor switching losses using a 9-region
%      piecewise model. Computes:
%         - instantaneous voltage and current waveforms for T1–T4
%         - instantaneous power loss p(t)
%         - RMS load current
%         - average switching and conduction loss per electrical cycle
%
%   6. Plots total power loss vs. firing angle.
%
%
%   OUTPUT
%   ------
%   The function does not return a variable directly but creates plots and
%   prints analysis results. Detailed loss values are computed internally
%   and plotted for inspection.
%   The power losses graph is optional and will only be created if the
%   values for triseI, tfallI, triseV, tfallV are provided
%
%
%   EXAMPLE
%   -------
%   ChargeFWR(110, 100, 12, 0.25, 400, "Wh", ...
%             Vforward=0.7, chargeTime=0.5, SoC=10, ...
%             tfallV=80e-6, triseV=50e-6, tfallI=80e-6, triseI=50e-6, ...
%             Ileak=0.05)
    arguments (Input)
        Vrms double
        f int32
        E double
        R double
        Ecap double
        EcapUnit string = "Wh";
        opts.chargeTime double = 0
        opts.SoC double = 0
        opts.Vforward double = 0
        opts.Ileak double = 0
        opts.triseV double = 0
        opts.tfallV double = 0
        opts.triseI double = 0
        opts.tfallI double = 0
    end
    f = double(f);
    T = 1/f;
    opts.Ileak = - opts.Ileak;
    Vrms = double(Vrms);
    
    % Check for valid input,E and Vforward values
    if(sqrt(2) * Vrms < E + 2 * opts.Vforward)
        disp("Source not strong enough to charge the Batter please select a source with Vrms > " + (E + 2* opts.Vforward)/sqrt(2))
        return
    elseif(R <= 0)
        disp("Please provide a valid value of R ex. R>0")
        return
    end

    % get exact roots
    [minFiringTime, maxFiringTime] = calcMinMaxFiringTimes(Vrms,E,f,opts.Vforward);
    
    % Source equation
    Vs = @(t) sqrt(2) * Vrms * sin(2*pi*f*t);

   % Time sampling array for firing 
    firingSamplingRate = 500;
    samplingRate = 500;
    avgCurrents = calcAvgCurrentValues(minFiringTime,maxFiringTime,firingSamplingRate,Vs,opts.Vforward,R,f,E,"FWR");

    % Plot avgCurrent w.r.t firingtime REMOVE NOT NEEDED
    %anglePlotter(minFiringTime,maxFiringTime,firingSamplingRate,avgCurrents,"Firing Angle","avgCurrent",f)

    % Calculate Charge time w.r.t firing angle
    chargeTimes = calcChargeTimes(E,avgCurrents,Ecap,EcapUnit);
    % Plot Charge time w.r.t firing angle
    anglePlotter(minFiringTime,maxFiringTime,firingSamplingRate,chargeTimes,"Firing Angle","chargeTime",f);

    % Calculate SoC w.r.t firing angle
    if (opts.chargeTime ~= 0)
        finalSoC = calcFinalSoC(opts.chargeTime,Ecap,avgCurrents,E,opts.SoC,EcapUnit);
        anglePlotter(minFiringTime,maxFiringTime,firingSamplingRate,finalSoC,"Firing Angle","finalSoC",f)
    end

    % Losses Calculations (9 regions for the thyristor)
    if(opts.tfallI == 0 || opts.tfallV == 0 || opts.triseI == 0 || opts.triseV == 0)
        
    else
        firingTimeSamples = linspace(minFiringTime,maxFiringTime,firingSamplingRate);
        powerLoss = zeros(1,firingSamplingRate);
        for i = 1:length(firingTimeSamples)
            powerLoss(i) = calcPowerLoss(firingTimeSamples(i), opts.tfallV, opts.triseV, opts.tfallI, opts.triseI, maxFiringTime, opts.Vforward, opts.Ileak, Vs, E, R, T, samplingRate);
        end
        anglePlotter(minFiringTime,maxFiringTime,samplingRate,powerLoss,"firing time","powerLoss",f)
    end

end

% -------------------------------------------------------------------------
%---------------------------HELPER FUNCTIONS-------------------------------
%--------------------------------------------------------------------------

% Calculates min and max firing times
function [minFiringTime, maxFiringTime] = calcMinMaxFiringTimes(Vrms,E,f,vForward)
    minFiringTime = (1/(2*pi*f))*asin((2*vForward+E)/(sqrt(2)*Vrms));
    maxFiringTime = (1/(2*pi*f))*(pi - asin((2*vForward+E)/(sqrt(2)*Vrms)));
end


% Calculates avg currents w.r.t firing angle range
function avgCurrents = calcAvgCurrentValues(minFiringTime,maxFiringTime,firingSamplingRate,Vs,vForward,R,f,E,circuitType)
    if(circuitType == "FWR")
        n = 2;
    else
        n = 1;
    end
    % Time sampling array for firing 
    firingSamples = linspace(minFiringTime,maxFiringTime,firingSamplingRate);
    % Calculating Average Current values for each firing angle
    avgCurrents = zeros(1,firingSamplingRate);
    I =@(t) (Vs(t) - E - 2*vForward)/R;
    for i = 1:firingSamplingRate
        avgCurrents(i) =n * f * integral(I,firingSamples(i),maxFiringTime);
    end
end

% Plots any function w.r.t angle
function anglePlotter(lowerLimit,upperLimit,samplingRate,data,xLabel,yLabel,f)
    lowerLimit = 2*f*lowerLimit*180;
    upperLimit = 2*f*upperLimit*180;
    samples = linspace(lowerLimit,upperLimit,samplingRate);
    figure
    plot(samples,data)
    xlabel(xLabel);
    ylabel(yLabel);
    angles = xticks;                            % get tick values
    labels = strcat(string(angles), "°");       % append degree
    xticklabels(labels);
end

%Calculates Charge time
function chargingTimes = calcChargeTimes(E,avgCurrents,Ecap,EcapUnit)
    if(EcapUnit == "Wh")
        Pavg = E .* avgCurrents;
        chargingTimes = Ecap ./ Pavg;
    elseif (EcapUnit == "Ah")
        chargingTimes = Ecap ./ avgCurrents;
    else    
        disp("Please use Wh or Ah as the input unit and try again");
        chargingTimes = zeros(size(avgCurrents));
    end
end


%calculate finalSoC
function finalSoC = calcFinalSoC(chargeTime,Ecap,avgCurrents,E,SoC,EcapUnit)
    SoC = SoC/100;

    if(EcapUnit == "Wh")
        Pavg = E .* avgCurrents;
        finalSoC = SoC + ((chargeTime .* Pavg) ./ Ecap);
        finalSoC(finalSoC > 1) = 1;
    elseif(EcapUnit == "Ah")
        finalSoC = SoC + ((chargeTime .* avgCurrents) ./ Ecap);
        finalSoC(finalSoC > 1) = 1;
    else
        finalSoC = zeros(size(avgCurrents));
        disp("Please use Wh or Ah as the input of the capacity type");
    end
end

% Linear function generator
function linearFunc = generateLinearFunc(startX,deltaX,startY,endY)
    linearFunc = @(t) startY + (endY - startY) * (t - startX) / (deltaX);
end


% Calculates rms Current value
function rms = calcRMS(x,y,f)
    rms = sqrt(f*trapz(x,y.^2));
end





% Generates piecewise equations for V and I and calculates the avg power
% loss for a given firingTime
function powerloss = calcPowerLoss(firingTime,tfallV,triseV,tfallI,triseI,maxFiringTime,Vforward,Ileak,Vs,E,R,T,samplingRate)

    % Genrating the Linear regions for I and V instantaneous
    %VfallT1,3
    VfallT1 = generateLinearFunc(firingTime,tfallV,(Vs(firingTime)-E)/2,Vforward);

    %VriseT1,3
    VriseT1 = generateLinearFunc(maxFiringTime,triseV,Vforward,(Vs(maxFiringTime + triseV) - E)/2);

    %IriseT1,3
    IriseT1 = generateLinearFunc(firingTime,triseI,Ileak,(Vs(firingTime + triseI) - E - 2 * Vforward)/R);

    %IfallT1,3
    IfallT1 = generateLinearFunc(maxFiringTime,tfallI,(Vs(maxFiringTime) - E - 2 * Vforward)/R,Ileak);


    %VfallT2,4
    VfallT2 = generateLinearFunc((T/2) + firingTime,tfallV,(-Vs((T/2) + firingTime)-E)/2,Vforward);

    %VriseT2,4
    VriseT2 = generateLinearFunc((T/2) + maxFiringTime,triseV,Vforward,(-Vs(maxFiringTime + triseV) - E)/2);

    %IriseT2,4
    IriseT2 = generateLinearFunc(T/2 + firingTime,triseI,Ileak,(-Vs((T/2) + firingTime + triseI) - E - 2 * Vforward)/R);

    %IfallT2,4
    IfallT2 = generateLinearFunc(T/2 + maxFiringTime,tfallI,(-Vs((T/2) + maxFiringTime) - E - 2 * Vforward)/R,Ileak);



    % Constructing full instantaneous V function for T1,3
    VT1 = @(t) ...
    (t < firingTime)              .* (Vs(t)-E)/2 + ...
    (t >= firingTime & t < firingTime + tfallV)    .* VfallT1(t) + ...
    (t >= firingTime + tfallV & t < maxFiringTime) .* Vforward + ...
    (t >= maxFiringTime & t < maxFiringTime + triseV) .* VriseT1(t) + ...
    (t >= maxFiringTime + triseV & t < ((T/2) + firingTime)) .* (Vs(t)-E)/2 + ...
    (t >= ((T/2) + firingTime) & t < ((T/2) + firingTime + tfallV)) .* (Vs(t) + VfallT2(t)) + ...
    (t >= ((T/2) + firingTime + tfallV) & t < ((T/2) + maxFiringTime)) .* (Vs(t) + Vforward) +...
    (t >= ((T/2) + maxFiringTime) & t < ((T/2) + maxFiringTime + triseV)) .* (Vs(t) + VriseT2(t)) + ...
    (t >= ((T/2) + maxFiringTime + triseV)) .* (Vs(t)-E)/2;

    % Constructing full instantaneous V function for T2,4
    VT2 = @(t) ...
    (t < firingTime)              .* (-Vs(t)-E)/2 + ...
    (t >= firingTime & t < firingTime + tfallV)    .* (VfallT1(t) - Vs(t)) + ...
    (t >= firingTime + tfallV & t < maxFiringTime) .* (-Vs(t) + Vforward) + ...
    (t >= maxFiringTime & t < maxFiringTime + triseV) .* (VriseT1(t)-Vs(t)) + ...
    (t >= maxFiringTime + triseV & t < ((T/2) + firingTime)) .* (-Vs(t)-E)/2 + ...
    (t >= ((T/2) + firingTime) & t < ((T/2) + firingTime + tfallV)) .* VfallT2(t) + ...
    (t >= ((T/2) + firingTime + tfallV) & t < ((T/2) + maxFiringTime)) .* Vforward +...
    (t >= ((T/2) + maxFiringTime) & t < ((T/2) + maxFiringTime + triseV)) .* VriseT2(t) + ...
    (t >= ((T/2) + maxFiringTime + triseV)) .* (-Vs(t)-E)/2;


    % Constructing full instantaneous I function for T1,3
    IT1 = @(t) ...
    (t < firingTime)              .* Ileak + ...
    (t >= firingTime & t < firingTime + triseI)    .* IriseT1(t) + ...
    (t >= firingTime + triseI & t < maxFiringTime) .* (Vs(t) - (2 * Vforward) - E)/R + ...
    (t >= maxFiringTime & t < maxFiringTime + tfallI) .* IfallT1(t) + ...
    (t >= maxFiringTime + tfallI) .* Ileak ;

    % Constructing full instantaneous I function for T2,4
    IT2 = @(t) ...
    (t < (T/2) + firingTime)              .* Ileak + ...
    (t >= (T/2) + firingTime & t < (T/2) + firingTime + triseI)    .* IriseT2(t) + ...
    (t >= (T/2) + firingTime + triseI & t < (T/2) + maxFiringTime) .* (-Vs(t) - (2 * Vforward) - E)/R + ...
    (t >= (T/2) + maxFiringTime & t < (T/2) + maxFiringTime + tfallI) .* IfallT2(t) + ...
    (t >= (T/2) + maxFiringTime + tfallI) .* Ileak ;

    t = linspace(0,T,samplingRate);
    % instantaneous Thyristor power T1 + T3
    pT1 = abs(2 .*VT1(t) .* IT1(t));

    % instantaneous Thyristor power T2 + T4
    pT2 = abs(2 .*VT2(t) .* IT2(t));
    
    % average summation of thyristor power losses
    pT = (1/T)*trapz(t,pT1+pT2);
   

    % Load Current
    IL = IT1(t) + IT2(t);


    % Power loss due to resistor
    pLossRes = (calcRMS(t,IL,1/T)).^2 *R;

    powerloss = pT + pLossRes;



end





