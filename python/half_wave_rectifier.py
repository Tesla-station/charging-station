import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def _norm_angle(theta: float) -> float:
    """Normalize angle to [0, 2*pi).
    Works for scalars and numpy arrays (uses modulo).
    """
    return np.mod(theta, 2 * np.pi)


def find_range_of_alphas_and_beta(E: float, Vm: float):
    """
    Return two angles (alpha1, alpha2) in radians in [0, 2*pi)
    that satisfy sin(theta) = E/Vm. Raises ValueError if |E/Vm|>1.
    The two solutions are arcsin(frac) and pi - arcsin(frac),
    normalized into [0, 2*pi).
    """
    frac = E / Vm
    if abs(frac) > 1.0:
        raise ValueError(f"|E/Vm| must be <= 1. Got {frac}")
    a = np.arcsin(frac)
    t1 = _norm_angle(a)
    t2 = _norm_angle(np.pi - a)
    return t1, t2


def _make_alpha_array(a1: float, a2: float, step: float = 0.01) -> np.ndarray:
    """Create a numpy array of angles from a1 to a2 (exclusive) with step.
    Handles wrap-around if a1 > a2 by concatenating [a1, 2pi) and [0, a2).
    """
    a1 = float(a1)
    a2 = float(a2)
    step = float(step)
    if step <= 0:
        raise ValueError("step must be positive")
    two_pi = 2 * np.pi
    if np.isclose(a1, a2):
        return np.array([a1])
    if a1 < a2:
        return np.arange(a1, a2, step, dtype=float)
    # wrap-around
    part1 = np.arange(a1, two_pi, step, dtype=float)
    part2 = np.arange(0.0, a2, step, dtype=float)
    if part1.size and part2.size:
        return np.concatenate((part1, part2))
    return part1 if part1.size else part2

def vs_t(x, Vm):
    return Vm*np.sin(x)

def i_t(x, alpha, beta, Vm, E, esr):
    if alpha <= x <= beta:
        return (Vm*np.sin(x)-E)/esr
    else:
        return 0.0  # for 2pi-(beta-alpha)
    
def i_t_including_leakage_simplified(x, alpha, beta, Vm, E, esr, ileak):
    if alpha <= x <= beta:
        return (Vm*np.sin(x)-E)/esr
    else:
        return ileak  # for 2pi-(beta-alpha)

def i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, tr, tf, freq):
    # Convert rise/fall times to angular ranges
    omega = 2*np.pi*freq
    d_theta_r = omega * tr   # rise-time converted to angle
    d_theta_f = omega * tf   # fall-time converted to angle

    natural_I = i_t(x, alpha, beta, Vm, E, esr)

    # Region 1: before α (Ileak)
    if x < alpha:
        return ileak

    # Region 2: α → α + d_theta_r (rising current to natural I)
    if alpha <= x < alpha + d_theta_r:
        m, b = line_from_points(alpha, ileak, alpha + d_theta_r, i_t(alpha + d_theta_r, alpha, beta, Vm, E, esr))
        return m*x + b  # RISE: ileak → natural

    # Region 3: The natural current region
    if alpha + d_theta_r <= x < beta:
        return natural_I

    # Region 4: β → β + theta_f (falling current back to ileak)
    if beta <= x < beta + d_theta_f:
        m, b = line_from_points(beta, i_t(beta, alpha, beta, Vm, E, esr), beta + d_theta_f, -ileak)
        return m*x + b # FALL: natural → -ileak

    # Region 5: after β + theta_f (-Ileak)
    return -ileak


def Vak_including_drop_simplified(x, alpha, beta, Vm, E, v_drop):
    if alpha <= x <= beta:
        return v_drop
    else:
        return Vm*np.sin(x)-E  # for 2pi-(beta-alpha)

def plot_Vak_with_t_simplified(x, alpha, beta, Vm, E, v_drop):
    y = [Vak_including_drop_simplified(xi, alpha, beta, Vm, E, v_drop) for xi in x]
    plt.plot(x, y)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Vak (V)")
    plt.title("Voltage across the load Vak vs Angle")
    plt.grid(True)
    plt.show()

def line_from_points(x1, y1, x2, y2):
    # Slope (m)
    m = (y2 - y1) / (x2 - x1)
    # Intercept (b)
    b = y1 - m * x1
    return m, b

def Vak_including_drop(x, alpha, beta, Vm, E, v_drop, tr, tf, freq):
    # Convert rise/fall times to angular ranges
    omega = 2*np.pi*freq
    d_theta_r = omega * tr   # rise-time converted to angle 
    d_theta_f = omega * tf   # fall-time converted to angle 

    natural_v = vs_t(x, Vm) - E

    # Region 1: before α (natural waveform)
    if x < alpha:
        return natural_v

    # Region 2: α → α + d_theta_f (falling voltage to v_drop)
    if alpha <= x < alpha + d_theta_f:
        m, b = line_from_points(alpha, vs_t(alpha, Vm)-E, alpha + d_theta_f, v_drop)
        return m*x + b  # FALL: natural → v_drop

    # Region 3: constant v_drop
    if alpha + d_theta_f <= x < beta:
        return v_drop

    # Region 4: β → β + theta_r (rising voltage back to natural)
    if beta <= x < beta + d_theta_r:
        m, b = line_from_points(beta, v_drop, beta + d_theta_r, vs_t(beta + d_theta_r, Vm)-E)
        return m*x + b # RISE: v_drop → natural

    # Region 5: after β + theta_r (natural waveform)
    return natural_v


def plot_Vak_with_t(x, alpha, beta, Vm, E, v_drop, tr, tf, freq):
    y = [Vak_including_drop(xi, alpha, beta, Vm, E, v_drop, tr, tf, freq) for xi in x]
    plt.plot(x, y)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Vak (V)")
    plt.title("Voltage across the load Vak vs Angle")
    plt.grid(True)  
    plt.show()

def plot_i_with_t(x, alpha, beta, Vm, E, esr, ileak, tr, tf, freq):
    y = [i_t_including_leakage(xi, alpha, beta, Vm, E, esr, ileak, tr, tf, freq) for xi in x]
    plt.plot(x, y)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Current i (A)")
    plt.title("Current i vs Angle")
    plt.grid(True)
    plt.show()

def compute_Iavg(Vm, E, esr, alpha, beta) -> float:
    i = lambda x: i_t(x, alpha, beta, Vm, E, esr)
    result, _ = quad(i, 0, 2*np.pi, limit=1000)
    Iavg = (1.0 / (2.0 * np.pi)) * result
    return Iavg

def compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, tr, tf, freq):
    integral, _ = quad(lambda x: i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, tr, tf, freq)**2, 0, 2.0 * np.pi, limit=1000)
    return np.sqrt(integral / (2.0 * np.pi))

def compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq, Irms):
    Vak = Vak_including_drop(x, alpha, beta, Vm, E, v_drop, tr, tf, freq)
    i = i_t_including_leakage(x, alpha, beta, Vm, E, esr, ileak, tr, tf, freq)
    return Vak * i + esr * i**2

def compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq, Irms):
    power_loss_func = lambda x: compute_power_loss(x, alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq, Irms)
    result, _ = quad(power_loss_func, 0, 2*np.pi, limit=1000)
    P_loss_avg = (1.0 / (2.0 * np.pi)) * result
    return P_loss_avg

def plot_power_loss_with_alpha(alpha_array, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq):
    y = []
    for alpha in alpha_array:
        Irms = compute_i_t_including_leakage_rms(alpha, beta, Vm, E, esr, ileak, tr, tf, freq)
        P_loss_avg = compute_average_power_loss(alpha, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq, Irms)
        y.append(P_loss_avg)
    plt.plot(alpha_array, y)
    plt.xlabel("Firing Angle (rad)")
    plt.ylabel("Average Power Loss (W)")
    plt.title("Average Power Loss vs Firing Angle")
    plt.grid(True)
    plt.show()

def plot_power_loss_with_t(x, Vm, E, esr, ileak, tr, tf, freq, alpha, beta):
    y = [compute_power_loss(xi, alpha, beta, Vm, E, esr, ileak, tr, tf, freq) for xi in x]
    plt.plot(x, y)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Power Loss (W)")
    plt.title("Power Loss vs Angle")
    plt.grid(True)
    plt.show()


def time_of_charging(Vrms, E, esr, C, step: float = 0.01, degrees: bool = True, show: bool = True,
                     charging_time_scale: str = 'linear', initial_soc: float = 0.0, required_charging_time: float = 0.0):
    Vm = np.sqrt(2.0) * Vrms
    alpha1, beta = find_range_of_alphas_and_beta(E, Vm)
    alphas = _make_alpha_array(alpha1, beta, step=step)
    
    # here we compute Iavg for each alpha
    Iavg = []
    for alpha in alphas:
        Iavg.append(compute_Iavg(Vm, E, esr, alpha, beta))
    Iavg = np.array(Iavg)

    # Compute charging time with respect to alphas (C / Iavg) when capacity C provided.
    # Only valid when Iavg > 0. Mask non-positive currents to NaN so they are not plotted.
    charging_time = C/Iavg
    if C is not None:
        # Broadcast: C scalar divided by Iavg array
        charging_time = np.full_like(Iavg, np.nan, dtype=float)
        positive = Iavg > 0
        charging_time[positive] = C / Iavg[positive]

    # Plot
    x = np.degrees(alphas) if degrees else alphas
    xlabel = "Firing angle (deg)" if degrees else "Firing angle (rad)"
    # Larger figure size so small numbers are visible
    fig, ax1 = plt.subplots(figsize=(14, 8))
    ax1.plot(x, Iavg, '-o', markersize=4, label='Iavg (A)', color='tab:blue')
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel('Iavg (A)', color='tab:blue', fontsize=12)
    ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=11)
    ax1.tick_params(axis='x', labelsize=11)
    ax1.grid(True, which='both', axis='both', linestyle='--', alpha=0.4)
    if charging_time is not None:
        ax2 = ax1.twinx()
        ax2.plot(x, charging_time, '-s', markersize=4, label='Charging time (s)', color='tab:orange')
        ylabel = 'Charging time (s)'
        if charging_time_scale and charging_time_scale != 'linear':
            ylabel += f' ({charging_time_scale})'
            ax2.set_yscale(charging_time_scale)
        ax2.set_ylabel(ylabel, color='tab:orange', fontsize=12)
        ax2.tick_params(axis='y', labelcolor='tab:orange', labelsize=11)
        # Auto-scale y-axis to fit all data with small padding
        ax2.margins(y=0.1)
    plt.title(f'Iavg and charging time vs firing angle (Vrms={Vrms} V, E={E} V, esr={esr} ohms, C={C} Ah)', fontsize=13, fontweight='bold')
    fig.tight_layout()

    if required_charging_time > 0.0 and initial_soc > 0.0:
        new_soc = (Iavg * required_charging_time  + (initial_soc/100.0)*C) / C
        new_soc = new_soc * 100.0  # convert to percentage
        new_soc = np.clip(new_soc, 0.0, 100.0)   # clip to 100%

        # Create second plot for final state of charge
        fig2, ax3 = plt.subplots(figsize=(14, 6))
        ax3.plot(x, new_soc, '-o', markersize=4, color='tab:green')
        ax3.set_xlabel(xlabel, fontsize=12)
        ax3.set_ylabel("Final State of Charge (%)", fontsize=12, color='tab:green')
        ax3.tick_params(axis='y', labelcolor='tab:green', labelsize=11)
        ax3.tick_params(axis='x', labelsize=11)
        ax3.grid(True, linestyle='--', alpha=0.4)
        plt.title(f"Final SOC vs Firing Angle\nInit SOC={initial_soc}%, Time={required_charging_time}h, C={C}Ah",fontsize=13,fontweight='bold')

        fig2.tight_layout()

    if show:
        plt.show()

    return alphas, Iavg, charging_time    

if __name__ == "__main__":
    Vrms = 60
    E = 12
    esr = 4.256
    capacity = 100
    v_drop = 1.2
    tr = 10*10**(-6)  # rise time in seconds
    tf = 5*10**(-6)  # fall time in seconds
    freq = 50  # frequency in Hz
    ileak = 0.02  # leakage current in Amperes
    Vm = np.sqrt(2)*Vrms
    alpha1, beta = find_range_of_alphas_and_beta(E, Vm)
    alphas = _make_alpha_array(alpha1, beta, step=0.01)
    x = np.linspace(0, 2*np.pi, 10000)
    # plot_Vak_with_t_simplified(x, alpha1, beta, Vm, E, v_drop)
    time_of_charging(Vrms, E, esr, capacity, step=0.01, charging_time_scale='log', initial_soc=30.0, required_charging_time=0.5)
    # plot_Vak_with_t(x, alpha1, beta, Vm, E, v_drop, tr, tf, freq)
    # plot_i_with_t(x, alpha1, beta, Vm, E, esr, ileak=ileak, tr=tr, tf=tf, freq=freq)
    # plot_power_loss_with_t(x, Vm, E, esr, ileak, tr, tf, freq, alpha1, beta)
    # plot_power_loss_with_alpha(alphas, beta, Vm, E, esr, ileak, v_drop, tr, tf, freq)
    # print(compute_i_t_including_leakage_rms(alpha1, beta, Vm, E, esr, ileak, tr, tf, freq))