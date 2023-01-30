from toolbox import myplot, myprocess
from rb_atom import RbAtom
from data import Data
from dynamics import dynamics
from pumping_config import *
from toolbox.constants import h, PRE, muB
import numpy as np
import scipy.optimize as optimize


def plot_energies(Bs, Rbx):
    atom = RbAtom(Rbx)
    list_energies = np.array([atom.energy_levels(B) for B in Bs]).T
    p = myplot.Plot()
    for list_energy in list_energies:
        p.line(Bs / GAUSS, list_energy / (h * PRE.G), format_='k')
    print(f'Plotting {RB_NAMES[Rbx]} Energy Levels')
    p.show(xlabel='B (Gauss)', ylabel='E / h (GHz)')


def plot_rf_freqs(Bs, Rbx):
    atom = RbAtom(Rbx)
    list_frequencies = np.array([atom.rf_frequencies(B) for B in Bs]).T
    p = myplot.Plot()
    for list_frequency in list_frequencies:
        p.line(Bs / GAUSS, list_frequency / PRE.M, format_='k')
    print(f'Plotting {RB_NAMES[Rbx]} RF Transition Frequencies')
    p.show(xlabel='B (Gauss)', ylabel='Frequency (GHz)')


def rf_frequency_to_B(target_frequency, Rbx):
    atom = RbAtom(Rbx)

    def B_to_difference(B):
        return np.average(atom.rf_frequencies(B)) - target_frequency

    return optimize.newton(B_to_difference, 1E-5)


def plot_dynamics(I):
    ts, Nks, Nis, Abs = dynamics(I)
    ts /= PRE.u
    p = myplot.Plot()
    for Nk in Nks:
        p.line(ts, Nk, format_='k')
    for Ni in Nis:
        p.line(ts, Ni, format_='r')
    print(f'Plotting Rb87 State Populations at I={I} W/m2')
    p.show(xlabel='Time ($\mu$s)', ylabel='Quantum State Populations', yrange=(0, 1), grid=True,
           legend=['Ground States'] + [None] * 7 + ['Excited States'] + [None] * 7)

    p = myplot.Plot()
    p.line(ts, Abs, format_='k')
    print(f'Plotting Rb87 Absorption at I={I} W/m2')
    p.show(xlabel='Time ($\mu$s)', ylabel='Absorption (a.u.)', yrange=(0, 1), grid=True)


def process_linear(number, show=False):
    set_ = Data.linear.sets[number]
    time, ch1, ch2 = set_.time(), set_.ch1(), set_.ch2()
    edge1, edge2, *_ = myprocess.index_edges(ch1)
    time, ch1, ch2 = time[edge1: edge2], ch1[edge1: edge2], ch2[edge1: edge2]

    # Linear fit to current ramp
    ch1_errors = [Data.linear.ch1_bin / 2] * len(time)
    current_func = myprocess.fit_func(myprocess.Function.line, time, ch1, ch1_errors)
    current, current_errors = current_func(time, Data.linear.sample_interval / 2)
    if show:
        p = myplot.Plot()
        p.line(time, ch1, format_='.r', marker_size=1)
        p.line(time, current, format_='k')
        p.show(xlabel='Time (s)', ylabel=Data.linear.ch1_label, legend=['data_assets Points', 'Linear Fit'])

    ch2 = myprocess.denoise(ch2, 0.1)
    peak_indices = myprocess.index_peaks(-ch2, distance=10, prominence=0.1)
    print(
        f'Transparency at Linear RF Frequency {set_.rf_frequency / PRE.k}kHz, Main Current {set_.main_coil_current}A\n'
        f'Peak Currents (A): {current[peak_indices]} +- {current_errors[peak_indices]}\n'
        f'Peak Transparencies (a.u.): {ch2[peak_indices]} +- {[Data.linear.ch2_bin / 2] * len(peak_indices)}\n')
    if show:
        p = myplot.Plot()
        p.line(current, ch2, format_='k')
        for i, y in zip(peak_indices, ch2[peak_indices]):
            p.line(current[i], y, format_='.r')
        p.show(xlabel=Data.linear.ch1_label, ylabel=Data.linear.ch2_label)


def process_quadratic(number, show=False):
    set_ = Data.quadratic.sets[number]
    time, ch1, ch2 = set_.time(), set_.ch1(), set_.ch2()
    edge1, edge2, *_ = myprocess.index_edges(ch1)
    time, ch1, ch2 = time[edge1: edge2], ch1[edge1: edge2], ch2[edge1: edge2]

    # Linear fit to current ramp
    ch1_errors = [Data.quadratic.ch1_bin / 2] * len(time)
    current_func = myprocess.fit_func(myprocess.Function.line, time, ch1, ch1_errors)
    current, current_errors = current_func(time, Data.quadratic.sample_interval / 2)
    if show:
        p = myplot.Plot()
        p.line(time, ch1, format_='.r', marker_size=1)
        p.line(time, current, format_='k')
        p.show(xlabel='Time (s)', ylabel=Data.quadratic.ch1_label, legend=['data_assets Points', 'Linear Fit'])

    offset, padding = 0.2, 0.05
    ch2 = myprocess.denoise(ch2, 0.1) + offset
    peak_indices = myprocess.index_peaks(-ch2, distance=10, prominence=0.05)
    print(f'Transparency at Quadratic RF Frequency {set_.rf_frequency / PRE.M}MHz, '
          f'Main Current {set_.main_coil_current}A, RF Amplitude {set_.rf_amplitude}V\n'
          f'Peak Currents (A): {current[peak_indices]} +- {current_errors[peak_indices]}\n'
          f'Peak Transparencies (a.u.): {ch2[peak_indices]} +- {[Data.quadratic.ch2_bin / 2] * len(peak_indices)}\n')
    if show:
        p = myplot.Plot()
        p.line(current, ch2, format_='k')
        for i, y in zip(peak_indices, ch2[peak_indices]):
            p.line(current[i], y, format_='.r')
        p.show(xlabel=Data.quadratic.ch1_label, ylabel=Data.quadratic.ch2_label)


def process_trans_rf(number, show=False):
    p = myplot.Plot()
    set_ = Data.trans_rf.sets[number]
    time, ch1, ch2 = set_.time(), set_.ch1(), set_.ch2()
    edge, *_ = myprocess.index_edges(ch1)
    time, ch1, ch2 = time[edge:], ch1[edge:], ch2[edge:]

    # Rabi fit to detector signal
    params, errors = myprocess.fit(myprocess.Function.rabi, time, ch2, [Data.trans_rf.ch2_bin / 2] * len(time))
    print(f'{RB_NAMES[set_.tuning]} Rabi Oscillations at RF Amplitude {set_.rf_amplitude}V\n'
          f'Decay Constant (1/s): {params[1]} +- {errors[1]}\n'
          f'Sinusoidal Decay Constant (1/s): {params[3]} +- {errors[3]}\n'
          f'Angular Frequency (1/s): {np.abs(params[4])} +- {np.abs(errors[4])}\n')
    if show:
        p.line(time / PRE.m, ch2, format_='.r', marker_size=2)
        p.line(time / PRE.m, myprocess.Function.rabi(time, *params), format_='k')
        p.show(xlabel='Time (ms)', ylabel=Data.linear.ch2_label, legend=['data_assets Points', 'Fit'])


def process_trans_light(number, show=False):
    set_ = Data.trans_light.sets[number]
    time, ch1, ch2 = set_.time(), set_.ch1(), set_.ch2()
    edges = myprocess.index_edges(ch1, prominence=0.1)
    edge1, edge2 = edges[0], edges[-1]
    time, ch1, ch2 = time[edge1: edge2], ch1[edge1: edge2], ch2[edge1: edge2]

    # Decay fit to signal
    params, errors = myprocess.fit(myprocess.Function.exp_decay, time, ch2, [Data.trans_rf.ch2_bin / 2] * len(time))
    print(f'{RB_NAMES[Data.trans_light.tuning]} Transient at Polarizer Angle {set_.polarizer_angle}°\n'
          f'Decay Constant (1/s): {params[1]} +- {errors[1]}\n')
    if show:
        p = myplot.Plot()
        p.line(time, ch2, format_='.r', marker_size=1)
        p.line(time, myprocess.Function.exp_decay(time, *params))
        p.show(xlabel='Time (s)', ylabel=Data.linear.ch2_label, legend=['data_assets Points', 'Decay Fit'])


class ProcessedData:
    class linear:  # all peaks
        fRF = [60E3, 80E3, 100E3, 120E3, 140E3, 160E3, 180E3, 300E3, 600E3, 1200E3, 800E3]
        Imain = [0, 0, 0, 0, 0, 0, 0, 52E-3, 105E-3, 200E-3, 200E-3]
        sigma_Isweep = 0.0004
        Isweep0 = [0.2924, 0.2943, 0.2942, 0.2939, 0.2934, 0.2941, 0.2935, None, None, None, None]
        Isweep87 = [0.4307, 0.4798, 0.5272, 0.5745, 0.6217, 0.6686, 0.7159, 0.2466, 0.1821, 0.2155, None]
        Isweep85 = [0.5028, 0.5740, 0.6444, 0.7148, 0.7858, 0.8562, 0.9268, 0.5987, 0.8836, None, 0.2151]

    class quadratic:  # only the 1V data_assets
        fRF85 = [6.35E6, 4.2E6, 8.26E6]
        Imain85 = [1.505, 1.008, 1.998]
        sigma_Isweeps85 = [2E-4, 8E-5, 2E-4]
        Isweeps85 = [[0.3391, 0.4312, 0.5229, 0.6146, 0.7072, 0.7998],
                     [0.41344, 0.45421, 0.49514, 0.53574, 0.57618, 0.61728],
                     [0.1230, 0.2804, 0.4383, 0.5947, 0.7545, 0.9137]]

        fRF87 = [9.3E6, 6.2E6, 12.4E6]
        Imain87 = [1.505, 1.008, 1.998]
        sigma_Isweeps87 = [2E-4, 7E-5, 2E-4]
        Isweeps87 = [[0.3206, 0.3812, 0.4414, 0.5002],
                     [0.25690, 0.28228, 0.30869, 0.33540],
                     [0.3840, 0.4891, 0.5959, 0.7020]]

    class trans_rf:  # all fits except 5V for Rb85
        VRF87 = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
        decay87 = [340.2, 482, 359.3, 472, 473, 473, 586, 576, 655]
        sigma_decay87 = [0.9, 2, 0.7, 2, 1, 1, 2, 3, 4]
        omega87 = [4055, 6028, 8106.4, 10396, 12007, 14205, 16015, 18278, 20423]
        sigma_omega87 = [1, 2, 0.7, 2, 1, 2, 3, 3, 4]

        VRF85 = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
        decay85 = [501.9, 372.4, 393.5, 384.2, 374.0, 455.5, 448.9, 449.2]
        sigma_decay85 = [0.8, 0.2, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3]
        omega85 = [3198, 4275.9, 5791.5, 7224.9, 8693.4, 9289.7, 10907.2, 12498.2]
        sigma_omega85 = [1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4]


def fit_BvsI(show=False):
    L = ProcessedData.linear
    Isweep = np.array([[L.Isweep0[i], L.Isweep87[i], L.Isweep85[i]]
                       for i, I in enumerate(L.Imain) if I == 0]).flatten()
    Bfield = np.array([[0, rf_frequency_to_B(L.fRF[i], Rb87), rf_frequency_to_B(L.fRF[i], Rb85)]
                       for i, I in enumerate(L.Imain) if I == 0]).flatten()
    Bsweep = myprocess.fit_func(myprocess.Function.line, Isweep, Bfield, [1E-5] * len(Bfield), absolute_sigma=False)
    if show:
        p = myplot.Plot()
        p.line(Isweep, Bsweep(Isweep, sigma=L.sigma_Isweep)[0] / PRE.u, format_='k')
        p.line(Isweep, Bfield / PRE.u, format_='.r')
        p.show(xlabel='Sweep Current (a.u.)', ylabel='Sweep Magnetic Field ($\mu$T)',
               legend=['data_assets Points', 'Linear Fit'])

    Imain = [52E-3, 105E-3, 200E-3, 52E-3, 105E-3, 200E-3]
    Bfield = np.subtract([rf_frequency_to_B(f, Rb87) for f in [300E3, 600E3, 1200E3]] +
                         [rf_frequency_to_B(f, Rb85) for f in [300E3, 600E3, 800E3]],
                         Bsweep([0.2466, 0.1821, 0.2155, 0.5987, 0.8836, 0.2151])[0])
    Bmain = myprocess.fit_func(myprocess.Function.line, Imain, Bfield, [1E-5] * len(Bfield), absolute_sigma=False)
    if show:
        p = myplot.Plot()
        p.line(Imain, Bmain(Imain)[0] / PRE.u, format_='k')
        p.line(Imain, Bfield / PRE.u, format_='.r')
        p.show(xlabel='Main Current (a.u.)', ylabel='Main Magnetic Field ($\mu$T)',
               legend=['data_assets Points', 'Linear Fit'])

    def B(Is, Im, sigma_Is=None, sigma_Im=None):
        Bs = Bsweep(Is, sigma=sigma_Is)
        Bm = Bmain(Im, sigma=sigma_Im)
        return Bs[0] + Bm[0], Bs[1] + Bm[1]

    return B


def fit_fRFvsB(show=False):
    L = ProcessedData.linear
    B = fit_BvsI()
    f87s = np.array([L.fRF[i] for i, I in enumerate(L.Isweep87) if I is not None])
    f85s = np.array([L.fRF[i] for i, I in enumerate(L.Isweep85) if I is not None])
    B87s = np.array([B(I, L.Imain[i], sigma_Is=L.sigma_Isweep)[0] for i, I in enumerate(L.Isweep87) if I is not None])
    B85s = np.array([B(I, L.Imain[i], sigma_Is=L.sigma_Isweep)[0] for i, I in enumerate(L.Isweep85) if I is not None])
    params87, errors87 = myprocess.fit(myprocess.Function.line, B87s, f87s, [1E-3] * len(f87s), absolute_sigma=False)
    params85, errors85 = myprocess.fit(myprocess.Function.line, B85s, f85s, [1E-3] * len(f85s), absolute_sigma=False)
    print(f'g87 = {params87[0] * h / muB} +- {errors87[0] * h / muB}\n'
          f'g85 = {params85[0] * h / muB} +- {errors85[0] * h / muB}\n'
          f'offset87 (Hz) = {params87[1]} +- {errors87[1]}\n'
          f'offset85 (Hz) = {params85[1]} +- {errors85[1]}')

    if show:
        p = myplot.Plot()
        p.line(B87s / PRE.u, f87s / PRE.M, format_='ob', marker_size=5)
        p.line(B87s / PRE.u, myprocess.Function.line(B87s, *params87) / PRE.M, format_='--b')
        p.line(B85s / PRE.u, f85s / PRE.M, format_='*r', marker_size=6)
        p.line(B85s / PRE.u, myprocess.Function.line(B85s, *params85) / PRE.M, format_='r')
        p.show(xlabel='Magnetic Field Along Main Axis ($\mu$T)', ylabel='Resonant Radio Frequency (MHz)',
               legend=['$^{87}$Rb data_assets Points', '$^{87}$Rb Linear Fit', '$^{85}$Rb data_assets Points',
                       '$^{85}$Rb Linear Fit'],
               xrange=[None, None])


def fit_fRabivsVRF(show=False):
    T = ProcessedData.trans_rf
    params87, errors87 = myprocess.fit(myprocess.Function.line, T.VRF87, T.omega87, T.sigma_omega87)
    params85, errors85 = myprocess.fit(myprocess.Function.line, T.VRF85, T.omega85, T.sigma_omega85)
    print(f'g87 = {params87[0] * h / muB} +- {errors87[0] * h / muB}\n'
          f'g85 = {params85[0] * h / muB} +- {errors85[0] * h / muB}\n'
          f'offset87 (Hz) = {params87[1]} +- {errors87[1]}\n'
          f'offset85 (Hz) = {params85[1]} +- {errors85[1]}')

    if show:
        p = myplot.Plot()
        p.line(T.VRF87, np.array(T.omega87) / PRE.k, format_='ob', marker_size=5)
        p.line(T.VRF87, myprocess.Function.line(T.VRF87, *params87) / PRE.k, format_='--b')
        p.line(T.VRF85, np.array(T.omega85) / PRE.k, format_='*r', marker_size=6)
        p.line(T.VRF85, myprocess.Function.line(T.VRF85, *params85) / PRE.k, format_='r')
        p.show(xlabel='Radio Frequency Amplitude (V)', ylabel='Rabi Frequency (kHz)',
               legend=['$^{87}$Rb data_assets Points', '$^{87}$Rb Linear Fit', '$^{85}$Rb data_assets Points',
                       '$^{85}$Rb Linear Fit'],
               xrange=[None, None])


if __name__ == '__main__':
    plot_energies(np.linspace(0, 0.1, 40), Rb85)
    plot_rf_freqs(np.linspace(0, 0.1, 40), Rb85)
    plot_dynamics(1E1)

    for n in Data.linear.sets:
        process_linear(n, show=True)
    for n in Data.quadratic.sets:
        process_quadratic(n, show=True)
    for n in Data.trans_rf.sets:
        process_trans_rf(n, show=True)
    for n in Data.trans_light.sets:
        process_trans_light(n, show=True)

    fit_fRFvsB(show=True)
    fit_fRabivsVRF(show=True)
