import os
from pathlib import Path
from pumping_config import Rb85, Rb87

DIR = Path(os.path.dirname(os.path.abspath(__file__)))
NAME = os.path.basename(DIR)
ASSETS_DIR = DIR / f'{NAME}_assets'

SET_NAME = 'ALL00{number:02d}'
CSV_NAME = 'F00{number:02d}CH{channel}.CSV'
TIME_COLUMN, VOLTAGE_COLUMN, N_COLUMNS = 3, 4, 6
POINTS_PER_SCALE = 25  # number of data_assets points per voltage scale

TEMPERATURE = 50 + 273
COUNT_POINTS = 2500


class LINEAR:
    folder = 'linear'
    labels = ['Sweep Current (a.u.)', 'Detector Signal (a.u.)']
    voltage_scales = [5E-1, 1E-1]
    sample_interval = 4E-3

    rf_amplitude = 1
    settings = {  # RF frequency, Main coil current
        37: (60E3, 0),
        38: (80E3, 0),
        39: (100E3, 0),
        40: (120E3, 0),
        41: (140E3, 0),
        42: (160E3, 0),
        43: (180E3, 0),
        44: (300E3, 52E-3),
        45: (600E3, 105E-3),
        47: (1200E3, 200E-3),
        48: (800E3, 200E-3)
    }


class QUADRATIC:
    folder = 'quadratic'
    labels = ['Sweep Current (a.u.)', 'Detector Signal (a.u.)']
    voltage_scales = [5E-1, 1]
    sample_interval = 1E-3

    settings = {  # RF frequency, Main coil current, RF amplitude
        52: (6.25E6, 1.505, 3),
        53: (6.25E6, 1.505, 1),
        54: (9.3E6, 1.505, 1),
        55: (9.3E6, 1.505, 3),
        58: (6.2E6, 1.008, 1),
        59: (6.2E6, 1.008, 3),
        60: (4.2E6, 1.008, 3),
        61: (4.2E6, 1.008, 1),
        62: (8.26E6, 1.998, 3),
        63: (8.26E6, 1.998, 1),
        64: (8.26E6, 1.998, 5),
        67: (12.4E6, 1.998, 3),
        68: (12.4E6, 1.998, 5),
        69: (12.4E6, 1.998, 1)
    }


class TRANS_LIGHT:
    folder = 'trans_light'
    labels = ['Square Wave Voltage (V)', 'Detector Signal (a.u.)']
    voltage_scales = [2, 5E-2]
    sample_interval = 1E-4

    rf_frequency = 300E3
    rf_amplitude = 5
    modulation_frequency = 5
    modulation_amplitude = 5
    tuning = Rb85
    settings = {0: 90, 1: 70, 2: 45, 3: 35, 4: 25, 5: 80, 6: 75, 7: 65, 8: 60, 9: 55, 10: 50}  # polarizer angle


class TRANS_RF:
    folder = 'trans_rf'
    labels = ['Square Wave Voltage (V)', 'Detector Signal (a.u.)']
    voltage_scales = [5, 2E-1]
    sample_interval = 4E-6

    rf_frequency = 300E3
    modulation_frequency = 5
    modulation_amplitude = 1
    settings = {  # Tuning, RF amplitude
        0: (Rb87, 2),
        1: (Rb87, 2.5),
        2: (Rb87, 3),
        3: (Rb87, 3.5),
        4: (Rb87, 4),
        5: (Rb87, 4.5),
        6: (Rb87, 5),
        33: (Rb85, 5),
        34: (Rb85, 4.5),
        35: (Rb85, 4),
        36: (Rb85, 3.5),
        37: (Rb85, 3),
        38: (Rb85, 2.5),
        39: (Rb85, 2),
        40: (Rb85, 1.5),
        41: (Rb85, 1),
        42: (Rb87, 1),
        43: (Rb87, 1.5)
    }


class CROSS_SECTION:
    labels = ['Temperature (C)', 'Detector Signal - Increasing Temp (a.u.)',
              'Detector Signal - Decreasing Temp (a.u.)']
    temperature_celsius = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    signal_inc = [1.89, 1.69, 1.63, 1.5, 1.32, 1.15, 0.96, 0.83, 0.71, 0.61, 0.55, 0.51, 0.49, 0.47, 0.46]
    signal_dec = [1.9, 1.71, 1.59, 1.44, 1.27, 1.09, 0.92, 0.78, 0.68, 0.59, 0.54, 0.5, 0.48, 0.47, 0.46]

    densities = [
        [290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400],  # temperature
        [3.3E15, 1.1E16, 2.9E16, 7.5E16, 1.8E17, 4.3E17, 8.3E17, 1.5E18, 3.7E18, 6.3E18, 1.2E19, 2.4E19]  # density
    ]
