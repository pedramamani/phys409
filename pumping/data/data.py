from toolbox import myextract
from data_config import *


class DataSet:
    def __init__(self, config, number):
        dir_set = ASSETS_DIR / config.folder / SET_NAME.format(number=number)
        self._ch1_filepath = str(dir_set / CSV_NAME.format(number=number, channel=1))
        self._ch2_filepath = str(dir_set / CSV_NAME.format(number=number, channel=2))
        self._time = self._ch1 = self._ch2 = None

    def time(self):
        if self._time is None:
            self._extract()
        return self._time

    def ch1(self):
        if self._ch1 is None:
            self._extract()
        return self._ch1

    def ch2(self):
        if self._ch2 is None:
            self._extract()
        return self._ch2

    def _extract(self):
        self._time, self._ch1 = myextract.extract(self._ch1_filepath, xcol=TIME_COLUMN, ycol=VOLTAGE_COLUMN)
        _, self._ch2 = myextract.extract(self._ch2_filepath, xcol=TIME_COLUMN, ycol=VOLTAGE_COLUMN)


class LinearSet(DataSet):
    def __init__(self, number):
        DataSet.__init__(self, LINEAR, number)
        self.rf_frequency, self.main_coil_current = LINEAR.settings[number]


class QuadraticSet(DataSet):
    def __init__(self, number):
        DataSet.__init__(self, QUADRATIC, number)
        self.rf_frequency, self.main_coil_current, self.rf_amplitude = QUADRATIC.settings[number]


class TransLightSet(DataSet):
    def __init__(self, number):
        DataSet.__init__(self, TRANS_LIGHT, number)
        self.polarizer_angle = TRANS_LIGHT.settings[number]


class TransRfSet(DataSet):
    def __init__(self, number):
        DataSet.__init__(self, TRANS_RF, number)
        self.tuning, self.rf_amplitude = TRANS_RF.settings[number]


class Linear:
    temperature = TEMPERATURE
    sample_interval = LINEAR.sample_interval
    ch1_label, ch2_label = LINEAR.labels
    ch1_bin, ch2_bin = [v / POINTS_PER_SCALE for v in LINEAR.voltage_scales]
    sets = {n: LinearSet(n) for n in LINEAR.settings}

    rf_amplitude = LINEAR.rf_amplitude


class Quadratic:
    sample_interval = QUADRATIC.sample_interval
    ch1_label, ch2_label = QUADRATIC.labels
    ch1_bin, ch2_bin = [v / POINTS_PER_SCALE for v in QUADRATIC.voltage_scales]
    sets = {n: QuadraticSet(n) for n in QUADRATIC.settings}


class TransLight:
    sample_interval = TRANS_LIGHT.sample_interval
    ch1_label, ch2_label = TRANS_LIGHT.labels
    ch1_bin, ch2_bin = [v / POINTS_PER_SCALE for v in TRANS_LIGHT.voltage_scales]
    sets = {n: TransLightSet(n) for n in TRANS_LIGHT.settings}

    rf_frequency = TRANS_LIGHT.rf_frequency
    rf_amplitude = TRANS_LIGHT.rf_amplitude
    modulation_frequency = TRANS_LIGHT.modulation_frequency
    modulation_amplitude = TRANS_LIGHT.modulation_amplitude
    tuning = TRANS_LIGHT.tuning


class TransRf:
    sample_interval = TRANS_RF.sample_interval
    ch1_label, ch2_label = TRANS_RF.labels
    ch1_bin, ch2_bin = [v / POINTS_PER_SCALE for v in TRANS_RF.voltage_scales]
    sets = {n: TransRfSet(n) for n in TRANS_RF.settings}

    rf_frequency = TRANS_RF.rf_frequency
    modulation_frequency = TRANS_RF.modulation_frequency
    modulation_amplitude = TRANS_RF.modulation_amplitude


class CrossSection:
    ch1_label, ch2_label, ch3_label = CROSS_SECTION.labels
    temperature_celsius = CROSS_SECTION.temperature_celsius
    signal_inc = CROSS_SECTION.signal_inc
    signal_dec = CROSS_SECTION.signal_dec
    densities = CROSS_SECTION.densities


class Data:
    linear = Linear
    quadratic = Quadratic
    trans_light = TransLight
    trans_rf = TransRf
    cross_section = CrossSection

    temperature = TEMPERATURE
    count_points = COUNT_POINTS
