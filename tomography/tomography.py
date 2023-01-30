import numpy as np
from constants import pi, PRE
from myplot import Plot

ALPHA_LEN, BETA_LEN = 2 ** 14 + 1, 2 ** 2 + 1  # 12, 1
ALPHA_MAX, BETA_MAX = np.deg2rad(45), np.deg2rad(80)  # 60, 60
ALPHA_STEP, BETA_STEP = 2 * ALPHA_MAX / ALPHA_LEN, 2 * BETA_MAX / BETA_LEN
ANGLE_FACTOR = 1.5


class Tomography:
    def __init__(self):
        self.sources = []
        self.scints = []
        self.area_factors = self.calc_area_factors()

    def add_source(self, position, activity=1):
        source = Source(np.array(position), activity)
        self.sources.append(source)
        return self

    def add_scint(self, position, profile, width, height):
        scint = Scint(np.array(position), profile, width, height)
        self.scints.append(scint)
        return self

    def grids(self, sources_center, sources_angle):
        grids = []
        for source in self.sources:
            x, y, z = source.position
            sin, cos = np.sin(sources_angle), np.cos(sources_angle)
            position = np.array([x * cos - y * sin, x * sin + y * cos, z]) + sources_center

            grid = np.ones((BETA_LEN, ALPHA_LEN))
            for scint in self.scints:
                grid = np.multiply(grid, source.project(position, scint))
            grids.append(grid)
        return grids

    def plot_grids(self, grids):
        for i, grid in enumerate(grids):
            p = Plot()
            extent_degrees = np.rad2deg([-ALPHA_MAX, ALPHA_MAX, -BETA_MAX, BETA_MAX])
            p.cmap(grid, label='Detection Probability (0-1)', extent=extent_degrees)
            p.show(xlabel='$\\alpha$ (°)', ylabel='$\\beta$ (°)',
                   title=f'Detection Probabilities for Source {i + 1}')

    def count(self, sources_center, sources_angle):
        count = 0
        grids = self.grids(sources_center, sources_angle)
        for source, grid in zip(self.sources, grids):
            grid = np.multiply(grid, self.area_factors)
            count += np.sum(grid) * source.activity
        return count

    def sinogram(self, sources_positions, sources_angles):
        counts = np.zeros((len(sources_positions), len(sources_angles)))
        for i, sources_position in enumerate(sources_positions):
            for j, sources_angle in enumerate(sources_angles):
                counts[i, j] = self.count(sources_position, sources_angle)
        return counts

    def calc_area_factors(self):
        ratios = np.zeros((BETA_LEN, ALPHA_LEN))
        for i in range(BETA_LEN):
            for j in range(ALPHA_LEN):
                alpha = (j - ALPHA_LEN // 2) * ALPHA_STEP
                beta = (-i + BETA_LEN // 2) * BETA_STEP
                A00 = pi / 2 - np.arccos(np.sin(alpha) * np.sin(beta))
                A01 = pi / 2 - np.arccos(np.sin(alpha + ALPHA_STEP) * np.sin(beta))
                A10 = pi / 2 - np.arccos(np.sin(alpha) * np.sin(beta + BETA_STEP))
                A11 = pi / 2 - np.arccos(np.sin(alpha + ALPHA_STEP) * np.sin(beta + BETA_STEP))
                ratios[i, j] = (A11 + A00 - A10 - A01) / (2 * pi)
        return ratios


class Source:
    def __init__(self, position, activity):
        self.position = position
        self.activity = activity

    def check(self, position, scint):
        x, y, z = scint.position - position
        if scint.is_left != (x < 0):
            raise SystemError('Source cannot be placed past a scintillator.')
        x, y, z = [-x, y, z] if (x < 0) else [x, -y, -z]
        profile = scint.profile if (x < 0) else np.flip(scint.profile, axis=0)

        alpha_max = max(-np.arctan2(y - scint.w / 2, x), np.arctan2(y + scint.w / 2, x))
        beta_max = max(-np.arctan2(z - scint.h / 2, x), np.arctan2(z + scint.h / 2, x))
        alpha_factor = (scint.dw / x) / (2 * ALPHA_MAX / ALPHA_LEN)
        beta_factor = (scint.dh / x) / (2 * BETA_MAX / BETA_LEN)

        if ALPHA_MAX < alpha_max:
            raise RuntimeError(f'Alpha range is limited. Increase to at least {np.rad2deg(alpha_max) + 1:.0f}°.')
        if BETA_MAX < beta_max:
            raise RuntimeError(f'Beta range is limited. Increase to at least {np.rad2deg(beta_max) + 1:.0f}°.')
        if alpha_factor < ANGLE_FACTOR:
            raise RuntimeError(f'Alpha resolution is limited. Increase grid length to at least '
                               f'{(ALPHA_LEN * ANGLE_FACTOR / alpha_factor) // 2 * 2 + 1:.0f}.')
        if beta_factor < ANGLE_FACTOR:
            raise RuntimeError(f'Beta resolution is limited. Increase grid length to at least '
                               f'{(BETA_LEN * ANGLE_FACTOR / beta_factor) // 2 * 2 + 1:.0f}.')
        return [x, y, z], profile

    def project(self, position, scint):
        [x, y, z], profile = self.check(position, scint)
        projection = np.zeros((BETA_LEN, ALPHA_LEN))
        dim_edges = np.array(profile.shape) + 1

        ih_edges, iw_edges = np.indices(dim_edges)
        w_edges = (iw_edges - scint.w_len / 2) * scint.dw
        h_edges = -(ih_edges - scint.h_len / 2) * scint.dh
        position_edges = np.dstack([np.zeros(dim_edges), w_edges, h_edges]) + [x, y, z]
        position_edges /= np.linalg.norm(position_edges, axis=2)[:, :, None]

        alpha_edges, beta_edges = np.arcsin(position_edges[:, :, 1]), np.arcsin(position_edges[:, :, 2])
        ialpha = np.round(alpha_edges / ALPHA_STEP + ALPHA_LEN / 2).astype(int)
        ibeta = np.round(-beta_edges / BETA_STEP + BETA_LEN / 2).astype(int)

        for i, row in enumerate(profile):
            for j, value in enumerate(row):
                i1_alpha, i2_alpha = min(ialpha[i, j], ialpha[i + 1, j]), max(ialpha[i, j + 1], ialpha[i + 1, j + 1])
                i1_beta, i2_beta = min(ibeta[i, j], ibeta[i, j + 1]), max(ibeta[i + 1, j], ibeta[i + 1, j + 1])
                projection[i1_beta: i2_beta + 1, i1_alpha: i2_alpha + 1] = value
        return projection


class Scint:
    def __init__(self, position, profile, width, height):
        if not (np.all(0 <= profile) and np.all(profile <= 1)):
            raise RuntimeError('Scintillator profile values must be between 0 and 1 (inclusive).')
        if profile.ndim != 2:
            raise RuntimeError('Scintillator profile must be 2-dimensional.')

        h_len, w_len = profile.shape
        x, y, z = position

        if w_len % 2 == 0 or h_len % 2 == 0:
            raise RuntimeError('Scintillator profile must have odd dimensions.')
        if x == 0:
            raise RuntimeError('Scintillator cannot be placed at an x position of 0.')

        self.is_left = x < 0
        self.profile = profile
        self.position = position
        self.w, self.h = width, height
        self.w_len, self.h_len = w_len, h_len
        self.dw, self.dh = width / w_len, height / h_len
