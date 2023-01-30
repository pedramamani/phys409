from pathlib import Path
import os
import cv2
import numpy as np
from skimage.transform import iradon
from scipy.ndimage import rotate

from data import data
from constants import pi, PRE
import myprocess
from tomography import Tomography
from myplot import Plot

DIR = Path(os.path.dirname(os.path.abspath(__file__)))
NAME = os.path.basename(DIR)
ASSETS_DIR = DIR / f'{NAME}_assets'

SCINT_X = 32.8 * PRE.m
D_APERTURE = 2 * 25.41 * PRE.m


def reconstruct(sinogram, filter_=True):
    filter_name = 'ramp' if filter_ else None
    image = iradon(sinogram, filter_name=filter_name)
    image = np.flip(image.T, axis=0)
    return image


def visualize_data(run_number, slice_=True, filter_=True):
    sinogram, positions, angles = data(run_number)

    y_mid = len(positions) // 2
    a_mid = len(angles) // 2
    center_index = myprocess.index_peaks(sinogram[:, 0], distance=y_mid, prominence=0.5)[0]
    sinogram = np.roll(sinogram, y_mid - center_index, axis=0)
    image = reconstruct(sinogram, filter_=filter_)
    image = myprocess.denoise(image, 0.4)
    image = rotate(image, 5, reshape=False)
    positions = positions - positions[y_mid]

    if slice_:
        sinogram_error = np.sqrt(sinogram)
        sinogram_low = sinogram - sinogram_error
        sinogram_high = sinogram + sinogram_error

        image_low = reconstruct(sinogram_low, filter_=filter_)
        image_low = myprocess.denoise(image_low, 0.4)
        image_low = rotate(image_low, 5, reshape=False)

        image_high = reconstruct(sinogram_high, filter_=filter_)
        image_high = myprocess.denoise(image_high, 0.4)
        image_high = rotate(image_high, 5, reshape=False)

        p = Plot()
        p.shade(positions / PRE.m, sinogram[:, a_mid], sinogram_low[:, a_mid], sinogram_high[:, a_mid])
        p.show(xlabel='y (mm)', ylabel='Counts', grid=True, title='Sino')
        p = Plot()
        p.shade(positions / PRE.m, image[y_mid, :], image_low[y_mid, :], image_high[y_mid, :])
        p.show(xlabel='x (mm)', ylabel='Intensity (a.u.)', grid=True, legend=['Value', 'Error Bound'])
    else:
        p = Plot()
        p.cmap(sinogram, flip=True, label='Counts', extent=[np.rad2deg(min(angles)), np.rad2deg(max(angles)),
                                                            min(positions) / PRE.m, max(positions) / PRE.m])
        p.show(xlabel='$\\varphi$ (°)', ylabel='y (mm)')
        p = Plot()
        p.cmap(image, flip=True, label='Intensity (a.u.)', extent=[min(positions) / PRE.m, max(positions) / PRE.m] * 2)
        p.show(xlabel='x (mm)', ylabel='y (mm)')


def simulate(source_positions, y_limit=(-0.01, 0.01), slice_=True, filter_=True):
    profile = cv2.imread(str(ASSETS_DIR / 'flat-21.png'), cv2.IMREAD_GRAYSCALE) / 255
    tomo = Tomography()
    tomo.add_scint([SCINT_X, 0, 0], profile, width=2 * PRE.m, height=D_APERTURE)
    tomo.add_scint([-SCINT_X, 0, 0], profile, width=2 * PRE.m, height=D_APERTURE)
    for position in source_positions:
        tomo.add_source([*position, 0], activity=1E5)

    positions = np.linspace(-y_limit, y_limit, 101)
    angles = np.linspace(0, pi, 51)
    sinogram = tomo.sinogram(positions, angles)

    y_mid = len(positions) // 2
    a_mid = len(angles) // 2
    center_index = myprocess.index_peaks(sinogram[:, 0], distance=y_mid, prominence=0.5)[0]
    sinogram = np.roll(sinogram, y_mid - center_index, axis=0)
    image = reconstruct(sinogram, filter_=filter_)

    try:
        f = image[:, y_mid]
        delta = 1
        while f[y_mid] / 2 < f[
            y_mid + delta]:  # (f[y_mid + 2 * delta] + f[y_mid] - 2 * f[y_mid + delta]) <= 2 * np.sqrt(f[y_mid])
            delta += 1
        print(f'{positions[y_mid + delta] * 2 / PRE.m:.2f}', end=', ')
    except IndexError as e:
        pass

    if slice_:
        sinogram_error = np.sqrt(sinogram)
        sinogram_low = sinogram - sinogram_error
        sinogram_high = sinogram + sinogram_error
        image_low = reconstruct(sinogram_low, filter_=filter_)
        image_high = reconstruct(sinogram_high, filter_=filter_)

        p = Plot()
        p.shade(positions / PRE.m, sinogram[:, a_mid], sinogram_low[:, a_mid], sinogram_high[:, a_mid])
        p.show(xlabel='y (mm)', ylabel='Counts', grid=True, title='Sino')
        p = Plot()
        p.shade(positions / PRE.m, image[y_mid, :], image_low[y_mid, :], image_high[y_mid, :])
        p.show(xlabel='x (mm)', ylabel='Intensity (a.u.)', grid=True)
    else:
        p = Plot()
        p.cmap(sinogram, flip=True, label='Counts', extent=[np.rad2deg(min(angles)), np.rad2deg(max(angles)),
                                                            min(positions) / PRE.m, max(positions) / PRE.m])
        p.show(xlabel='$\\varphi$ (°)', ylabel='y (mm)')
        p = Plot()
        p.cmap(image, flip=True, label='Intensity (a.u.)', extent=[min(positions) / PRE.m, max(positions) / PRE.m] * 2)
        p.show(xlabel='x (mm)', ylabel='y (mm)', title='Reconstructed Image')

def fit_res():
    p = Plot()
    ws = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]
    fwhm = [1.6, 2.6, 2.0, 2.4, 2.9, 3.0, 3.5, 4.3, 4.2]
    fwhm_sim = [0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 1.0, 1.0, 1.1]

    p.errbar(ws, fwhm_sim, yerr=0.05, format_='ok', marker_size=4, z_order=100)
    p.errbar(ws, fwhm, yerr=0.2, format_='^b', z_order=100)

    ws.pop(1)
    fwhm_sim.pop(1)
    fwhm.pop(1)
    line = myprocess.Function.line
    params_sim, errors_sim, chi2_sim = myprocess.fit(line, ws, fwhm_sim, sigma=[0.05] * len(ws))
    params, errors, chi2 = myprocess.fit(line, ws, fwhm, sigma=[0.2] * len(ws))
    p.line(ws, line(ws, *params_sim), format_='k')
    p.line(ws, line(ws, *params), format_='--b')
    print(params_sim, errors_sim, chi2_sim)
    print(params, errors, chi2)

    p.show(xlabel='Slit width (mm)', ylabel='Peak FWHM (mm)', grid=True, legend=['Simulation', 'Experiment'])



if __name__ == '__main__':
    positions = [[0, 0]]  # single source
    # positions = [[-sep / 2, 0], [sep / 2, 0]]  # double source
    # positions = [[1, 2], [-1, 0], [-1, 1]]  # multi source

    # simulate(np.array(positions) * PRE.m, p, y_limit=10 * PRE.m, slice_=True, filter_=True, param=2, format_='k')
    # simulate(np.array(positions) * PRE.m, p, y_limit=10 * PRE.m, slice_=True, filter_=True, param=6, format_='--b')
    visualize_data(101, slice_=True, filter_=True)
    visualize_data(109, slice_=True, filter_=True)
