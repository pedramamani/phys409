import os
import csv
import numpy as np
from pathlib import Path
from constants import PRE

DIR = Path(os.path.dirname(os.path.abspath(__file__)))
NAME = os.path.basename(__file__).split('.')[0]
ASSETS_DIR = DIR / f'{NAME}_assets'


def data(run_number):
    list_counts = []
    list_angle = []

    directory = ASSETS_DIR / f'run{run_number:03d}'
    file_names = [f for f in os.listdir(directory) if f.endswith('.dat')]
    for i, file_name in enumerate(file_names):
        name_chunks = file_name.split()
        angle = float(name_chunks[4].replace('_', '.'))

        with open(directory / file_name, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            list_y = []
            list_count = []

            reader.__next__()
            for row in reader:
                list_y.append(float(row[0]))
                list_count.append(float(row[1]))

        list_counts.append(list_count)
        list_angle.append(np.abs(angle))

    list_angle, list_counts = zip(*sorted(zip(list_angle, list_counts)))
    return np.array(list_counts).T, np.array(list_y) * PRE.m, np.deg2rad(list_angle)
