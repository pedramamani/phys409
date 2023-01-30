# Positron Emission Tomography

## Which run is which?

The sheet in [this file](https://docs.google.com/spreadsheets/d/13XCm0C1i5AAMUfErssmMvvk31R0u61NUlUm9d0b30Ho/edit?usp=sharing) provides information on the various data runs.

## Format of individual data files

Data is collected in a file with the following naming format:

[DATE], [TIME], [SCAN#], [ANGLE IN DEGREES] Deg.dat

For example:

Aug-01-20, 1_30 PM, Scan1, -50_4 Deg.dat

The data within these files is of the following format:

Distance (mm)	Counts
0.0	130.0
0.2	178.0
0.4	170.0
0.6	176.0
0.8	184.0
1.0	167.0
.
.
.

where the distance is the lateral distance that the source has been moved.
