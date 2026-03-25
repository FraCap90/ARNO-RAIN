ARNO-RAIN Dataset – README

This document provides a brief overview of the dataset structure, file contents, and basic usage instructions.

Folder Structure

The dataset is organized into folders containing observed and simulated rainfall time series for each station, along with supporting metadata and parameter files.

1) stations_metadata.csv
Contains geographic and technical information for each station, including:
- Start and end year of record 
- Station ID
- Station name
- Longitude and latitude (decimal degrees)
- Elevation (m a.s.l.)
- Mean annual rainfall (mm)

2) IDF_parameters.csv
Provides the parameters of the regional IDF (Intensity–Duration–Frequency) curves used for calibrating the multifractal disaggregation model. Specifically, it includes:
- Scaling exponents (n)
- MRC parameters (Cβ, CLN, and D) resulting from model calibration

3) Time series folders (per station)
For each station, the dataset includes:
- Observed rainfall time series at 15-minute and daily resolution (2001–2020)
- Synthetic daily rainfall series (2000 years) generated using the CoSMoS-2s model
- Synthetic 15-minute rainfall series (2000 years) obtained via multifractal canonical disaggregation

File Formats
All files are provided in CSV format
Time series are organized with timestamps and corresponding rainfall values

Basic Usage
Use stations_metadata.csv to identify station characteristics and locations
Use IDF_parameters.csv for model calibration and reproducibility of the disaggregation procedure
Time series files can be directly imported into common analysis environments (e.g., R and Python)

Notes
The dataset is designed for benchmarking, evaluation, and intercomparison of stochastic rainfall generators and disaggregation methods
Synthetic series (2000 years) enable robust analysis of variability and extreme rainfall beyond the limits of observed records