# Charge Transfer Distance and Transferred Charge Calculation

This repository provides a Python script to calculate the **Charge Transfer Distance (DCT)** and the **Transferred Charge (qCT)** from Gaussian cube files. It is designed to facilitate research in computational chemistry by automating the calculation of important charge transfer parameters from quantum chemical results.

## Features

- Computes **Charge Transfer Distance (DCT)** between ground and excited states.
- Computes **Transferred Charge (qCT)**.
- Outputs **centroids** of positive and negative density changes for further visualization.
- Supports parallel execution using **multithreading** for faster computation on large cube files.

## Dependencies

To run this script, you need the following dependencies:

- Python 3.7+
- NumPy
- `concurrent.futures` (included with Python 3)
- Gaussian software for generating cube files

You can install the Python dependencies using:

```bash
pip install numpy
