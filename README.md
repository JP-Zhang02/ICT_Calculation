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
```

## Gaussian Cube File Preparation with `cubegen`

The script works with **Gaussian cube files** for ground and excited state electron densities. These files can be generated using **Gaussian software**. To generate cube files, you can use the following commands:

### Step 1: Run the Calculation

First, perform a density functional theory (DFT) or time-dependent DFT (TD-DFT) calculation using Gaussian to generate the checkpoint (.chk) file.

### Step 2: Generate Cube Files

Once you have the checkpoint file, use **`cubegen`** to generate the cube files for both the ground and excited states:

```bash
cubegen 0 density=scf ground_state.chk ground_state.cube -3 h
cubegen 0 density=transition excited_state.chk excited_state.cube -3 h
```

- `density=scf`: Specifies the type of density to extract (e.g., ground state).
- `density=transition`: Specifies the density for the excited state.
- `ground_state.chk` / `excited_state.chk`: The checkpoint files from Gaussian calculations.
- `ground_state.cube` / `excited_state.cube`: Output cube files.

Make sure that the cube files for both ground and excited states have **identical grid sizes** for successful calculation.

## Usage

The Python script uses command line arguments to specify input files and an optional output file. Here's how you can use it:

### Running the Script

```bash
python calculate_dct_qct.py <ground_state_file> <excited_state_file> [--output_file <output_file>]
```

### Example Command

```bash
python calculate_dct_qct.py ground_state.cube excited_state.cube --output_file results.txt
```

- `<ground_state_file>`: Path to the cube file representing the ground state density.
- `<excited_state_file>`: Path to the cube file representing the excited state density.
- `[--output_file <output_file>]`: (Optional) Path to the output file for storing results. If not provided, defaults to `results.txt`.

### Outputs

The script outputs:

1. **Charge Transfer Distance (DCT)**: The distance between the centroids of positive and negative charge densities.
2. **Transferred Charge (qCT)**: The total amount of charge transferred during the excitation.
3. **Centroids of Positive and Negative Density Changes** (`R+` and `R-`): Coordinates of positive and negative density centroids.

The results are also saved in the specified output file in a structured format for easy reference and later use.

### Example Output

```
=============================================
Results for files: 
Ground State File: ground_state.cube
Excited State File: excited_state.cube
Charge Transfer Distance (DCT): 12.34 Å
Transferred Charge (qCT): 0.78 e
Centroid of Positive Density Changes (R+): [2.34 1.23 0.98]
Centroid of Negative Density Changes (R-): [5.67 4.56 3.45]
=============================================
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### MIT License

```text
MIT License

Copyright (c) 2024 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Contributions

Contributions, bug reports, and feature requests are welcome! Please open an issue or submit a pull request.

## Acknowledgments

This script uses Gaussian cube files as input, which are generated by **Gaussian software**. We acknowledge the contributions of the computational chemistry community for developing and maintaining such powerful tools.

Feel free to reach out if you have questions or need further clarification regarding the usage of this script.
