import numpy as np
import logging
import argparse
from concurrent.futures import ThreadPoolExecutor

# 配置日志，记录详细的信息到文件 'dct_qct_calculation.log' 中
logging.basicConfig(filename='dct_qct_calculation.log',
                    filemode='w',
                    level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Function to read a Gaussian cube file and extract the origin, grid size, grid spacing, and density data.
def read_cube_file(filename):
    """
    Reads a cube file and extracts the density grid data along with grid dimensions and spacings.
    The cube file contains electron density information.

    Parameters:
    filename (str): Path to the cube file.

    Returns:
    origin (np.array): The origin of the grid (in Cartesian coordinates).
    grid_size (np.array): The number of grid points along each axis.
    grid_spacing (np.array): The grid spacing in angstroms along each axis.
    density_data (np.array): The electron density values on the grid.
    """
    logging.info(f"Reading cube file: {filename}")
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"Cube file {filename} not found.")
        raise

    # The 3rd line contains information about the number of atoms and the origin of the grid.
    atom_count_line = lines[2].split()
    num_atoms = int(atom_count_line[0])  # Number of atoms in the system (not used here)
    origin = np.array([float(x) for x in atom_count_line[1:]])

    logging.debug(f"Number of atoms: {num_atoms}")
    logging.debug(f"Origin of the grid: {origin}")

    # The next three lines contain the grid size and grid spacing along each axis.
    grid_size = []
    grid_spacing = []
    for i in range(3):
        grid_line = lines[3 + i].split()
        grid_size.append(int(grid_line[0]))  # The number of grid points along this axis
        grid_spacing.append([float(x) for x in grid_line[1:]])
        logging.debug(f"Grid size along axis {i}: {grid_size[-1]}, Grid spacing: {grid_spacing[-1]}")

    grid_size = np.array(grid_size)  # Convert grid size to NumPy array
    grid_spacing = np.array(grid_spacing)  # Convert grid spacing to NumPy array

    # The remaining lines contain the electron density values.
    data_start = 6 + num_atoms  # Skip header and atomic positions
    density_data = []
    for line in lines[data_start:]:
        density_data.extend([float(x) for x in line.split()])

    # Convert the list of density values into a NumPy array and reshape it to match the grid dimensions.
    density_data = np.array(density_data).reshape(grid_size)

    logging.info(f"Successfully read density data from cube file: {filename}")
    logging.debug(f"Density data shape: {density_data.shape}")
    
    return origin, grid_size, grid_spacing, density_data

# Function to compute the density difference Δρ(r) between ground and excited states.
def compute_density_difference(ground_density, excited_density):
    """
    Computes the difference between the excited state and ground state densities.
    This represents the charge transfer upon excitation.

    Parameters:
    ground_density (np.array): The electron density of the ground state.
    excited_density (np.array): The electron density of the excited state.

    Returns:
    delta_rho_plus (np.array): The positive density changes (density increase).
    delta_rho_minus (np.array): The negative density changes (density depletion).
    """
    # Compute the density difference Δρ(r)
    delta_rho = excited_density - ground_density
    
    # Split the density into positive (density increase) and negative (density depletion) parts
    delta_rho_plus = np.where(delta_rho > 0, delta_rho, 0)  # Positive density changes
    delta_rho_minus = np.where(delta_rho < 0, delta_rho, 0)  # Negative density changes
    
    logging.info("Computed density differences between ground and excited states.")
    logging.debug(f"Delta rho (positive): Sum = {np.sum(delta_rho_plus):.6f}")
    logging.debug(f"Delta rho (negative): Sum = {np.sum(delta_rho_minus):.6f}")
    
    return delta_rho_plus, delta_rho_minus

# Function to calculate the centroid of a density distribution.
def calculate_centroid(density, origin, grid_size, grid_spacing):
    """
    Calculates the centroid (center of mass) of the given electron density distribution.

    Parameters:
    density (np.array): The electron density distribution on the grid.
    origin (np.array): The origin of the grid.
    grid_size (np.array): The number of grid points along each axis.
    grid_spacing (np.array): The grid spacing along each axis.

    Returns:
    centroid (np.array): The centroid coordinates in angstroms.
    """
    total_density = np.sum(density)
    
    if total_density == 0:
        logging.warning("Total density is zero. Returning zero vector for centroid.")
        return np.zeros(3)

    # Create a 3D grid of indices that represents the position of each grid point
    grid_indices = np.indices(grid_size)  # shape: (3, nx, ny, nz)

    # Initialize the grid positions array
    grid_positions = np.zeros((3,) + tuple(grid_size))  # shape: (3, nx, ny, nz)

    # Calculate the Cartesian coordinates of each grid point based on the origin and spacing
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for i in range(3):  # x, y, z axes
            futures.append(executor.submit(lambda idx: origin[idx] + grid_indices[idx] * grid_spacing[idx][idx], i))

        for i, future in enumerate(futures):
            grid_positions[i] = future.result()

    # Calculate the centroid for each axis by weighting the positions by the density values
    centroid = np.zeros(3)
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for i in range(3):
            futures.append(executor.submit(lambda idx: np.sum(grid_positions[idx] * density) / total_density, i))

        for i, future in enumerate(futures):
            centroid[i] = future.result()

    logging.debug(f"Calculated centroid: {centroid}")
    return centroid

# Function to compute the charge transfer distance (DCT) and transferred charge (qCT).
def compute_dct_and_qct(delta_rho_plus, delta_rho_minus, origin, grid_size, grid_spacing):
    """
    Computes the charge transfer distance (DCT) and the amount of transferred charge (qCT).

    Parameters:
    delta_rho_plus (np.array): The positive density changes (density increase).
    delta_rho_minus (np.array): The negative density changes (density depletion).
    origin (np.array): The origin of the grid.
    grid_size (np.array): The number of grid points along each axis.
    grid_spacing (np.array): The grid spacing along each axis.

    Returns:
    dct (float): The charge transfer distance in angstroms.
    qct (float): The amount of charge transferred (in electron charge units).
    centroid_plus (np.array): The centroid coordinates of positive density changes.
    centroid_minus (np.array): The centroid coordinates of negative density changes.
    """
    centroid_plus = calculate_centroid(delta_rho_plus, origin, grid_size, grid_spacing)
    centroid_minus = calculate_centroid(delta_rho_minus, origin, grid_size, grid_spacing)
    
    dct = np.linalg.norm(centroid_plus - centroid_minus)
    qct = np.sum(delta_rho_plus)  # Assuming charge conservation, qCT should be same for + and -

    logging.info("Computed DCT and qCT.")
    logging.debug(f"DCT (Charge Transfer Distance): {dct:.6f} Å")
    logging.debug(f"qCT (Transferred Charge): {qct:.6f} e")
    logging.debug(f"Centroid (R+): {centroid_plus}")
    logging.debug(f"Centroid (R-): {centroid_minus}")
    
    return dct, qct, centroid_plus, centroid_minus

# Main function to execute the calculations.
def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Calculate Charge Transfer Distance (DCT) and Transferred Charge (qCT) from Gaussian cube files.")
    parser.add_argument("ground_state_file", help="Path to the cube file representing the ground state density")
    parser.add_argument("excited_state_file", help="Path to the cube file representing the excited state density")
    parser.add_argument("--output_file", default="results.txt", help="Path to the output file for storing results (default: results.txt)")

    args = parser.parse_args()

    # 获取命令行参数
    ground_state_file = args.ground_state_file
    excited_state_file = args.excited_state_file
    result_file = args.output_file

    logging.info("Starting DCT and qCT calculation.")

    # Load the ground state and excited state density data from cube files
    try:
        origin, grid_size, grid_spacing, ground_density = read_cube_file(ground_state_file)
        _, _, _, excited_density = read_cube_file(excited_state_file)
    except Exception as e:
        logging.error(f"Error while reading cube files: {e}")
        print(f"Error while reading cube files: {e}")
        return

    # Verify if the grid sizes match for both cube files
    if ground_density.shape != excited_density.shape:
        logging.error("Mismatch in grid size between ground state and excited state cube files.")
        print("Mismatch in grid size between ground state and excited state cube files.")
        return

    # Compute the density differences Δρ+ and Δρ-
    delta_rho_plus, delta_rho_minus = compute_density_difference(ground_density, excited_density)

    # Compute DCT (charge transfer distance), qCT (transferred charge), and centroids (R+ and R-)
    dct, qct, centroid_plus, centroid_minus = compute_dct_and_qct(delta_rho_plus, delta_rho_minus, origin, grid_size, grid_spacing)

    # Print the results and log them
    print(f"Charge Transfer Distance (DCT): {dct:.2f} Å")
    print(f"Transferred Charge (qCT): {qct:.2f} e")
    print(f"Centroid of Positive Density Changes (R+): {centroid_plus}")
    print(f"Centroid of Negative Density Changes (R-): {centroid_minus}")

    logging.info(f"Charge Transfer Distance (DCT): {dct:.2f} Å")
    logging.info(f"Transferred Charge (qCT): {qct:.2f} e")
    logging.info(f"Centroid of Positive Density Changes (R+): {centroid_plus}")
    logging.info(f"Centroid of Negative Density Changes (R-): {centroid_minus}")

    # Write results to the result file
    try:
        with open(result_file, 'a') as f:
            f.write("=============================================\n")
            f.write(f"Results for files: \nGround State File: {ground_state_file}\nExcited State File: {excited_state_file}\n")
            f.write(f"Charge Transfer Distance (DCT): {dct:.2f} Å\n")
            f.write(f"Transferred Charge (qCT): {qct:.2f} e\n")
            f.write(f"Centroid of Positive Density Changes (R+): {centroid_plus}\n")
            f.write(f"Centroid of Negative Density Changes (R-): {centroid_minus}\n")
            f.write("=============================================\n\n")
        print(f"Results successfully saved to {result_file}")
    except Exception as e:
        logging.error(f"Error while writing results to file: {e}")
        print(f"Error while writing results to file: {e}")

# Entry point of the script
main()
