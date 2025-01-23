import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(
    prog='plot_results_2d',
    description='Plot 2D benchmark and calculation results')
parser.add_argument("out",
                    help="Path to folder containing the calculation data.")
parser.add_argument("-d", "--data_path", help="Path to folder with benchmark data",
                    default="../data")
args = parser.parse_args()

# Read grid information
with open(f"{args.out}/grid_info.txt", "r") as f:
    lx, ly = map(float, f.readline().split())
    nx, ny = map(int, f.readline().split())

print(f"Tentando abrir: {os.path.join(args.out, 'field_e.txt')}")
# Read density data
de = np.loadtxt(os.path.join(args.out, "field_e.txt"))
di = np.loadtxt(os.path.join(args.out, "density_i.txt"))

# Create grid
x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
X, Y = np.meshgrid(x, y)

# Plot electron density
plt.figure(figsize=(10, 8))
plt.pcolormesh(X, Y, de.T, shading='auto')
plt.colorbar(label='Electron Density')
plt.title('Electron Density')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig(f"{args.out}/electron_density_2d.png")
plt.close()

# Plot ion density
plt.figure(figsize=(10, 8))
plt.pcolormesh(X, Y, di.T, shading='auto')
plt.colorbar(label='Ion Density')
plt.title('Ion Density')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig(f"{args.out}/ion_density_2d.png")
plt.close()

print(f"Plots saved in {args.out}")