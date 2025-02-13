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

grid_info_file = os.path.join(args.out, "grid_info.txt")
if not os.path.exists(grid_info_file):
    raise FileNotFoundError(f"Arquivo {grid_info_file} não encontrado. Certifique-se de que o simulation salve o grid_info.txt.")
with open(grid_info_file, "r") as f:
    lx, ly = map(float, f.readline().split())
    nx, ny = map(int, f.readline().split())

x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
X, Y = np.meshgrid(x, y)

def plot_field(filename, title, label, cmap='viridis', output_name=None):
    data = np.loadtxt(os.path.join(args.out, filename))
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(X, Y, data.T, shading='auto', cmap=cmap)
    plt.colorbar(label=label)
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('Y')
    if output_name is None:
        output_name = f"{filename.split('.')[0]}_2d.png"
    plt.savefig(os.path.join(args.out, output_name))
    plt.close()

plot_field("density_e.txt", "Electron Density", "Electron Density")
plot_field("density_i.txt", "Ion Density", "Ion Density")

plot_field("phi_field.txt", "Electric Potential", "Potential (V)")

plot_field("electric_field_x.txt", "Electric Field (X Component)", "E_x (V/m)")
plot_field("electric_field_y.txt", "Electric Field (Y Component)", "E_y (V/m)")


vel_e = np.loadtxt(os.path.join(args.out, "velocity_e.txt"))
vel_i = np.loadtxt(os.path.join(args.out, "velocity_i.txt"))

def plot_velocity_histograms(vel_data, species_label, bins=100):
    components = ['vx', 'vy', 'vz']
    for i, comp in enumerate(components):
        plt.figure(figsize=(8,6))
        plt.hist(vel_data[:, i], bins=bins, density=True, alpha=0.7)
        plt.title(f"{species_label} Velocity Distribution ({comp})")
        plt.xlabel(f"{comp} (m/s)")
        plt.ylabel("Probability Density")
        plt.savefig(os.path.join(args.out, f"{species_label.lower()}_velocity_hist_{comp}.png"))
        plt.close()

plot_velocity_histograms(vel_e, "Electron")
plot_velocity_histograms(vel_i, "Ion")

plt.figure(figsize=(8,6))
plt.hist2d(vel_e[:, 0], vel_e[:, 1], bins=100, density=True, cmap='plasma')
plt.colorbar(label='Probability Density')
plt.title("Electron Velocity Distribution (vx vs vy)")
plt.xlabel("vx (m/s)")
plt.ylabel("vy (m/s)")
plt.savefig(os.path.join(args.out, "electron_velocity_hist2d_vx_vy.png"))
plt.close()

print(f"Todos os plots foram salvos em {args.out}")