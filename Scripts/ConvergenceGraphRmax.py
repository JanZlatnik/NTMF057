import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import os
import toml

title = r'Grid length convergence test'

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 10
SIZE = 6

rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Latin Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = ["Latin Modern Roman"]
plt.rcParams['axes.titlepad'] = 10
plt.rcParams['axes.labelpad'] = 10
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=MEDIUM_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=LEGEND_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)

def load_energies(base_path):
    energies_dict = {}
    settings_dict = {}
    run_folders = sorted([f for f in os.listdir(base_path) if f.startswith("Run#")],
                        key=lambda x: int(x.split("#")[1]))

    for folder in run_folders:
        energy_file = os.path.join(base_path, folder, "Task 1", "energies.txt")
        settings_file = os.path.join(base_path, folder, "settings.toml")
        try:
            data = np.loadtxt(energy_file, comments='#', usecols=(0, 1))
            energies_dict[folder] = data

            if os.path.exists(settings_file):
                settings = toml.load(settings_file)
                settings_dict[folder] = settings
            print(f"Loaded {folder}: {len(data)} energies")
        except Exception as e:
            print(f"Error loading {energy_file}: {str(e)}")

    return energies_dict, settings_dict

def calculate_differences(energies_dict):
    last_run = max(energies_dict.keys(), key=lambda x: int(x.split("#")[1]))
    reference_energies = energies_dict[last_run]

    differences = {}
    for run, energies in energies_dict.items():
        if run == last_run:
            continue

        min_length = min(len(energies), len(reference_energies))

        diff = np.abs(energies[:min_length, 1] - reference_energies[:min_length, 1])
        differences[run] = (energies[:min_length, 0], diff)
        print(f"Comparing {run} with {last_run}: using {min_length} energies")

    return differences, last_run

def plot_convergence(differences, last_run, settings_dict, legend_key="rmax", output_name="convergence"):
    fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    xmax = settings_dict[last_run][legend_key] if legend_key in settings_dict[last_run] else "N/A"
    legend_title = rf"$r_\mathrm{{max}} = {xmax}\,\mathrm{{au}}$"

    for run, (x, diff) in differences.items():
        run_number = int(run.split("#")[1])
        value = settings_dict[run].get(legend_key, "N/A")
        label = rf"${value}\,\mathrm{{au}}$"
        ax.semilogy(x, diff, '-', label=label)

    ax.set_xlabel(r'$n$')
    ax.set_ylabel(r'$|\Delta E_n|\,(\mathrm{a.u.})$')
    ax.set_title(title)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(False)
    ax.legend(title=legend_title, facecolor='None', edgecolor='None',ncol=2)

    fig.savefig(f'{output_name}.pdf', format='pdf', bbox_inches='tight')
    plt.show()

def main():
    convergence_test_path = r"D:\\Univerzita\\Numerické metody\\Radial_problem\\Scripts\\ConvergenceTestRmax"
    legend_key = "rmax"

    print("Loading energy data...")
    energies_dict, settings_dict = load_energies(convergence_test_path)

    print("\nCalculating energy differences...")
    differences, last_run = calculate_differences(energies_dict)

    print("\nGenerating plot...")
    plot_convergence(differences, last_run, settings_dict, legend_key, "rmax_graph")

if __name__ == "__main__":
    main()
