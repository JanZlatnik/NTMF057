import os
import shutil
import subprocess
import toml
import numpy as np


def get_project_root():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(current_dir)


def create_and_run_test(run_number, rmax, points_per_au, exe_path, settings_path, output_dir):
    folder_name = os.path.join(output_dir, f"Run#{run_number}")
    os.makedirs(folder_name, exist_ok=True)

    with open(settings_path, "r") as f:
        settings = toml.load(f)

    settings["rmax"] = float(rmax)
    settings["nr"] = int(points_per_au * (rmax - settings["rmin"]))

    new_settings_path = os.path.join(folder_name, "settings.toml")
    with open(new_settings_path, "w") as f:
        toml.dump(settings, f)

    exe_name = os.path.basename(exe_path)
    shutil.copy2(exe_path, os.path.join(folder_name, exe_name))

    try:
        result = subprocess.run([os.path.join(folder_name, exe_name)],
                                cwd=folder_name,
                                check=True,
                                capture_output=True,
                                text=True)
        if result.returncode != 0:
            print(f"Error in Run#{run_number}:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Failed to run executable in Run#{run_number}:")
        print(e.stderr)
    except Exception as e:
        print(f"Unexpected error in Run#{run_number}:")
        print(str(e))


def main():
    project_root = get_project_root()

    exe_path = os.path.join(project_root, "target", "release", "Radial_problem.exe")
    settings_path = os.path.join(project_root, "settings.toml")
    output_dir = os.path.join(project_root, "Scripts", "ConvergenceTestNR")

    rmax_values = [12 for _ in range(9)]
    base_points_per_au = 100
    points_increase_factor = 2.0**(1/2)

    for i, rmax in enumerate(rmax_values, 1):
        points_per_au = base_points_per_au * (points_increase_factor ** i)
        print(f"\nCreating Run#{i}")
        print(f"rmax: {rmax:.1f} au")
        print(f"Points per au: {points_per_au:.0f}")
        create_and_run_test(i, rmax, points_per_au, exe_path, settings_path, output_dir)


if __name__ == "__main__":
    main()