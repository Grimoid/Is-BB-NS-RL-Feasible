import os
import numpy as np
import h5py

def read_results(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    regret_lines = []
    restart_lines = []
    avg_time = None
    stddev_time = None
    reading_regret = False
    reading_restarts = False

    for line in lines:
        if line.startswith("Regret:"):
            reading_regret = True
            reading_restarts = False
            continue
        elif line.startswith("Restarts:"):
            reading_regret = False
            reading_restarts = True
            continue
        elif line.startswith("MTPE:"):  # Mean time per experiment
            reading_regret=False
            reading_restarts=False
            avg_time = float(lines[lines.index(line) + 1].strip())
        elif line.startswith("STDEVPE:"):  # Standard deviation per experiment
            reading_regret=False
            reading_restarts=False
            stddev_time = float(lines[lines.index(line) + 1].strip())

        if reading_regret:
            regret_lines.append(line.strip().split())
        elif reading_restarts:
            restart_lines = line.strip().split()

    Regret = np.array(regret_lines, dtype=float)
    Restarts = np.array(restart_lines, dtype=int)

    # Process the file name for saving
    name = file_name.split('_')

    for i in range(len(name)):
        if name[i] == 'tweak':
            name[i+1] = str(int(float(name[i+1])))
        if name[i] == 'xi':
            name[i+1] = "0" + name[i+1].replace("0", "")
    
    name = "_".join(name)[:-4]  # Remove the .txt extension

    # Save the data to an HDF5 file
    with h5py.File(name, 'w') as hf:
        hf.create_dataset("Regret", data=Regret)
        hf.create_dataset("Restarts", data=Restarts)
        hf.create_dataset("AvgTime", data=avg_time)
        hf.create_dataset("StdDevTime", data=stddev_time)

    return Regret, Restarts, avg_time, stddev_time


def calculate_statistics(Regret, Restarts):
    regfinal = np.mean(Regret[:, -1])
    stdevfinal = np.std(Regret[:, -1])
    restarts = np.mean(Restarts)

    return regfinal, stdevfinal, restarts

def process_directory(directory, title="Normal_Geo", T=70000, N=4000, K=5):
    files = os.listdir(directory)
    txt_files = [f for f in files if f.endswith(".txt")]
    for file_name in txt_files:
        if title in file_name and f"T_{T}" in file_name and f"N_{N}" in file_name and f"K_{K}" in file_name:
            full_path = os.path.join(directory, file_name)
            print(f"Processing file: {full_path}")
            Regret, Restarts,avg_time, stddev_time = read_results(full_path)
            regfinal, stdevfinal, restarts = calculate_statistics(Regret, Restarts)
            print(f"File: {file_name}")
            print(f"Mean final regret is {regfinal}")
            print(f"Standard deviation of the final regret is {stdevfinal}")
            print(f"Mean number of restarts is {restarts}\n")
            print(f"Mean duration of one run {avg_time}\n")
            print(f"Standard deviation of one run {stddev_time}\n")

        else:print("File not found, change title, T, N, K.")


for Problem in ["Normal", "Worst","Normal_Unif","Worst_Unif"]:
    t_list=[1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000]
    N=4000
    for T in t_list:
        try:
            directory = f"./cpp_results_{Problem}_{T}_{N}/"
            process_directory(directory, title=Problem+"_Geo", T=T,N=N)
        except:print("Problem with: ", Problem, T)