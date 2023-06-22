import numpy as np
import os

# Set seed for reproducibility
np.random.seed(277)

# Outpath
#outpath = "/nfs/research/icortes/mkalyva/ClonalSimResults/MPN_Paper/ET1/Run_1_100_Sims"
outpath = "/Users/mkalyva/DATA/MSI_evolution_2022/MPN_trees/ET1"
obs_tree = "/Users/mkalyva/DATA/MSI_evolution_2022/MPN_trees/ET1/cancer_tree_tree_ET1.tree.nwk"
make_dir = f'mkdir -p {outpath}'
make_logs = f'mkdir -p {outpath}/logs'
os.system(make_dir)
os.system(make_logs)

# Function to run simulations
def run_simulation():
    # Draw s from a uniform distribution on (0, 2).
    s = np.random.uniform(1.9, 2)

    # Draw N from 10X, where X is uniformly distributed on (1, 9).
    X = np.random.uniform(1, 9)
    N = round(10 * X)

    # Draw L from round(Y), where Y is a Gaussian with mean 35 and std 5.
    L = round(np.random.normal(35, 5))
    while L < 2:  # redraw L until L >= 2
        L = round(np.random.normal(35, 5))

    # Draw g uniformly on the range 2, ..., L.
    g = round(np.random.randint(2, L+1))

    # Set k
    k = 22

    # Calculate mutation rate
    total_length_of_patient_tree = 723
    mutation_rate = total_length_of_patient_tree / (L - 1)

    # Set epsilon threshold
    epsilon_threshold = 0.0225

    # Printing the values
    print(f"s = {s}")
    print(f"N = {N}")
    print(f"L = {L}")
    print(f"g = {g}")
    print(f"k = {k}")
    print(f"mutation_rate = {mutation_rate}")
    #print(f"epsilon_threshold = {epsilon_threshold}")
    print(f"Outpath = {outpath}")


    # Submit the job
    bash_command = f'python /Users/mkalyva/GitHub/ColoSim/clonalSim.py --N {N} --generations {L} --mut_samples {k} --s {s} --mu {mutation_rate} --disease {g} --output_path {outpath} --observed_data_path {obs_tree}'
    os.system(bash_command)

    # Return parameters
    return s, N, L, g, k, mutation_rate

# Number of simulations
num_simulations = 4

# Run simulations
for sim_num in range(num_simulations):
    # Perform simulation...
    s, N, L, g, k, mutation_rate = run_simulation()

