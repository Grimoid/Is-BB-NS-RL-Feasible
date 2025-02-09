#!/bin/bash

# Load necessary modules

# Define the base parameters for xi_list
xi_list=(
    "0.3"
    "0.4"
    "0.5"
    "0.6"
    "0.7"
    "0.8"
)

# Define the values for T_list
T_list_values=(
    #"1000"
    #"2000"
    #"5000"
    #"10000"
    #"15000"
    #"20000"
    #"25000"
    #"30000"
    #"35000"
    #"40000"
    #"45000"
    #"50000"
    #"55000"
    #"60000"
    #"65000"
    #"70000"
    #"75000"
    #"80000"
    #"85000"
    #"90000"
    #"95000"
    "100000"
)

# Define the PB values
PB_list=(
    "Worst_Unif"
    "Normal_Unif"
     "Worst"
     "Normal"
)

# Ensure base output and error directories exist
mkdir -p ./out ./error

# Initialize an empty array to hold all the combined parameters
params=()

# Loop through each value of the T_list, PB_list, and xi_list, and combine them
for T_value in "${T_list_values[@]}"
do
    for PB in "${PB_list[@]}"
    do
        # Create output and error directories specific to T_value and PB
        mkdir -p ./out/out_T_${T_value}_${PB} ./error/error_T_${T_value}_${PB}

        for xi_value in "${xi_list[@]}"
        do
            params+=("$T_value $PB $xi_value")
        done
    done
done

# Calculate the number of instances based on the number of parameter sets
num_instances=${#params[@]}

# Base name for tmux sessions
session_base="experiment_session"

# Function to find the next available tmux session index
find_next_available_index() {
    index=0
    while tmux has-session -t "$index" 2>/dev/null; do
        index=$((index + 1))
    done
    echo $index
}

# Create and start tmux sessions
for param in "${params[@]}"
do
    T_value=$(echo $param | cut -d ' ' -f1)
    PB=$(echo $param | cut -d ' ' -f2)
    xi_value=$(echo $param | cut -d ' ' -f3)

    # Generate a unique program name
    program_name="program_${T_value}_${PB}_${xi_value}"

    # Compile the C++ files into the executable with the appropriate options
    g++ -O3 -march=native -fopenmp -std=c++11 parallel_exps.cpp algorithms.cpp problems.cpp globals.cpp BALG.cpp MASTER.cpp UCB.cpp -o $program_name 2> ./error/error_T_${T_value}_${PB}/compile_error_${xi_value}.err

    # Check if the compilation was successful
    if [ $? -ne 0 ]; then
        echo "Compilation failed for T=$T_value, PB=$PB, xi=$xi_value. Skipping this instance."
        continue
    fi

    # Prepare output and error file paths
    output_file="./out/out_T_${T_value}_${PB}/output_T${T_value}_xi${xi_value}.out"
    error_file="./error/error_T_${T_value}_${PB}/error_T${T_value}_xi${xi_value}.err"

    next_index=$(find_next_available_index)
    session_name="${next_index}"

    tmux new-session -d -s $session_name
    tmux send-keys -t $session_name "./$program_name $T_value $xi_value $PB > $output_file 2> $error_file" C-m
    echo "Session $session_name started with T=$T_value, xi=$xi_value, PB=$PB."
done

echo "$num_instances instances have been started in tmux sessions."
