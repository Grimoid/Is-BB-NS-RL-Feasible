def plot_final_regret_vs_xi(problems, T, N, updateCSV="on"):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import h5py
    from collections import OrderedDict
    import seaborn as sns
    import matplotlib
    import matplotlib.ticker as ticker
    import math  # Import math module for ceiling function

    # Close any previous plots
    plt.close('all')

    # Set seaborn style
    sns.set_style("whitegrid")

    # Set default font properties after seaborn
    matplotlib.rcParams['font.size'] = 12
    matplotlib.rcParams['axes.titlesize'] = 23
    matplotlib.rcParams['axes.labelsize'] = 19
    matplotlib.rcParams['xtick.labelsize'] = 15
    matplotlib.rcParams['ytick.labelsize'] = 15
    matplotlib.rcParams['legend.fontsize'] = 17
    matplotlib.rcParams['figure.titlesize'] = 26
    matplotlib.rcParams['axes.titleweight'] = 'bold'
    matplotlib.rcParams['axes.labelweight'] = 'bold'
    matplotlib.rcParams['figure.titleweight'] = 'bold'

    # Constants
    tweak = 1
    K = 5
    beta_type = 'beta'

    # Policies and names
    policies = [
        "MASTER_ALGO",
        "RANDALG_p_0_05",
        "RR_OPT",
        "QCD_PLUS_KLUCB",
        "QCD_PLUS_UCB",
        "GLRKLUCB"
    ]

    final_names = ["MASTER", "RR_p.05", "RR", "QCD+klUCB", "QCD+UCB", "GLRklUCB"]

    # Use specific colors, line styles, and markers
    colors = ['magenta', 'blue', 'red', 'limegreen', 'black', 'darkviolet']
    line_styles = ['-', '--', '-.', 'solid', (0, (3, 1, 1, 1)), (0, (1, 1))]
    markers = ['o', 's', 'D', '^', 'v', 'x']

    # Initialize xi_values
    xi_values = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    # Create the figure and axes for subplots
    fig_rows = 2
    fig_cols = 2
    fig, axs = plt.subplots(fig_rows, fig_cols, figsize=(13, 9))
    axs = axs.flatten()

    # Initialize a dictionary to hold data for all problems
    data_all = {}

    # Loop over problems to create subplots
    for idx, PB in enumerate(problems):
        # Base filename for the current problem
        fname = f"cpp_results_{PB}_{T}_{N}/{PB}_Geo_{beta_type}_tweak_{tweak}"

        # Initialize dictionary to hold data for the current problem
        data = {policy: {'xi': [], 'regfinal': [], 'stdevfinal': []} for policy in policies}

        # Loop over xi_values
        for xi in xi_values:
            xi = round(xi, 2)  # Ensure xi is rounded appropriately
            for policy in policies:
                # Construct filename
                name = fname + f"_{policy}_T_{T}_N_{N}_K_{K}_xi_{xi}"

                # Try to open the file
                try:
                    with h5py.File(name, 'r') as hf:
                        Regret = np.array(hf.get('Regret'))
                except FileNotFoundError:
                    print(f"File {name} not found.")
                    continue  # Skip to next iteration if file not found

                # Compute regfinal and stdevfinal
                regfinal = np.mean(Regret[:, -1])
                stdevfinal = np.std(Regret[:, -1])

                # Append data
                data[policy]['xi'].append(xi)
                data[policy]['regfinal'].append(regfinal)
                data[policy]['stdevfinal'].append(stdevfinal)

        # Store data for the current problem
        data_all[PB] = data

        # Plot regfinal vs xi for each policy in the current subplot
        ax = axs[idx]

        for idx_policy, policy in enumerate(policies):
            xi_list = data[policy]['xi']
            reg_list = data[policy]['regfinal']
            std_list = data[policy]['stdevfinal']

            if len(xi_list) == 0:
                continue  # Skip if no data for this policy

            # Sort the data for plotting
            xi_list, reg_list, std_list = zip(*sorted(zip(xi_list, reg_list, std_list)))

            xi_list = list(xi_list)
            reg_list = list(reg_list)
            std_list = list(std_list)

            # Plot with error bars (still plotting vs xi)
            ax.errorbar(xi_list, reg_list, yerr=std_list, label=final_names[idx_policy], color=colors[idx_policy],
                        linestyle=line_styles[idx_policy], marker='o', capsize=3, markersize=6, linewidth=2)

        # Set axis labels and title
        if PB in ['Normal_Unif', 'Worst_Unif']:
            # Calculate N_c for each xi
            N_c_values = [int(math.ceil(T ** (1 - xi))) for xi in xi_values]
            # Set x-axis label to 'Number of Change Points'
            ax.set_xlabel(r'Number of Change Points')
            # Set x-ticks to xi_values
            ax.set_xticks(xi_values)
            # Map xi_values to N_c
            xi_labels = [str(N_c) for N_c in N_c_values]
            # Set the x-tick labels to the N_c values
            ax.set_xticklabels(xi_labels)
        else:
            ax.set_xlabel(r'$\xi$')

        if idx % fig_cols == 0:
            ax.set_ylabel('Final Dynamic Regret')  # Set y-axis label only on the first column for clarity

        if PB == 'Normal':
            ax.set_title('Geometric Uniform')
        elif PB == 'Worst':
            ax.set_title('Geometric Worst')
        elif PB == 'Normal_Unif':
            ax.set_title('Deterministic Uniform')
        elif PB == 'Worst_Unif':
            ax.set_title('Deterministic Worst')

        # Format y-axis with scientific notation if needed
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((0, 0))
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.get_offset_text().set_fontsize(13)

        # Adjust tick parameters
        ax.tick_params(axis='both', which='major')

        # Add grid
        ax.grid(alpha=0.3)

    # Remove any unused subplots if problems < number of subplots
    if len(problems) < len(axs):
        for idx in range(len(problems), len(axs)):
            fig.delaxes(axs[idx])

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0.00, 1, 1])
    plt.subplots_adjust(bottom=0.135, top=0.87, hspace=0.35)


    # Create a single legend for all subplots
    # Collect all handles and labels from subplots with data
    handles_labels = [ax.get_legend_handles_labels() for ax in axs[:len(problems)]]
    handles, labels = zip(*handles_labels)
    handles = sum(handles, [])
    labels = sum(labels, [])
    from collections import OrderedDict
    by_label = OrderedDict(zip(labels, handles))

    # Create a single legend at the bottom with frame (legend box)
    legend = fig.legend(by_label.values(), by_label.keys(), loc='lower center', ncol=len(policies),
                        prop={'weight': 'bold'})

    # Update the linewidth for the lines in the legend
    for line in legend.get_lines():
        line.set_linewidth(4)  # Increase linewidth of legend lines

    # Set the suptitle
    plt.suptitle(rf'Results for Horizon T={T}, {N} Runs')

    # Save the plot
    plt.savefig(f'Robustness_T_{T}_N_{N}_suppl.pdf', dpi=300, bbox_inches='tight')
    plt.show(block=False)
    plt.pause(0.01)
    plt.close('all')

    # Optionally, save the results to a CSV file
    if updateCSV == "on":
        df_list = []
        for PB in problems:
            data = data_all[PB]
            for idx_policy, policy in enumerate(policies):
                xi_list = data[policy]['xi']
                reg_list = data[policy]['regfinal']
                std_list = data[policy]['stdevfinal']
                if len(xi_list) == 0:
                    continue
                df_temp = pd.DataFrame({
                    'Problem': PB,
                    'Policy': final_names[idx_policy],
                    'xi': xi_list,
                    'Regret': reg_list,
                    'StdDev': std_list
                })
                df_list.append(df_temp)
        df_all = pd.concat(df_list, ignore_index=True)
        df_all.to_csv(f"Regret_vs_xi_All_Problems_T_{T}_N_{N}.csv", index=False)


# Define the list of problems
problems_list = ["Normal", "Worst", "Normal_Unif", "Worst_Unif"]

# Call the function with specific T and N values
for T in [1000,2000]+list(range(5000, 105000, 5000)):
    plot_final_regret_vs_xi(problems=problems_list, T=T, N=4000)