def plot_dynamic_regret(PB_NAME, T,N):
    import numpy as np
    import matplotlib.pyplot as plt
    import h5py
    import seaborn as sns  # For accessing colorblind palette
    import matplotlib
    import matplotlib.ticker as ticker

    # Close any previous plots
    plt.close('all')

    # Import seaborn and set style
    sns.set_style("whitegrid")

    # Set default font properties after seaborn
    matplotlib.rcParams['font.size'] = 12
    matplotlib.rcParams['axes.titlesize'] = 28
    matplotlib.rcParams['axes.labelsize'] = 23
    matplotlib.rcParams['xtick.labelsize'] = 16
    matplotlib.rcParams['ytick.labelsize'] = 17
    matplotlib.rcParams['legend.fontsize'] = 21
    matplotlib.rcParams['figure.titlesize'] = 32
    matplotlib.rcParams['axes.titleweight'] = 'bold'
    matplotlib.rcParams['axes.labelweight'] = 'bold'
    matplotlib.rcParams['figure.titleweight'] = 'bold'

    # Values of xi to iterate over
    xi_values = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    # Problem parameters
    PB_values = [PB_NAME]  # List of problems
    K = 5
    tweak = 1
    beta_type = 'beta'

    # Policies and styles
    policies = [
        "MASTER_ALGO",
        "RANDALG_p_0_05",
        "RR_OPT",
        "QCD_PLUS_KLUCB",
        "QCD_PLUS_UCB",
        "GLRKLUCB"
    ]

    colors = ['magenta', 'blue', 'red', 'limegreen', 'black', 'darkviolet']
    line_styles = ['-', '--', '-.', 'solid', (0, (3, 1, 1, 1)), (0, (1, 1))]
    markers = ['o', 's', 'D', '^', 'v', 'x']
    final_names = ["MASTER", "RR_p.05", "RR", "QCD+klUCB", "QCD+UCB", "GLRklUCB"]

    lP = len(policies)

    # Create a custom ScalarFormatter subclass
    from matplotlib.ticker import ScalarFormatter

    class FixedScalarFormatter(ScalarFormatter):
        def _set_format(self):
            self.format = '%.1f'  # Force format to one decimal place

    # Create the figure and axes for subplots
    fig, axs = plt.subplots(2, 3, figsize=(19, 9))
    axs = axs.flatten()

    # Loop over problems and xi values
    for iPB, PB in enumerate(PB_values):
        for idx, xi in enumerate(xi_values):

            ax = axs[idx]

            # Plot each policy for this value of xi and PB
            for imeth in range(lP):
                policy = policies[imeth]
                fname = f"cpp_results_{PB}_{T}_{N}/" + PB + "_Geo_" + str(beta_type) + "_tweak_" + str(tweak)
                name = fname + f"_{policy}_T_{T}_N_{N}_K_{K}_xi_{xi}"

                # Load the data
                try:
                    with h5py.File(name, 'r') as hf:
                        Regret = np.array(hf.get('Regret'))
                except FileNotFoundError:
                    print(f"File not found: {name}")
                    continue

                tsave = np.floor(np.linspace(1, T, num=500))
                reg_mean = np.mean(Regret, axis=0)
                reg_std = np.std(Regret, axis=0)

                # Plot with shaded area for standard deviation
                ax.fill_between(tsave, reg_mean - reg_std, reg_mean + reg_std, alpha=0.2, color=colors[imeth])

                # Add label only once
                if idx == 0:
                    ax.plot(tsave, reg_mean, label=final_names[imeth], color=colors[imeth],
                            linestyle=line_styles[imeth], linewidth=2)
                else:
                    ax.plot(tsave, reg_mean, color=colors[imeth],
                            linestyle=line_styles[imeth], linewidth=2)

            # Style the subplot
            if PB == "Normal":
                ax.set_title(fr"$\xi$ = {xi}")
            elif PB == "Normal_Unif":
                ax.set_title(fr"$N_C$ = {int(np.ceil(T ** (1 - xi)))}")
            elif PB == "Worst":
                ax.set_title(fr"$\xi$ = {xi}")
            elif PB == "Worst_Unif":
                ax.set_title(fr"$N_C$ = {int(np.ceil(T ** (1 - xi)))}")
            else:
                # For other PB values, default title
                ax.set_title(fr"$\xi$ = {xi}")
            ax.grid(alpha=0.3)
            if idx > 2:
                ax.set_xlabel("Environment Steps")
            if idx == 0 or idx == 3:
                ax.set_ylabel("Dynamic Regret")

            # Set the x-axis limits to [0, T]
            ax.set_xlim([0, T])

            # Use ScalarFormatter with scientific notation for the x-axis
            formatter_x = ScalarFormatter(useMathText=True)
            formatter_x.set_powerlimits((0, 0))
            ax.xaxis.set_major_formatter(formatter_x)
            ax.xaxis.get_offset_text().set_fontsize(13)

            # Use custom FixedScalarFormatter for the y-axis
            formatter_y = FixedScalarFormatter(useMathText=True)
            formatter_y.set_powerlimits((0, 0))
            ax.yaxis.set_major_formatter(formatter_y)
            ax.yaxis.get_offset_text().set_fontsize(13)

            # Optionally, adjust tick parameters
            ax.tick_params(axis='both', which='major')

    # Create unified title
    if PB_NAME == "Normal":
        plt.suptitle(f"Results for Geometric Uniform, T={T}, {N} Runs")
    if PB_NAME == "Normal_Unif":
        plt.suptitle(f"Results for Deterministic Uniform, T={T}, {N} Runs")
    if PB_NAME == "Worst":
        plt.suptitle(f"Results for Geometric Worst, T={T}, {N} Runs")
    if PB_NAME == "Worst_Unif":
        plt.suptitle(f"Results for Deterministic Worst, T={T}, {N} Runs")

    # Collect all handles and labels from subplots
    handles_labels = [ax.get_legend_handles_labels() for ax in axs.flat]
    handles, labels = zip(*handles_labels)
    handles = sum(handles, [])
    labels = sum(labels, [])
    from collections import OrderedDict
    by_label = OrderedDict(zip(labels, handles))

    # Create a single legend at the bottom
    legend = fig.legend(by_label.values(), by_label.keys(), loc='lower center', ncol=len(policies),
                        prop={'weight': 'bold'})

    # Update the linewidth for the lines in the legend
    for line in legend.get_lines():
        line.set_linewidth(4)  # Increase linewidth of legend lines

    # Adjust the layout to fit everything cleanly
    plt.tight_layout(rect=[0, 0.00, 1, 1])
    plt.subplots_adjust(bottom=0.16, top=0.87, hspace=0.34)

    # Save the plot
    plt.savefig(f"D_Reg_comp_{PB_NAME}_{T}_{N}_suppl.pdf", dpi=300)
    plt.show()

    #plt.show(block=False)
    #plt.pause(0.01)
    #plt.close('all')
N=4000
for PB in ["Normal", "Normal_Unif", "Worst", "Worst_Unif"]:
    for T in [1000,2000]+list(range(5000, 105000, 5000)):
        plot_dynamic_regret(PB, T, N)