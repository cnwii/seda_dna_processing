#!/usr/bin/env python

# Parameters
N_bp=25      # Number of terminal bases to analyze
TRIM=1       # Bases to trim before analysis

# Function to count A/T/C/G at 5' and 3' ends
Count_ATCG() {
  FILE=$1
  N_bp=$2
  TRIM=$3

  awk -v n_bp="$N_bp" -v trim="$TRIM" \
      -v out5="${FILE/.fastq}_5_end_freq" \
      -v out3="${FILE/.fastq}_3_end_freq" '
    BEGIN {
      print "Position_from_5end\tA_freq\tT_freq\tC_freq\tG_freq\tTotal" > out5
      print "Position_from_3end\tA_freq\tT_freq\tC_freq\tG_freq\tTotal" > out3
    }
    NR % 4 == 2 {
      seq = $0
      len = length(seq)
      for (i = trim + 1; i <= n_bp; i++) {
        base5 = substr(seq, i, 1)
        base3 = substr(seq, len - i + 1, 1)
        count5[i][base5]++
        count3[i][base3]++
        total5[i]++
        total3[i]++
      }
    }
    END {
      for (i = trim + 1; i <= n_bp; i++) {
        printf "%d\t%d\t%d\t%d\t%d\t%d\n", i, count5[i]["A"]+0, count5[i]["T"]+0, count5[i]["C"]+0, count5[i]["G"]+0, total5[i]+0 >> out5
        printf "%d\t%d\t%d\t%d\t%d\t%d\n", i, count3[i]["A"]+0, count3[i]["T"]+0, count3[i]["C"]+0, count3[i]["G"]+0, total3[i]+0 >> out3
      }
    }
  ' "$FILE"
}

# Run function on all fastq files in directory
for FILE in *.fastq; do
  echo "Processing $FILE..."
  Count_ATCG "$FILE" "$N_bp" "$TRIM"
done

echo "Calculate nucleotide frequencies complete."


import pandas as pd

# Load the frequency tables for 3' and 5' end damage
file_3_end = "ERR5729614_filtered_25bp_100k_3_end_freq"
file_5_end = "ERR5729614_filtered_25bp_100k_5_end_freq"

# Try reading with tab separator
df_3 = pd.read_csv(file_3_end, sep="\t")
df_5 = pd.read_csv(file_5_end, sep="\t")

# Show the first five rows 
df_3.head(), df_5.head()

# binomtest() performs a binomial hypothesis test
# in this case, it tests whether the number of "successes" (i.e., DNA base calls) in a fixed number of independent Bernoulli trials is consistent with a specified probability of success 
# C-->T transitions at reads ends are enriched and G-->A transitions on the reverse strand are elevated
from scipy.stats import binomtest


import matplotlib.pyplot as plt

# Subset the 5' end data to positions 2 to 25
df_5_plot = df_5[(df_5['Position_from_5end'] >= 2) & (df_5['Position_from_5end'] <= 25)].copy()

# Subset the 3' end data to positions 2 to 25
df_3_plot = df_3[(df_3['Position_from_3end'] >= 2) & (df_3['Position_from_3end'] <= 25)].copy()

# Create side-by-side plots like in mapDamage (smiley plot layout)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# 5' end plot: T_freq vs position
ax1.plot(df_5_plot['Position_from_5end'], df_5_plot['T_freq'], marker='o', color='firebrick')
ax1.set_xlabel("Position from 5′ end")
ax1.set_ylabel("Frequency")
ax1.set_title("5′ End (C→T transitions)")
ax1.grid(True)

# 3' end plot: A_freq vs position (in reversed x to mirror orientation)
ax2.plot(df_3_plot['Position_from_3end'], df_3_plot['A_freq'], marker='o', color='darkblue')
ax2.set_xlabel("Position from 3′ end")
ax2.set_title("3′ End (G→A transitions)")
ax2.grid(True)

# Flip the x-axis of the 3' end so both plots mirror inward (smiley layout)
ax2.invert_xaxis()

plt.suptitle("Smiley Plot: 5′ and 3′ End Damage Patterns (Sample: ERR5729614)", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()


import os
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pymc as pm
import arviz as az

def run_damage_model(sample_name, positions, success_counts, fail_counts, label, color):
    N_total = success_counts + fail_counts

    with pm.Model() as model:
        alpha = pm.Beta('alpha', alpha=1, beta=1)
        beta = pm.HalfNormal('beta', sigma=1)
        p = alpha * pm.math.exp(-beta * (positions - 1))
        obs = pm.Binomial('obs', n=N_total, p=p, observed=success_counts)

        trace = pm.sample(
            draws=500,
            tune=250,
            target_accept=0.9,
            cores=1,
            progressbar=False,
            return_inferencedata=True
        )

    posterior_samples = trace.posterior.stack(samples=("chain", "draw"))
    alpha_samples = posterior_samples["alpha"].values
    beta_samples = posterior_samples["beta"].values

    pos_range = np.arange(2, 26)
    mean_p, hdi_lower, hdi_upper = [], [], []

    for i in pos_range:
        probs = alpha_samples * np.exp(-beta_samples * (i - 1))
        mean_p.append(np.mean(probs))
        hdi = az.hdi(probs, hdi_prob=0.95)
        hdi_lower.append(hdi[0])
        hdi_upper.append(hdi[1])

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(pos_range, mean_p, color=color, label=f"{label} estimate")
    ax.fill_between(pos_range, hdi_lower, hdi_upper, color=color, alpha=0.3, label="95% Credible Interval")
    ax.set_xlabel("Position from read end")
    ax.set_ylabel("Damage rate")
    ax.set_title(f"{label} decay model – {sample_name}")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    return fig

# Pair 5' and 3' end files
files_5 = sorted(glob("/mnt/c/Users/slwrigh4/ASU Dropbox/Sterling Wright/PC/Desktop/PSU/Romania/03_UPDATED_ANALYSES_COMPLETE_DATASET/07_ChangePoint/freq/*_5_end_freq"))
files_3 = sorted(glob("/mnt/c/Users/slwrigh4/ASU Dropbox/Sterling Wright/PC/Desktop/PSU/Romania/03_UPDATED_ANALYSES_COMPLETE_DATASET/07_ChangePoint/freq/*_3_end_freq"))

paired_samples = [
    (os.path.basename(f5).replace("_5_end_freq", ""), f5, f3)
    for f5, f3 in zip(files_5, files_3)
    if os.path.basename(f5).replace("_5_end_freq", "") == os.path.basename(f3).replace("_3_end_freq", "")
]

# Output PDF
output_pdf = "/mnt/c/Users/slwrigh4/ASU Dropbox/Sterling Wright/PC/Desktop/PSU/Romania/03_UPDATED_ANALYSES_COMPLETE_DATASET/07_ChangePoint/freq/mcmc_damage_summary.pdf"
with PdfPages(output_pdf) as pdf:
    for sample, file_5, file_3 in paired_samples:
        try:
            df_5 = pd.read_csv(file_5, sep="\t")
            df_3 = pd.read_csv(file_3, sep="\t")

            # 5' end
            df5sub = df_5[df_5['Position_from_5end'].between(2, 10)]
            pos5 = df5sub['Position_from_5end'].values
            T = df5sub['T_freq'].values
            C = df5sub['C_freq'].values
            fig_ct = run_damage_model(sample, pos5, T, C, label="C→T (5′)", color="firebrick")
            pdf.savefig(fig_ct)
            plt.close(fig_ct)

            # 3' end
            df3sub = df_3[df_3['Position_from_3end'].between(2, 10)]
            pos3 = df3sub['Position_from_3end'].values
            A = df3sub['A_freq'].values
            G = df3sub['G_freq'].values
            fig_ga = run_damage_model(sample, pos3, A, G, label="G→A (3′)", color="darkblue")
            pdf.savefig(fig_ga)
            plt.close(fig_ga)

        except Exception as e:
            print(f"Error processing {sample}: {e}")


import os
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pymc as pm
import arviz as az

# Updated plotting function with optional axis flipping and cleaned sample name
def run_damage_model(sample_name, positions, success_counts, fail_counts, label, color, flip_x=False):
    N_total = success_counts + fail_counts

    with pm.Model() as model:
        alpha = pm.Beta('alpha', alpha=1, beta=1)
        beta = pm.HalfNormal('beta', sigma=1)
        p = alpha * pm.math.exp(-beta * (positions - 1))
        obs = pm.Binomial('obs', n=N_total, p=p, observed=success_counts)

        trace = pm.sample(
            draws=500,
            tune=250,
            target_accept=0.9,
            cores=1,
            progressbar=False,
            return_inferencedata=True
        )

    posterior_samples = trace.posterior.stack(samples=("chain", "draw"))
    alpha_samples = posterior_samples["alpha"].values
    beta_samples = posterior_samples["beta"].values

    pos_range = np.arange(2, 26)
    mean_p, hdi_lower, hdi_upper = [], [], []

    for i in pos_range:
        probs = alpha_samples * np.exp(-beta_samples * (i - 1))
        mean_p.append(np.mean(probs))
        hdi = az.hdi(probs, hdi_prob=0.95)
        hdi_lower.append(hdi[0])
        hdi_upper.append(hdi[1])

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(pos_range, mean_p, color=color, label=f"{label} estimate")
    ax.fill_between(pos_range, hdi_lower, hdi_upper, color=color, alpha=0.3, label="95% Credible Interval")
    ax.set_xlabel("Position from read end")
    ax.set_ylabel("Damage rate")
    
    # Clean up sample name
    sample_clean = sample_name.replace("filtered_25bp_100k_", "")
    ax.set_title(f"{label} decay model – {sample_clean}")

    if flip_x:
        ax.set_xlim(ax.get_xlim()[::-1])  # reverse the x-axis

    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    return fig

"Function redefined successfully with required modifications."


# Convert the script into a Jupyter Notebook format
import nbformat as nbf

# Define notebook content
nb = nbf.v4.new_notebook()
cells = []

# Markdown cell
cells.append(nbf.v4.new_markdown_cell("# DamageBayes: MCMC-Based Ancient DNA Damage Assessment\nThis notebook runs PyMC-based Bayesian models to estimate postmortem DNA damage patterns across samples."))

# Imports
cells.append(nbf.v4.new_code_cell("""\
import os
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pymc as pm
import arviz as az
"""))

# Parameters and setup
cells.append(nbf.v4.new_code_cell("""\
# Set input/output directories
input_dir = "/mnt/c/Users/slwrigh4/ASU Dropbox/Sterling Wright/PC/Desktop/PSU/Romania/03_UPDATED_ANALYSES_COMPLETE_DATASET/07_ChangePoint/freq"
output_pdf = os.path.join(input_dir, "mcmc_damage_summary_combined.pdf")
output_csv = os.path.join(input_dir, "global-damage-results.csv")
"""))

# Define modeling function
cells.append(nbf.v4.new_code_cell("""\
def run_damage_model(sample_name, positions, success_counts, fail_counts, label, color, flip_x=False):
    N_total = success_counts + fail_counts

    with pm.Model() as model:
        alpha = pm.Beta('alpha', alpha=1, beta=1)
        beta = pm.HalfNormal('beta', sigma=1)
        p = alpha * pm.math.exp(-beta * (positions - 1))
        obs = pm.Binomial('obs', n=N_total, p=p, observed=success_counts)

        trace = pm.sample(
            draws=500,
            tune=250,
            target_accept=0.9,
            cores=1,
            progressbar=False,
            return_inferencedata=True
        )

    # Posterior summaries
    summary = az.summary(trace, var_names=["alpha", "beta"], round_to=3)
    alpha_mean = summary.loc["alpha", "mean"]
    beta_mean = summary.loc["beta", "mean"]

    # Posterior predictions
    posterior_samples = trace.posterior.stack(samples=("chain", "draw"))
    alpha_samples = posterior_samples["alpha"].values
    beta_samples = posterior_samples["beta"].values

    pos_range = np.arange(2, 26)
    mean_p, hdi_lower, hdi_upper = [], [], []

    for i in pos_range:
        probs = alpha_samples * np.exp(-beta_samples * (i - 1))
        mean_p.append(np.mean(probs))
        hdi = az.hdi(probs, hdi_prob=0.95)
        hdi_lower.append(hdi[0])
        hdi_upper.append(hdi[1])

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(pos_range, mean_p, color=color, label=f"{label} estimate")
    ax.fill_between(pos_range, hdi_lower, hdi_upper, color=color, alpha=0.3, label="95% Credible Interval")
    ax.set_xlabel("Position from read end")
    ax.set_ylabel("Damage rate")
    ax.set_title(f"{label} – {sample_name.replace('filtered_25bp_100k_', '')}")

    if flip_x:
        ax.set_xlim(ax.get_xlim()[::-1])

    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    
    return fig, alpha_mean, beta_mean
"""))

# Run modeling across samples
cells.append(nbf.v4.new_code_cell("""\
# Gather files
files_5 = sorted(glob(os.path.join(input_dir, "*_5_end_freq")))
files_3 = sorted(glob(os.path.join(input_dir, "*_3_end_freq")))

paired_samples = [
    (os.path.basename(f5).replace("_5_end_freq", ""), f5, f3)
    for f5, f3 in zip(files_5, files_3)
    if os.path.basename(f5).replace("_5_end_freq", "") == os.path.basename(f3).replace("_3_end_freq", "")
]

results = []

# PDF output
with PdfPages(output_pdf) as pdf:
    for sample, file_5, file_3 in paired_samples:
        try:
            df_5 = pd.read_csv(file_5, sep="\\t")
            df_3 = pd.read_csv(file_3, sep="\\t")

            df5sub = df_5[df_5['Position_from_5end'].between(2, 10)]
            pos5 = df5sub['Position_from_5end'].values
            T = df5sub['T_freq'].values
            C = df5sub['C_freq'].values
            fig_ct, alpha_ct, beta_ct = run_damage_model(sample, pos5, T, C, label="C→T (5′)", color="firebrick")

            df3sub = df_3[df_3['Position_from_3end'].between(2, 10)]
            pos3 = df3sub['Position_from_3end'].values
            A = df3sub['A_freq'].values
            G = df3sub['G_freq'].values
            fig_ga, alpha_ga, beta_ga = run_damage_model(sample, pos3, A, G, label="G→A (3′)", color="darkblue", flip_x=True)

            # Combine and save
            fig_combined, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
            for ax, fig in zip([ax1, ax2], [fig_ct, fig_ga]):
                for child in fig.get_children():
                    if isinstance(child, plt.Axes):
                        for artist in child.get_children():
                            ax.add_artist(artist)
            fig_combined.suptitle(sample.replace("filtered_25bp_100k_", ""), fontsize=14, fontweight='bold')
            pdf.savefig(fig_combined)
            plt.close('all')

            # Append summary
            results.append({
                "Sample": sample.replace("filtered_25bp_100k_", ""),
                "Alpha_CtoT_Mean": round(alpha_ct, 3),
                "Beta_CtoT_Mean": round(beta_ct, 3),
                "Alpha_GtoA_Mean": round(alpha_ga, 3),
                "Beta_GtoA_Mean": round(beta_ga, 3)
            })
        except Exception as e:
            print(f"Error processing {sample}: {e}")
"""))

# Save results to CSV
cells.append(nbf.v4.new_code_cell("""\
df_results = pd.DataFrame(results)
df_results.to_csv(output_csv, index=False)
df_results.head()
"""))

# Save notebook
nb['cells'] = cells
notebook_path = "/mnt/c/Users/slwrigh4//Desktop/PSU/Romania/03_UPDATED_ANALYSES_COMPLETE_DATASET/07_ChangePoint/Damage-Bayes_Batch_Analysis.ipynb"
with open(notebook_path, 'w') as f:
    nbf.write(nb, f)

notebook_path
