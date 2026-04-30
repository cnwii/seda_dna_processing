#!/usr/bin/env python
import os
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pymc as pm
import arviz as az


input_dir = "/hpcfs/users/a1226686/PROJECTS/Bushfire/3_SEQUENCES/deduplicated/damabayes"

output_pdf = os.path.join(input_dir, "damagebayes_with_smiley_plots.pdf")
output_csv = os.path.join(input_dir, "global-damage-results.csv")


def run_damage_model(positions, success_counts, fail_counts):
    N_total = success_counts + fail_counts

    with pm.Model() as model:
        alpha = pm.Beta("alpha", alpha=1, beta=1)
        beta = pm.HalfNormal("beta", sigma=1)

        p = alpha * pm.math.exp(-beta * (positions - 1))

        pm.Binomial(
            "obs",
            n=N_total,
            p=p,
            observed=success_counts
        )

        trace = pm.sample(
            draws=1000,
            tune=1000,
            target_accept=0.95,
            chains=4,
            cores=1,
            progressbar=False,
            return_inferencedata=True
        )

    summary = az.summary(trace, var_names=["alpha", "beta"])

    alpha_mean = summary.loc["alpha", "mean"]
    beta_mean = summary.loc["beta", "mean"]

    posterior_samples = trace.posterior.stack(samples=("chain", "draw"))
    alpha_samples = posterior_samples["alpha"].values
    beta_samples = posterior_samples["beta"].values

    pos_range = np.arange(2, 26)

    mean_p = []
    hdi_lower = []
    hdi_upper = []

    for i in pos_range:
        probs = alpha_samples * np.exp(-beta_samples * (i - 1))
        mean_p.append(np.mean(probs))

        hdi = az.hdi(probs, hdi_prob=0.95)
        hdi_lower.append(hdi[0])
        hdi_upper.append(hdi[1])

    return pos_range, mean_p, hdi_lower, hdi_upper, alpha_mean, beta_mean


files_5 = sorted(glob(os.path.join(input_dir, "*_5_end_freq")))
files_3 = sorted(glob(os.path.join(input_dir, "*_3_end_freq")))

paired_samples = []

for f5 in files_5:
    sample = os.path.basename(f5).replace("_5_end_freq", "")
    f3 = os.path.join(input_dir, sample + "_3_end_freq")

    if f3 in files_3:
        paired_samples.append((sample, f5, f3))
    else:
        print(f"Warning: no matching 3′ file found for {sample}")


results = []

with PdfPages(output_pdf) as pdf:

    for sample, file_5, file_3 in paired_samples:

        alpha_ct = beta_ct = alpha_ga = beta_ga = None

        try:
            print(f"Processing {sample}...")

            df_5 = pd.read_csv(file_5, sep="\t")
            df_3 = pd.read_csv(file_3, sep="\t")

            sample_clean = sample.replace("filtered_25bp_100k_", "")

            # Raw smiley plot data: positions 2–25
            df_5_plot = df_5[df_5["Position_from_5end"].between(2, 25)].copy()
            df_3_plot = df_3[df_3["Position_from_3end"].between(2, 25)].copy()

            # DamageBayes model data: positions 2–10
            df5sub = df_5[df_5["Position_from_5end"].between(2, 10)]
            df3sub = df_3[df_3["Position_from_3end"].between(2, 10)]

            pos5 = df5sub["Position_from_5end"].values
            T = df5sub["T_freq"].values
            C = df5sub["C_freq"].values

            pos3 = df3sub["Position_from_3end"].values
            A = df3sub["A_freq"].values
            G = df3sub["G_freq"].values

            pos5_model, mean_ct, hdi_ct_lower, hdi_ct_upper, alpha_ct, beta_ct = run_damage_model(
                pos5,
                T,
                C
            )

            pos3_model, mean_ga, hdi_ga_lower, hdi_ga_upper, alpha_ga, beta_ga = run_damage_model(
                pos3,
                A,
                G
            )

            fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey="row")

            ax1, ax2 = axes[0]
            ax3, ax4 = axes[1]

            # -------------------------
            # Row 1: raw smiley plots
            # -------------------------

            ax1.plot(
                df_5_plot["Position_from_5end"],
                df_5_plot["T_freq"],
                marker="o",
                color="firebrick"
            )
            ax1.set_title("Raw 5′ T frequency")
            ax1.set_xlabel("Position from 5′ end")
            ax1.set_ylabel("Raw frequency")
            ax1.grid(True)

            ax2.plot(
                df_3_plot["Position_from_3end"],
                df_3_plot["A_freq"],
                marker="o",
                color="darkblue"
            )
            ax2.set_title("Raw 3′ A frequency")
            ax2.set_xlabel("Position from 3′ end")
            ax2.invert_xaxis()
            ax2.grid(True)

            # -------------------------
            # Row 2: DamageBayes model
            # -------------------------

            ax3.plot(
                pos5_model,
                mean_ct,
                marker="o",
                color="firebrick",
                label="C→T 5′ model"
            )
            ax3.fill_between(
                pos5_model,
                hdi_ct_lower,
                hdi_ct_upper,
                color="firebrick",
                alpha=0.3
            )
            ax3.set_title("DamageBayes 5′ C→T decay")
            ax3.set_xlabel("Position from 5′ end")
            ax3.set_ylabel("Modelled damage rate")
            ax3.grid(True)
            ax3.legend()

            ax4.plot(
                pos3_model,
                mean_ga,
                marker="o",
                color="darkblue",
                label="G→A 3′ model"
            )
            ax4.fill_between(
                pos3_model,
                hdi_ga_lower,
                hdi_ga_upper,
                color="darkblue",
                alpha=0.3
            )
            ax4.set_title("DamageBayes 3′ G→A decay")
            ax4.set_xlabel("Position from 3′ end")
            ax4.invert_xaxis()
            ax4.grid(True)
            ax4.legend()

            fig.suptitle(sample_clean, fontsize=16, fontweight="bold")
            plt.tight_layout(rect=[0, 0, 1, 0.96])

            pdf.savefig(fig)
            plt.close(fig)

        except Exception as e:
            print(f"Error processing {sample}: {e}")

        finally:
            results.append({
                "Sample": sample.replace("filtered_25bp_100k_", ""),
                "Alpha_CtoT_Mean": round(alpha_ct, 3) if alpha_ct is not None else None,
                "Beta_CtoT_Mean": round(beta_ct, 3) if beta_ct is not None else None,
                "Alpha_GtoA_Mean": round(alpha_ga, 3) if alpha_ga is not None else None,
                "Beta_GtoA_Mean": round(beta_ga, 3) if beta_ga is not None else None
            })


df_results = pd.DataFrame(results)
df_results.to_csv(output_csv, index=False)

print(f"Processed {len(paired_samples)} paired samples")
print(f"PDF saved to: {output_pdf}")
print(f"CSV saved to: {output_csv}")
