# This script demonstrates the PERCEPT method for determining protein abundance thresholds
# It first shows a single iteration with visualization to explain the method
# Then performs multiple iterations to determine robust threshold values

### [1] SET UP WORKSPACE ###
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
import matplotlib.pyplot as plt
from pathlib import Path

# Set project root directory and paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
figures_dir = project_root / "data" / "figures" / "percept"
results_dir = project_root / "data" / "processed" / "results"

# Load data and define constants
df = pd.read_excel(input_dir / "4_pg_chivas_filtered.xlsx")
control_timepoints = ['screen', 'prior', '1mth', '3mth']
patients = [2, 9, 10, 13, 17, 20, 25, 26, 30, 32, 33, 43, 55, 57, 59, 61, 64, 68, 71, 75]

### [2] DEFINE FUNCTIONS ###
def percept(m0, m1, F, p):
    return m0 + ((m0 - m1) / -(F**p))

def apply_percept(data, hypothetical_mean, penalty):
    tval, pval = ttest_1samp(data, popmean=hypothetical_mean, nan_policy='omit')
    sample_mean = np.mean(data)
    scaled_value = percept(m0=hypothetical_mean, m1=sample_mean, F=penalty, p=pval)
    return scaled_value, pval

def perform_analysis():
    # Calculate ratios for all patients
    ratios_df = pd.DataFrame()
    for patient in patients:
        control_samples = [col for col in df.columns if col.startswith(f"{patient}_") and col.split('_')[1] in control_timepoints]
        if len(control_samples) >= 2:
            sample1, sample2 = np.random.choice(control_samples, 2, replace=False)
            ratio = df[sample1] / df[sample2]
            ratios_df[f'{patient}_ratio'] = ratio

    # Convert to log2-fold change
    log2_fc_df = np.log2(ratios_df)
    numeric_cols = log2_fc_df.columns

    # Apply PERCEPT scaling
    results = log2_fc_df[numeric_cols].apply(lambda row: apply_percept(row.dropna(), m0, F), axis=1)
    log2_fc_df['PERCEPT_scaled'] = results.apply(lambda x: x[0])

    # Calculate thresholds
    thresholds = {
        '95': (np.percentile(log2_fc_df['PERCEPT_scaled'].dropna(), 2.5),
               np.percentile(log2_fc_df['PERCEPT_scaled'].dropna(), 97.5)),
        '98': (np.percentile(log2_fc_df['PERCEPT_scaled'].dropna(), 1),
               np.percentile(log2_fc_df['PERCEPT_scaled'].dropna(), 99)),
        '99': (np.percentile(log2_fc_df['PERCEPT_scaled'].dropna(), 0.5),
               np.percentile(log2_fc_df['PERCEPT_scaled'].dropna(), 99.5))
    }
    
    return log2_fc_df, thresholds

### [3] DEMONSTRATE METHOD WITH TWO ITERATIONS ###
# Set parameters
m0 = 0  # hypothetical mean
F = 10 * len(patients)  # penalty factor

def plot_percept_distribution(iteration_number):
    # Perform analysis
    demo_results, demo_thresholds = perform_analysis()
    
    # Create plot
    plt.figure(figsize=(12, 6))
    plt.hist(demo_results['PERCEPT_scaled'].dropna(), bins=50, edgecolor='black')
    
    # Add threshold lines and labels
    colors = {'95': 'r', '98': 'y', '99': 'g'}
    for ci, (lower, upper) in demo_thresholds.items():
        # Add lines and labels for lower threshold
        plt.axvline(lower, color=colors[ci], linestyle='dashed', linewidth=2, label=f'{ci}% CI')
        plt.text(lower, plt.gca().get_ylim()[1], f'{lower:.3f}', 
                va='bottom', ha='right', color=colors[ci])
        
        # Add lines and labels for upper threshold
        plt.axvline(upper, color=colors[ci], linestyle='dashed', linewidth=2)
        plt.text(upper, plt.gca().get_ylim()[1], f'{upper:.3f}', 
                va='bottom', ha='right', color=colors[ci])

    plt.xlim(-1.5, 1.5)
    plt.title(f'Distribution of PERCEPT-scaled Control Sample Ratios\nIteration {iteration_number}')
    plt.xlabel('PERCEPT-scaled Ratio')
    plt.ylabel('Frequency')
    plt.legend()
    plt.tight_layout()
    
    # Save and show plot
    plt.savefig(figures_dir / f"percept_demonstration_iter{iteration_number}.svg", bbox_inches='tight')
    plt.show()

# Run two iterations
plot_percept_distribution(1)
plot_percept_distribution(2)

### [4] PERFORM MULTIPLE ITERATIONS ###
num_iterations = 100
all_thresholds = []

for i in range(num_iterations):
    _, thresholds = perform_analysis()
    all_thresholds.append([
        thresholds['95'][0], thresholds['95'][1],
        thresholds['98'][0], thresholds['98'][1],
        thresholds['99'][0], thresholds['99'][1]
    ])

# Create results DataFrame
results_df = pd.DataFrame(all_thresholds, 
                         columns=['95% CI Lower', '95% CI Upper',
                                 '98% CI Lower', '98% CI Upper',
                                 '99% CI Lower', '99% CI Upper'])

# Calculate summary statistics
summary = results_df.agg(['mean', 'std', 'min', 'max'])

# Save results
results_df.to_excel(results_dir / 'percept_threshold_iterations.xlsx', index=False)
summary.to_excel(results_dir / 'percept_threshold_summary.xlsx')

print("\nFinal threshold values (mean ± std):")
for ci in ['95%', '98%', '99%']:
    lower_mean = summary.loc['mean', f'{ci} CI Lower']
    lower_std = summary.loc['std', f'{ci} CI Lower']
    upper_mean = summary.loc['mean', f'{ci} CI Upper']
    upper_std = summary.loc['std', f'{ci} CI Upper']
    print(f"\n{ci} CI:")
    print(f"Lower: {lower_mean:.4f} ± {lower_std:.4f}")
    print(f"Upper: {upper_mean:.4f} ± {upper_std:.4f}")