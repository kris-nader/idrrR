
# idrrR: Integrated Drug Response Resource Tool

`idrrR` is an R package designed to harmonize and analyze pharmacogenomic drug response data. It serves as the computational engine behind the **iDRR (Integrated Drug Response Resource)** project, enabling researchers to calculate robust consensus drug response curves and compare their own experimental data against a large-scale integrated database.

## About the Project

Pharmacogenomic datasets (such as GDSC, CTRP, and others) often show discordant results for the same drug-cell line pairs due to variations in experimental protocols, assay types, and analytical methods. **iDRR** addresses this inconsistency by integrating multiple datasets to derive a **"Consensus Curve"**-a unified, robust estimate of drug sensitivity that mitigates noise and outliers.

**Key Features:**
*   **Consensus Modeling:** fits a 4-parameter logistic (4-PL) curve to integrated data from multiple studies.
*   **omprehensive Metrics:** Implements the novel Drug Response Score as well as Standard Error of Estimate (SEE).
*   **Benchmarking:** Allows users to overlay their own data on consensus curves to quantify deviation (Standard Error of Estimate).


---

## Installation

You can install `idrrR` directly from the source. Ensure you have the required dependencies:

```r
# Install idrrR
devtools::install_github("kris-nader/idrrR")
```

---

## Quick Start Guide

### 1. Load the Library
Loading the library automatically initializes the consensus database (~70MB) for fast lookups.

```r
library(idrrR)
```

### 2. View Consensus Curves
Query the database for any Drug-Cell Line pair to see the consensus fit.

```r
# Plot Paclitaxel response in MCF7 cells
result <- plot_consensus_from_db(
  drug_name = "Paclitaxel", 
  cell_line = "MCF7"
)

# Render the plot
print(result$plot)

# Inspect the calculated metrics
print(result$metrics)
```

### 3. Analyze Your Own Data
You can analyze your own screening data using the `analyze_own_data` function.

**Input Requirements:**
*   **Columns:** `Drug_Name`, `Cell_Line`, `Dose` (µM), `Inhibition` (%).
*   **Scale:** `Inhibition` should be 0-100 (or 0-1 fraction, auto-detected).
    *   *Note: Convert Viability to Inhibition (100 - Viability) before analysis.*

```r
# Create sample data
my_data <- data.frame(
  Drug_Name = rep("TestDrug", 5),
  Cell_Line = rep("TestCell", 5),
  Dose = c(0.001, 0.01, 0.1, 1, 10),
  Inhibition = c(5, 15, 50, 85, 95)
)

# Analyze and calculate metrics
results <- analyze_own_data(my_data)

# View results (returns a list with $metrics and $plots)
print(results$metrics)
print(results$plots[[1]])
```

### 4. Benchmark Against Consensus
Compare your experimental results with the established consensus to validate your findings.

```r
# Overlay your points on the global consensus
comparison <- plot_custom_on_consensus(
  drug_name = "Paclitaxel",
  cell_line = "MCF7",
  dose = c(0.001, 0.01, 0.1, 1, 10),  # Your doses
  inhibition = c(0, 10, 45, 80, 95),  # Your inhibition values
  dataset_name = "My_Lab_Data"
)

# Visual comparison
print(comparison$plot)

# Quantify deviation (Standard Error of Estimate)
print(comparison$SEE_custom_vs_consensus)
```

### 5. Quantifying Deviation (SEE)
When benchmarking, the tool calculates the **Standard Error of Estimate (SEE)** to quantify how much your data deviates from the consensus model.

**Formula:**
$$ SEE = \sqrt{ \frac{\sum (y_{custom} - y_{consensus})^2}{n - p} } $$
Where:
*   $y_{custom}$: Your observed inhibition values.
*   $y_{consensus}$: Inhibition predicted by the consensus curve at your doses.
*   $n$: Number of data points.
*   $p$: Number of model parameters.

**Interpretation:**
*   **Low SEE (< 10):** High consistency. Your data closely matches the aggregated consensus.
*   **High SEE (> 35):** Significant deviation. This may indicate experimental differences (assay type, passage number) or true biological variation.

#### Example: Validating Consistency
To demonstrate, we can generate "perfect" data points from the consensus curve and then introduce noise to see how SEE changes.

```r
# 1. Get Consensus Parameters for Paclitaxel x MCF7
consensus <- plot_consensus_from_db("Paclitaxel", "MCF7")
params <- consensus$metrics
slope <- params$SLOPE
ic50 <- params$IC50ABS
max_resp <- params$Max_Response

# 2. Generate "Perfect" Data Points (lying exactly on the curve)
doses <- c(0.001, 0.01, 0.1, 1, 10)
# Formula: y = 0 + (max - 0) / (1 + 10^(slope * (log10(ic50) - log10(dose))))
perfect_inhib <- max_resp / (1 + 10^(slope * (log10(ic50) - log10(doses))))

# Benchmark Perfect Data
res_perfect <- plot_custom_on_consensus("Paclitaxel", "MCF7", doses, perfect_inhib, "Perfect_Fit")
print(paste("SEE (Perfect):", round(res_perfect$SEE_custom_vs_consensus, 4)))
# Output: SEE (Perfect): 0

# 3. Benchmark "Skewed" Data (Simulating experimental error)
skewed_inhib <- perfect_inhib
skewed_inhib[3] <- skewed_inhib[3] + 25 # Add 25% deviation to one point

res_skewed <- plot_custom_on_consensus("Paclitaxel", "MCF7", doses, skewed_inhib, "Skewed_Fit")
print(paste("SEE (Skewed):", round(res_skewed$SEE_custom_vs_consensus, 4)))
# Output: SEE (Skewed): ~17.5 (Significant deviation)
```
