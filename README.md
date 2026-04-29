# ElasticityPdPvDerModel

This repository contains a two-stage distribution-system planning workflow:

1. Python scripts preprocess load and PV time series into monthly scenarios.
2. AMPL models solve the PV allocation problem and then the distribution planning problem with the PV results as input.

The code is organized around a single feeder test case and produces scenario CSV files, solver outputs, and figures.

## Overview

The main flow is:

1. Run `data_pross.py` to cluster the raw load and PV traces into monthly scenario curves and probabilities.
2. Use `pv_decision.mod` in AMPL to solve the PV allocation problem and generate `pv-solution.csv`.
3. Run `der_decision.run` or `MSD_der_decision.run` in AMPL to solve the network planning model using the PV allocation results.
4. Run `plots.py` and `analysis.py` to generate figures from the exported results.

## Requirements

- Python 3.10+.
- Python packages: `pandas`, `numpy`, `matplotlib`, and `scikit-learn`.
- AMPL.
- A MILP/NLP solver compatible with AMPL, configured in the run files as `cplex`.
- The AMPL CSV interface library referenced by the models as `amplcsv.dll`.

## Repository Layout

- `data_pross.py` prepares the monthly load and PV scenario CSVs used by the AMPL models.
- `pv_decision.mod` solves the PV investment/allocation problem.
- `model.mod` and `MSD_model.mod` define the distribution-system planning models.
- `der_decision.run` runs the single-month distribution planning workflow.
- `MSD_der_decision.run` runs the monthly/scenario version of the planning workflow.
- `plots.py` builds the scenario plot figure in `figures/`.
- `analysis.py` post-processes solver outputs into publication-style figures.
- `input-data/` stores raw inputs and preprocessing outputs.
- `results/` stores solver outputs.
- `figures/` stores generated plots.

## Quick Start

### 1. Prepare the Python environment

Install the Python dependencies and confirm that the input files are present in `input-data/`:

- `load_clean.csv`
- `ninja_pv_52.1305_4.9992_uncorrected.csv`
- `nodes.csv`
- `profiles.csv`

### 2. Generate the scenario inputs

Run `data_pross.py` to build these files in `input-data/`:

- `parametros_load_curvas.csv`
- `parametros_load_probs.csv`
- `parametros_pv_curvas.csv`
- `parametros_pv_probs.csv`

These files contain the monthly scenario centroids and their probabilities for the AMPL models.

### 3. Solve the PV allocation problem

Open `pv_decision.mod` in AMPL and solve it with the required `input-data/` files loaded. The distribution planning scripts expect a PV allocation file named `pv-solution.csv` with the columns:

- `N`
- `PV_allocation`
- `PV_inverter`
- `PV_panels`

If you already have numbered PV outputs in `results/` such as `pv-solution_1.csv`, rename or copy the desired file to `pv-solution.csv` before running the next step.

### 4. Solve the distribution planning model

Run one of the AMPL scripts:

- `der_decision.run` for the base planning model.
- `MSD_der_decision.run` for the monthly/scenario planning variant.

Both scripts read the PV allocation file and solve the network planning model with `cplex`.

### 5. Generate plots and analysis outputs

- Run `plots.py` to generate `figures/curvas_scenarios.svg`.
- Run `analysis.py` after the AMPL solve to generate substation power plots in `figures/`.

## Inputs and Outputs

### Inputs

- `input-data/load_clean.csv` and `input-data/ninja_pv_52.1305_4.9992_uncorrected.csv` are the raw traces used by `data_pross.py`.
- `input-data/nodes.csv` and `input-data/profiles.csv` provide the network and profile data used by the AMPL models.
- `input-data/parametros_*.csv` are generated scenario files consumed by `MSD_der_decision.run` and the plotting scripts.

### Outputs

- `results/se_caseA.csv` is the main substation output used by `analysis.py`.
- `results/data_voltage.csv` is the voltage output written by the planning models when the corresponding table export is enabled.
- `figures/curvas_scenarios.svg` is the scenario summary plot from `plots.py`.
- `figures/se_power_delivery_*.svg` are the substation delivery plots from `analysis.py`.

## Notes

- The AMPL run files currently reference `pv-solution.csv`, but the repository contains numbered files such as `results/pv-solution_1.csv` through `results/pv-solution_20.csv`. Make sure the file name matches before solving.
- `analysis.py` also references `Voltage_WithoutAllocation_.csv`, which is not present in the repository. Either provide that file or update the script before running the voltage section.
- The model files are configured for `cplex`; if you use another solver, update the `option solver` lines in the AMPL run files.

## Typical Workflow

1. Prepare the Python environment.
2. Run `data_pross.py`.
3. Solve `pv_decision.mod` in AMPL.
4. Run `der_decision.run` or `MSD_der_decision.run` in AMPL.
5. Run `plots.py` and `analysis.py` to generate figures.