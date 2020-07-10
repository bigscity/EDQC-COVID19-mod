# **EDQC-COVID19-mod**

A  modified SEIR model which introduced dynamic parameters to incoporate the effects of targeted interventions (TIs) for COVID-19.

## Requirements

All required packages are list in `environment.yml` which is available from [Anaconda](https://anaconda.org/). Install with command on CMD or terminal:

```
conda env create -f environment.yml
```

## Simulating scenarios with or without some of Tis

Run `model_simulation.py` after all requirements installed. Both figures and CSV data would be generated in `./simuresult/` within minutes.

- Figures are produced in PDF format, with Figure 4A, 4B each a single file and Figure 4C-4H in one single file `intervention_simulation.pdf` with multiple plots.
- The detailed results are listed as CSV files, which are named by the format `resultsimuXYZ.csv`. `X, Y, Z` are dummy variables indicating whether a targeted intervention is implemented (1) or not (0).
  - `X` represents the implementation of localized lockdown.
  - `Y` represents the implementation of close-contact tracing.
  - `Z` represents the implementation of community-based nucleic acid testing (NAT) .
  - For example, the result file for scenario without any TIs is `resultsimu000.csv`; the result file for scenario with only intervention of close-contact tracing is `resultsimu010.csv`.
- For each CSV file, the following columns are contained
  - `Date`  specifies the date of rows.
  - `reported_I` actually reported number of cumulative cases by date of onset.
  - `conbined_intervention_I`  estimated number of cumulative cases with all TIs implemented, *i.e.*, the actual scenario in Beijing.
  - `simulation_I` estimated number of cumulative cases under simulation scenario with or without some of TIs, which is defined by the name of CSV.
  - Both `conbined_intervention_I` and `simulation_I` are median value of the simulation, the IQR is also provided in `conbined_intervention_I_up(/down)` and  `simulation_I(/down)`.

## Structure of this repository 

- `./modeldata/` contains necessary data for model.
- `./objs/` contains parameters of the model.
- `./simuresult/` contains output of simulation.
- `./model_simulation.py` the main code of simulation.
- `./utils.py` contains som necessary funtion of model.
- `./Pearson_corr_check.R` code for computing the Pearson's correlation coefficients between estimated cumulative numbers and reported cumulative numbers. (should run in `R v3.6.0`)
- `./environment.yml` lists the packages of the environment.

## License

This repository is licensed under the GNU General Public License (GPL).