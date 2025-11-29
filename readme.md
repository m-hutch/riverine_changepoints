# Riverine Changepoints

This repository reviews methods of identifying critical changepoints in USGS river gage data.

## Files

`data_explor.R` - an introduction to the example data set: river gages in Kerr County, Texas and their discharge values in the summer of 2025.

`utils.R` - utility functions used for pulling and cleaning data, and plotting.

`analytical_bocpd` - directory of files related to an analytical Bayesian online changepoint detection (BOCPD) solution using a gamma-gamma conjugate prior.

`numerical_bocpd` - directory of files related to a numerical BOCPD solution using `cmdrstan` and informed priors.

`ensemble_algorithm_bocpd` - directory of files related to an example using `Rbeast` package using uninformed priors.

`presentation` - directory with presentation slides and notes
