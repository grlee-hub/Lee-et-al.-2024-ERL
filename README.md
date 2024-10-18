# Statistical seasonal prediction of Arctic sea ice concentration based on spatiotemporal anomaly persistent method

## Overview

This repository contains the code for the S-EOF (Season-reliant Empirical Orthogonal Function) model, a statistical approach for predicting Arctic sea ice concentration. The model leverages the slow-changing properties of sea ice, incorporating a form of "memory" into its predictions.

## Repository Structure

### /data

Input data for the model consists of two main sources:

1. HadISST sea ice concentration (sic) data
2. ERA5 mean sea level pressure (msl) data

These datasets can be downloaded from their respective official sources. Please ensure you have the necessary permissions to use these datasets.

### /function

The model utilizes the persistent nature of sea ice by calculating sea ice area from sea ice concentration. This calculation is performed using a custom function stored in `siecalcTLL.ncl`. This function is called within `model.ncl` to compute the sea ice area.

### /output and /output_msl

The model can be run with or without atmospheric correction, which attempts to simulate rapid changes in sea ice due to atmospheric circulation. Results are stored in two separate directories:

- `/output`: Contains results without atmospheric correction
- `/output_msl`: Contains results with atmospheric correction

The `optMSLP` option in `model.ncl` controls whether atmospheric correction is applied.

## Example Output

For demonstration purposes, we have included example output for a model initialized in December 2021. You can find these results in the appropriate output directory based on whether atmospheric correction was applied.

## Usage

To use this model:

1. Ensure you have the required input data in the `/data` directory.
2. Adjust parameters in `model.ncl` as needed, including the `optMSLP` option for atmospheric correction.
3. Run `model.ncl` to generate predictions.
4. Check the appropriate output directory for results.

## Citation

If you use this model in your research, please cite our paper:

```
Lee, G. R., Woo, S. H., Baek, E. H., Kim, J. H., Kim, B. M., & Jeong, J. H. (2024). Statistical seasonal prediction of Arctic sea ice concentration based on spatiotemporal anomaly persistent method. Environmental Research Letters.
```

## Contact

For questions or issues regarding this model, please conatact to jjeehoon@gmail.com or grlee0313@gmail.com
