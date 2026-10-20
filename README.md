# Weibull CUSUM Control Charts under Type-I Progressive Censoring

A simulation and real-data study of likelihood-ratio CUSUM control charts for monitoring Weibull lifetimes when test data is progressively censored. The project covers maximum likelihood estimation under two censoring schemes, design of optimal removal plans, bootstrap calibration of control limits, and run-length analysis under parameter shifts.

## Motivation

In life testing, units are often removed from the test before failing, to save time or free up equipment. This produces censored data, where some lifetimes are known exactly and others are only known to exceed a removal time. A monitoring chart built for complete data will misjudge the process under censoring, so both estimation and chart calibration need to account for how the data was collected.

## Censoring Schemes

- **Type-I progressive censoring (PC):** at fixed inspection times $\tau_1, \dots, \tau_k$, a prespecified number $R_j$ of surviving units is randomly removed. Units still running after $\tau_k$ are censored there.
- **Type-I progressive hybrid censoring (PHC):** the same scheme, but the experiment stops early once a cap on the number of failures is reached (70% of units in this study).

## Methodology

1. **Weibull fitting and goodness of fit:** fit two-parameter Weibull models to both datasets and check fit with the Kolmogorov-Smirnov test.
2. **Censored MLE:** derive the censored Weibull log-likelihood and estimate shape and scale under PC and PHC data, optimizing on the log scale for stability.
3. **Optimal censoring design:** estimate the Fisher information of each candidate removal plan by Monte Carlo, then compare plans using A-optimality (trace of the inverse information) and D-optimality (determinant of the inverse information).
4. **Bootstrap control limits (Phase I):** estimate in-control parameters from reference data and choose the CUSUM decision interval $h$ by parametric bootstrap so that the in-control average run length (ARL₀) is close to 370.
5. **Run-length analysis:** simulate in-control ARL₀ and SDRL₀, and out-of-control ARL₁ under shifts in the scale ($\beta_1 = \delta \beta_0$) and shape ($\alpha_1 = \gamma \alpha_0$) parameters.
6. **Real-data charts (Phase II):** apply the calibrated CUSUM to the held-out portion of each dataset.

The CUSUM statistic accumulates the Weibull log-likelihood ratio for an upward scale shift:

$$S_i = \max(0,\ S_{i-1} + W_i - k), \qquad W_i = \log \frac{f(x_i; \alpha, \beta_1)}{f(x_i; \alpha, \beta_0)}$$

and signals when $S_i > h$.

## Datasets

| Dataset | Source | n | Fitted shape α̂ | Fitted scale β̂ | KS p-value |
|---|---|---|---|---|---|
| Guinea pig survival times | Bjerkedal (1960) | 72 | 1.730 | 103.49 | 0.109 |
| Ball bearing failure times | Lawless (1982) | 23 | 2.102 | 81.87 | 0.617 |

Both datasets are consistent with a Weibull model at the 5% level.

## Key Results

- **Censored MLE:** under both PC and PHC, the estimates stay close to the full-data fit (for example, α̂ = 1.58 and β̂ = 100.27 under PC for the guinea pig data, versus 1.73 and 103.49 from the complete data).
- **Optimal design:** equal removals at each inspection time minimized both A- and D-optimality criteria for both datasets, outperforming front-heavy and back-heavy plans.
- **Bootstrap calibration:** decision intervals of h = 2.05 and h = 2.57 achieved ARL₀ of about 391 and 382, close to the target of 370.
- **Scale shift detection:** the chart detects increases in scale quickly, with ARL₁ dropping from 85.6 at a 25% shift to 8.0 when the scale doubles.

| Scale shift δ | 1.25 | 1.50 | 1.75 | 2.00 |
|---|---|---|---|---|
| ARL₁ | 85.6 | 25.8 | 12.4 | 8.0 |

## Limitations

- The chart is designed for upward scale shifts and is essentially blind to shape increases: under shape shifts it rarely signals within the simulated stream. A separate or two-sided chart would be needed to monitor the shape parameter.
- The ARL₀ / SDRL₀ table currently reuses the same simulation settings for every subgroup configuration, so its rows differ only by Monte Carlo noise.

## Figures

- `fig1_histogram_fits.png`: histograms with fitted Weibull densities
- `fig2_arl_plots.png`: ARL₁ and SDRL profiles under scale and shape shifts
- `fig3_real_data_cusum.png`: Phase II CUSUM charts on both datasets

## Requirements

Python 3 with `numpy`, `scipy`, and `matplotlib`.

## Author

Chahat Ahuja, Undergraduate Research Assistant at IIT Delhi under Prof. Neeraj Joshi
