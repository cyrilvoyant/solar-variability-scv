# Solar variability and forecastability: scientific use contract

## Authority and purpose

This contract describes `scv.m` and `forecastability.m` at source commit `8cc62e8d5fbf58e418026a63be426589b93863fd`. It documents existing code; it does not deploy a WebSemantic adapter. Read the README and CITATION.cff for the study and full attribution.

The stochastic coefficient of variation measures departure from a clear-sky irradiance curve. Forecastability F also accounts for residual autocorrelation at the requested horizon. Neither is an observed accuracy score of a trained forecasting model.

## Inputs, units and context

- `I`, `I_clr`: real irradiance vectors of equal length, aligned on the same regularly sampled time grid and in the same unit, typically W/m². The caller supplies the clear-sky curve: the functions do not calculate it from location or weather.
- `d`: finite scalar diffuse-floor fraction, 0 <= d < 1, default 0.2. Dimensionless. It scales the metric; it does not clip observations.
- `T`: finite integer number of samples per day, required for F.
- `n`: finite integer forecast horizon in samples, 1 <= n <= floor(T/2), required for F. Always report the sampling interval and convert n to elapsed time explicitly.

Declare the clear-sky model, period, sampling interval, observation quality controls and missing-data treatment. Reject complex inputs before calling. Do not mix energy with irradiance or silently resample data. Keep missing timestamps as missing positions in the regular grid: removing them changes lag meaning.

## Equations and outputs

`sCV = sqrt(2)/abs(d-1) * sqrt(mean((I-I_clr)^2))/sqrt(mean(I_clr^2))`, computed on jointly available samples. The two mean-square quantities use the same sample mask. Outputs are sCV (dimensionless) and `info`: `rmse_clear_sky`, `rms_clear_sky` (input irradiance unit), `d` and retained count `n`. This `info.n` is a sample count, not the forecast horizon argument.

For F, the residual is `I-I_clr`. For each lag from n to floor(T/2), the code computes Pearson correlation between the two available residual vectors. `rho_max = max(abs(rho))`; `F = 1+sCV*(rho_max-1)`. Additional info fields are `scv`, `rho_max`, `lag_of_max`, `rho`, `lags`. Lag values are in samples; correlations and F are dimensionless.

**sCV = 1 is a reference level, not a maximum. F is not a probability.** For real inputs and defined correlations its theoretical range is [1-sCV, 1], which can include negative values. There is no clipping. A persistent diffuse-floor sky has sCV = sqrt(2), despite being temporally constant relative to that floor. Do not call all high-sCV sites rapidly fluctuating.

## Missing values and limits

sCV removes paired NaN observations, rejects Inf, mismatched lengths, no available pairs, an all-zero available clear-sky curve and invalid d. Identical zero nighttime samples leave sCV unchanged, but the time grid still matters to F. For F, the series length must exceed floor(T/2)+2. Correlations with fewer than three pairs or zero variance are undefined; the function errors if all lags are undefined. When sCV is exactly zero, F=1 and correlation/lag-of-maximum fields are NaN, with empty correlation/lag arrays.

Use long representative records for site characterization (the study uses annual records). Short examples below are teaching checks. Clear-sky bias, sampling choices and missingness affect interpretation. F does not establish operational forecast accuracy, confidence intervals or causality. The 68-station SIAR data are not included; repository demonstrations do not reproduce that station experiment.

## Minimal MATLAB examples

```matlab
s = scv([400 500 760 200], [400 800 800 400]);
% Approximately 0.5070; synthetic irradiances in W/m².
[F, info] = forecastability(500*ones(1,48),500*ones(1,48),24,1);
% F=1, info.scv=0, info.rho_max=NaN: clear-sky identity case.
```

The second example is an explicit idealized constant clear-sky test, not a realistic daily trajectory. Do not infer field forecastability from it. A hand-assumed rho_max=0.6 gives F=1+s*(0.6-1), but that is an illustration, not correlation estimated from the four input samples.

## Citation

C. Voyant, A. Julien, M. Despotovic, G. Notton, L. A. Garcia-Gutierrez, C. F. Nicolosi, P. Blanc and J. Bright, *Stochastic coefficient of variation: Assessing the variability and forecastability of solar irradiance*, Renewable Energy 256 (2026), 123913. https://doi.org/10.1016/j.renene.2025.123913

Software concept: https://doi.org/10.5281/zenodo.23562553 ; v1.0.0 archive: https://doi.org/10.5281/zenodo.23562554 . Cite the exact commit or version used; these later documentation files are not automatically in the existing archive. See repository LICENSE for software terms.
