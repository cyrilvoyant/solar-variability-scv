# Solar variability and forecastability: scientific use contract

## Authority and purpose

This contract describes `scv.m` and `forecastability.m` at source commit `8cc62e8d5fbf58e418026a63be426589b93863fd`. It documents existing code; it does not deploy a WebSemantic adapter. Read the README and CITATION.cff for the study and full attribution.

Scientific authority checked: the author manuscript [arXiv:2506.21807v1](https://arxiv.org/abs/2506.21807v1), Sections 2–4 and 6, especially Eqs. 9–13. The publisher version has not been compared line by line. The method first separates the deterministic solar trend from stochastic departures, then uses temporal dependence to characterize forecastability: understand the resource before choosing a more complex forecasting model.

The stochastic coefficient of variation measures departure from a clear-sky irradiance curve. Forecastability F also accounts for residual autocorrelation at the requested horizon. Neither is an observed accuracy score of a trained forecasting model.

## Inputs, units and context

- `I`, `I_clr`: real irradiance vectors of equal length, aligned on the same regularly sampled time grid and in the same unit, typically W/m². The caller supplies the clear-sky curve: the functions do not calculate it from location or weather.
- `d`: finite scalar diffuse-floor fraction, 0 <= d < 1, default 0.2. Dimensionless. It scales the metric; it does not clip observations.

Section 4.2 adopts **d=0.2 for the reported metric and cross-site comparisons**. Changing d is a sensitivity analysis or explicitly labeled sCV_d, not a silent replacement of the paper's convention. The physical framework assumes irradiance between the diffuse-floor and clear-sky envelopes, neglecting rare extreme attenuation and overirradiance. The functions do not enforce those physical assumptions.
- `T`: finite integer number of samples per day, required for F.
- `n`: finite integer forecast horizon in samples, 1 <= n <= floor(T/2), required for F. Always report the sampling interval and convert n to elapsed time explicitly.

Declare the clear-sky model, period, sampling interval, observation quality controls and missing-data treatment. Reject complex inputs before calling. Do not mix energy with irradiance or silently resample data. Keep missing timestamps as missing positions in the regular grid: removing them changes lag meaning.

## Equations and outputs

`sCV = sqrt(2)/abs(d-1) * sqrt(mean((I-I_clr)^2))/sqrt(mean(I_clr^2))`, computed on jointly available samples. The two mean-square quantities use the same sample mask. Outputs are sCV (dimensionless) and `info`: `rmse_clear_sky`, `rms_clear_sky` (input irradiance unit), `d` and retained count `n`. This `info.n` is a sample count, not the forecast horizon argument.

For F, the residual is `I-I_clr`. For each lag from n to floor(T/2), the code computes Pearson correlation between the two available residual vectors. `rho_max = max(abs(rho))`; `F = 1+sCV*(rho_max-1)`. Additional info fields are `scv`, `rho_max`, `lag_of_max`, `rho`, `lags`. Lag values are in samples; correlations and F are dimensionless.

**The paper constructs a normalized [0,1] interpretation under its working assumptions; the raw implementation does not enforce that interval. F is not a probability.** The normalization uses equiprobable clear and diffuse-floor states and representative clear-sky moments for the selected subset (Section 4.1, Eqs. 6–10). Outside that reference scenario, substitution of a persistent diffuse-floor sky into Eq. 10 gives sqrt(2): this is a mathematical boundary check, not a reported experimental result. Section 4.1 and the conclusion explicitly discuss the constant-cloud-cover limitation. Do not silently clip values or present this exception as the paper's intended operating regime. For defined real correlations F lies in [1-sCV,1]; an out-of-regime sCV>1 can therefore permit negative F and calls for checking the assumptions.

## Missing values and limits

sCV removes paired NaN observations, rejects Inf, mismatched lengths, no available pairs, an all-zero available clear-sky curve and invalid d. Identical zero nighttime samples leave sCV unchanged, but the time grid still matters to F. For F, the series length must exceed floor(T/2)+2. Correlations with fewer than three pairs or zero variance are undefined; the function errors if all lags are undefined. When sCV is exactly zero, F=1 and correlation/lag-of-maximum fields are NaN, with empty correlation/lag arrays.

Use long representative records for site characterization (the study uses annual records). Short examples below are teaching checks. Clear-sky bias, sampling choices and missingness affect interpretation. F does not establish operational forecast accuracy, confidence intervals or causality. The 68-station SIAR data are not included; repository demonstrations do not reproduce that station experiment.

Section 4.1 recommends at least one year for the long-period assumptions and a minimum time step of a few minutes to reduce sensor-response and high-frequency overirradiance effects. It warns that time steps above one hour can undermine the representative-subset approximation used in normalization. These are scientific interpretation limits, not runtime checks. Section 4.4 explicitly clarifies that rho_max uses absolute autocorrelation, which matches `max(abs(rho))` in the code. F depends on horizon even though sCV does not; do not describe both as horizon-independent.

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
