# Spatial-Mark-Resight
Various spatial mark resight MCMC samplers in nimble allowing for all sample types including "unknown marked status". Handles known and unknown number of marked individuals and optional telemetry.

If you want to use individual ID covariates, there is a different Github project for that: Spatial-Mark-Resight-IDcov.

This code contains contributions from Glenn Stauffer (though any problems should be blamed on me).

11/27/24:

I added

1) "DA2" version with a Poisson observation model for single sessions. This uses a reversible-jump MCMC approach that was already being used in the Multisession versions.
In this version, the prior distribution for N is Poisson instead of Binomial, and you can update only a subset of the z's on each iteration, say 25%. I can make a negative binomial version of this, or generally any count model can be used.

2) "Poisson DA2 Marginal", which marginalizes over individual ID instead of updating the latent capture history. This includes results from Herliansyah et al.
to speed up the N/z and s updates. This approach only works with Poisson count models (not negative binomial, etc.)

https://link.springer.com/article/10.1007/s13253-023-00598-3

12/15/24:

Added "Poisson Dcov DA2 Marginal" that allows spatial density covariates.