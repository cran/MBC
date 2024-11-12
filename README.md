### Description

[MBC](https://cran.r-project.org/package=MBC) is an R package for calibrating
and applying univariate and multivariate bias correction algorithms for climate
model simulations of multiple climate variables. Three iterative multivariate
methods are supported: (i) MBC Pearson correlation (`MBCp`), (ii) MBC rank
correlation (`MBCr`), and (iii) MBC N-dimensional probability density function
transform (`MBCn`). The first two, `MBCp` and `MBCr` (Cannon, 2016), match
marginal distributions and inter-variable dependence structure. Dependence
structure can be measured either by the Pearson correlation (`MBCp`) or by the
Spearman rank correlation (`MBCr`). The energy distance score (`escore`) is
recommended for model selection. The third, `MBCn` (Cannon, 2018), which
operates on the full multivariate distribution, is more flexible and can be
considered to be a multivariate analogue of univariate quantile mapping. All
aspects of the observed distribution are transferred to the climate model
simulations. In each of the three methods, marginal distributions are corrected
by the univariate change-preserving quantile delta mapping (`QDM`) algorithm
(Cannon et al., 2015). Finally, an implementation of the Rank Resampling for
Distributions and Dependences (`R2D2`) method introduced by Vrac (2018) is also
included.

### References 

Cannon, A.J., 2018. Multivariate quantile mapping bias correction: An
N-dimensional probability density function transform for climate model
simulations of multiple variables. Climate Dynamics, 50(1-2):31-49.
doi:10.1007/s00382-017-3580-6

Cannon, A.J., 2016. Multivariate bias correction of climate model output:
Matching marginal distributions and inter-variable dependence structure. Journal
of Climate, 29:7045-7064. doi:10.1175/JCLID-15-0679.1

Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction of simulated
precipitation by quantile mapping: How well do methods preserve relative changes
in quantiles and extremes? Journal of Climate, 28:6938-6959.
doi:10.1175/JCLI-D-14-00754.1

Francois, B., M. Vrac, A.J. Cannon, Y. Robin, and D. Allard, 2020. Multivariate
bias corrections of climate simulations: Which benefits for which losses? Earth
System Dynamics, 11:537-562. doi:10.5194/esd-11-537-2020

Vrac, M., 2018. Multivariate bias adjustment of high-dimensional climate
simulations: the Rank Resampling for Distributions and Dependences (R2D2) bias
correction. Hydrology and Earth System Sciences, 22:3175-3196.
doi:10.5194/hess-22-3175-2018