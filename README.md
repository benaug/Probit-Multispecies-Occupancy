# Probit-Multispecies-Occupancy
Probit Multispecies Occupancy models using the parameter expanded data augmentation approach of Dorazio et al. 2025

This is a bare bones example of fitting probit multispecies occupancy models in Nimble using parameter expansion as proposed by 
Dorazio et al. (2025). Dorazio et al. follows up on Tobler et al. (2019) who used the model without parameter expansion and found it can mix very
poorly and have convergence problems, especially as the number of species and sites increase. Parameter expansion does really improve
our ability to use this model reliably.

Dorazio et al.:
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70168

Tobler et al.:
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2754

Like Dorazio, I am using the parameter expanded data augmentation (PXDA) through the Huang-Wand prior for the expanded precision matrix.
I implemented Dorazio's w and z full conditional updates in Nimble but the z full conditionals are wrong because they assume independence, 
so I am not using that approach. Therefore, I don't include the z full conditionals in this repository. This mistake produced the results
in the manuscript that the estimates of the correlation parameters were shrunk towards 0. This does not happen when using a correct
MCMC algorithm, or at least it is much reduced (only looked at 1 simulation scenario so far).

I include 3 approaches:

1. The Tobler et al. data augmentation version in files with "DA Wishart".

2. A parameter expanded data augmentation version using the Huang Wand prior that Dorazio uses. Here, I use the default block 
RW updates for the multivariate normal at each site (same as Tobler et al.), but with 10 tries per iteration instead of 1 (you can modify this).
Multiple tries really does improve mixing, but I am unsure what would be the optimal number of tries.
These files have "PXDA HuangWand" in the file name.

3. A parameter expanded data augmentation version using the Huang Wand prior that Dorazio uses. Here, I propose the multivariate normals
at each site from their full conditionals like Dorazio does, but condition on the current occupancy states. Then, I use a simple RW update for each multivariate
normal indices separately that allows the occupancy states to change. One could also just retain the block RW updates.
These files have "PXDA HuangWand V2" in the file name.

I've done a single simulation scenario and for both 2 and 3, the parameter estimates and HPD coverage look good except for some rare cases of nonconvergence which
is due to the likelihood becoming flat when the expected occupancy for a species approaches 1. I will try to improve this further by introducing
constraints to prevent this from happening. Also, version 1 is actually much more efficient than version 2, at least in the scenario I considered
with 10 species and 250 sites. The truncated MVN full conditionals mix better per iteration, but are massively slower, resulting in a net loss of
the effective sample size per unit time.

Models are currently set up with occupancy and detection intercepts only. I will add them later, but be careful if you try to do this yourself,
the full conditional updates for B will need to be modified or not used if you add covariate effects.

Dorazio's code can be found here: https://github.com/RobertDorazio/MultispeciesOccupancy/tree/main