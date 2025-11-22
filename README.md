# Probit-Multispecies-Occupancy
Probit Multispecies Occupancy models using the parameter expanded data augmentation approach of Dorazio et al. 2025

This is a bare bones example of fitting probit multispecies occupancy models in Nimble largely following the approach proposed by 
Dorazio et al. 2025. Like Dorazio, I am using the parameter expanded data augmentation (PXDA) through the Huang-Wand prior for the expanded precision matrix.
I implemented Dorazio's w and z updates in Nimble (not in this repository) but I think the z full conditionals are wrong because they assume independence.
Instead, I provide a version (V2) that uses Dorazio's w updates conditioned on the z states and then use a MH update for w's where the z states are not observed.
I will improve on this later.

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70168

This paper follows up on Tobler et al. (2019) who used the model without parameter expansion which can mix very poorly especially as the number of species increases:

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2754

I provide the Tobler et al. data augmentation version in files with "DA Wishart" and the Dorazio parameter expanded 
data augmentation version in files with "PXDA HuangWand" since Dorazio used the Huang-Wand hierarchical prior for the 
precision matrix. "PXDA HaungWand" implements PXDA with a block random walk for w which is an improvement over Tobler et al.,
but w still does not mix well. "PXDA HaungWand V2" also includes Dorazio's full conditional w updates for observed z states
and independent RW updates for z's that are latent to allow z states to change.

Models are currently set up with occupancy and detection intercepts only. I will add them later, but be careful if you try to do this yourself.

I have PXDA set up with full conditional updates for B, but the custom conjugate update is only coded for the intercept only model.
For 1 data set, the mixing was similar, but the full conditionals are faster to compute. I noticed a case of one of the betas drifting
off towards positive infinity using an RW update with a very diffuse prior. This shouldn't happen with the conjugate update, or you can
set bounds on the scale of w that keep mean psi between, say, 0.001 and 0.999.

Finally, I haven't formally tested this code, but it appears to be working correctly. 

I haven't used Dorazio's code, but it can be found here: https://github.com/RobertDorazio/MultispeciesOccupancy/tree/main