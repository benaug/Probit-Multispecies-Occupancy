# Probit-Multispecies-Occupancy
Probit Multispecies Occupancy models using the parameter expanded data augmentation approach of Dorazio et al. 2025

This is a bare bones example of fitting probit multispecies occupancy models in Nimble largely following the approach proposed by 
Dorazio et al. 2025. Like Dorazio, I am using the parameter expanded data augmentation (PXDA) through the Huang-Wand prior for the expanded precision matrix.
I implemented Dorazio's w and z updates in Nimble (not in this repository) but I think the z full conditionals are wrong because they assume independence.
Instead, I provide a version that uses Dorazio's w updates conditioned on the z states and then use a MH update for w's where the z states are not observed.
I will improve on this later.

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70168

This paper follows up on Tobler et al. (2019) who used the model without parameter expansion which can mix very poorly especially as the number of species increases:

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2754

I provide the Tobler et al. data augmentation version in files with "DA Wishart" and the Dorazio parameter expanded 
data augmentation version in files with "PXDA HuangWand" since Dorazio used the Huang-Wand hierarchical prior for the 
precision matrix. "PXDA HaungWand" implements PXDA with a block random walk for w which is an improvement over Tobler et al.,
but w still does not mix well. "PXDA HaungWand V2" also includes Dorazio's full conditional w updates for observed z states
and independent RW updates for z's that are latent to allow z states to change.

Models are currently set up with intercepts only. Covariates can be added by following Dorazio.

Finally, I haven't formally tested this code, but it appears to be working correctly. 
The regular DA version can take a very long time to converge or may not converge. 
I expect PXDA improves mixing more with more species. 
Generally, you will need a lot of sites for multispecies occupancy to work well.

I haven't used Dorazio's code, but it can be found here: https://github.com/RobertDorazio/MultispeciesOccupancy/tree/main