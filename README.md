# Probit-Multispecies-Occupancy
Probit Multispecies Occupancy models using parameter expanded data augmentation approach of Dorazio et al. 2025 (mostly)

This is a bare bones example of fitting probit multispecies occupancy models in nimble using the approach proposed by Dorazio et al. 2025.

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70168

This paper follows up on Tobler et al. (2019) who used the model without parameter expansion which can mix very poorly especially as the number of species increases:

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2754

I provide the Tobler et al. data augmentation version in files with "DA Wishart" and the Dorazio parameter expanded data augmentation
version in files with "PXDA HuangWand"
since Dorazio used the Huang-Wand hierarchical prior for the precision matrix. Nimble is able to provide the conjugate
update for the precision matrix and I use a custom update for the conjugate update of the auxiliary parameters.

Currently set up with intercepts only. Covariates can be added by following Dorazio.

Dorazio features currently not incorporated:
1) Full conditional updates for w and z (I'm using block RW on w/z together)
2) Full conditional updates for B (I'm using RW or slice)
3) Full conditional updates for p (detection) (I'm using RW)

I think this parameter expansion and Huang-Wand prior get us most of the way to the efficiency of Dorazio's approach.
Need a truncated multivariate Normal random vector generator for w full conditionals. An advantage of nimble is that it is
fully compiled to c++.

Finally, I haven't formally tested this code, but it appears to be working correctly. The regular DA version can take a
very long time to converge or may not converge. I expect PXDA improves mixing more with more species.
Generally, you will need a lot of sites for multispecies occupancy to work well.

I haven't tried Dorazio's code, but it can be found here:
https://github.com/RobertDorazio/MultispeciesOccupancy/tree/main