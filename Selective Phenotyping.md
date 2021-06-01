
# **E-GP: selective phenotyping**

author: Germano Costa-Neto

date: May 2021

______________________________________________________________________________



> As a toy-example, we made available a set of 100 genotypes over 4 environments.

#### **Functions and data**

> Functions for Enviromic Assembly came from the [EnvRtype-package](https://github.com/allogamous/EnvRtype). The Supplementary data can be downloaded [here](https://github.com/gcostaneto/EGP/blob/main/Genomic_Enviromic_Example%202.rar)


```{r, eval=FALSE}
# packages
require(STPGA)
require(superheat)

# source codes
source('https://raw.githubusercontent.com/gcostaneto/EGP/main/supplementary_codes_2.R')

# supplementary data
load('Genomic_Enviromic.RData')  # genomic relatedness for A effects (100 x 100)



```


#### **Definition of the Effective Number of Observations**

> First create the genotype-environment kernel by the kronecker product between K_ET and K_G;

```{r, eval=FALSE}
superheat(K_E) # environmic kernel builted for some 4 environments
superheat(K_G) # genomic kernel with 100 genotypes

K_GE = kronecker(K_E,K_G,make.dimnames = T)
colnames(K_GE) # names of genotypes x environment combinations

```

> Then idenditify the effective number of genotype-environment combinations.
> This algoritm is inspired in Misztal (2016), in which the effect number of genotypes using only K_G.
> By default, fraction is = 98% (0.98).

```{r, eval=FALSE}
SVD = N_GE(Kernel = K_GE,plot = F,fraction = .98,svd.print = T)
ntrain = SVD$Ne        # N_GE
svdGEI = K_GE %*% SVD$svd.v # SVD-based GEI matrix for optimize the genetic algorithm
rownames(svdGEI) = rownames(K_GE)
```


#### **Genetic Algorithm for Selecting MET Subsets**

**Parameters**:

> 100 iteraions (niterations);

> 80% of mutation for each generated solution (mutprob = 0.80)

> 5 solutions selected as elite parents for next generations (nelite = 5);

> PEVmean criteria (errorstat = 'PEVMEAN')

```{r, eval=FALSE}
TS <- GenAlgForSubsetSelectionNoTest(P = as.matrix(svdGEI), ntoselect = ntrain,
                                     nelite = 5,mutprob = .8, niterations = 100, plotiters = F, 
                                     lambda = 1e-5, errorstat = "PEVMEAN", mc.cores = 4)[[1]]


TS # training set (super-optmized) consisted by some genotypes at some environments
superheat(TS)

saveRDS(object = TS, file = paste0(ntrain,"_optmized_SET"))

```

#### **Final Considerations**

> * in the file TS *, it is presented genotype_by_environement combinations to compose a super-optimized training set for genomic prediction.

> * we recommend to check for your dataset which optimization parameters lead to better selective schemes;

> * we also strong recommend the use of different enviromic assembly according to the environmental complexity of you will be dealing with.

> * If you found some concerns or have any suggestions, please contact us.
