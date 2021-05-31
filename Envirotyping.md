
# **E-GP: Enviromic Assembly**

author: Germano Costa-Neto

date: May 2021

______________________________________________________________________________



> As a toy-example, we made available a set of 100 genotypes over 4 environments.

#### **Functions and data**

> Functions for Enviromic Assembly came from the [EnvRtype-package](https://github.com/allogamous/EnvRtype)


```{r, eval=FALSE}

library(devtools)
install_github('allogamous/EnvRtype')
require(EnvRtype)

# and other functions
load('plot_T.RData')
load('phenoMET.RData')
load('K_G.RData')
load('leafcurve.RData')
load('apar.RData')
load('covbt.RData')
load('SoilClass.RData')
```




#### **Remote data collection of daily weather**

```{r, eval=FALSE}
env.i = c("1_AN","1_PI","2_AN","2_PI")   # id of locations
lat = c(-22.875,-22.705,-22.875,-22.705) # latitude coordinates
lon = c(-47.997,-47.637,-47.997,-47.637) # longitude coordinates
plant.date = c("2016-01-26","2016-01-21","2017-01-12","2017-01-10") # planting dates
harv.date = c('2016-08-01',"2016-07-14","2017-07-25","2017-07-15")  # harvest datte

df.clim <- get_weather(env.id = env.i,lat = lat,lon = lon,start.day = plant.date,end.day = harv.date,country = 'BRA') 

```

> The process of collecting and processing T_matrix were detailed in the previous markdown file (E-GP: Enviromic Assembly).


#### **Processing raw-data**

> * Data processing are done through three functions from EnvRtype (computation of radiation parameters, temperature-related factors and atmospheric demands).

```{r, eval=FALSE}
df.clim = param_radiation(env.data = df.clim, merge = TRUE)
df.clim = param_temperature(env.data = df.clim,Tbase1 = 8,Tbase2 = 45,Topt1 = 30,Topt2 = 37,merge = T)
df.clim = param_atmospheric(env.data = df.clim,merge=TRUE)
```

> * Additionally, leaf.cuve function estimates a canopy growing of a grain crop (e.g., maize) based on a 3-segmented function as a funciton of growing degree-days (GDD) or days after emergence. For that, the user must provide at least three paremeters (initial values, maximum values and end values). We use the same function to estimate the leaf area index (LAI) and the crop coefficient of evaptoranspiration (kc)

```{r, eval=FALSE}
require(plyr)
require(reshape2)

# crop coefficient as a function of end kc (0.3, 1.2 and 0.35) in days after emergence
kc=ddply(df.clim,.(env),function(x) leaf.curve(.endKc = c(.3,1.2,.35),
                                             .lenghPhase = c(15,36,66,90),.dae = x$daysFromStart,.digits=2))

# leaf area index (LAI) across crop development, LAI initial = 0.7; LAI max = 3.0; LAI end = 2
LAI =ddply(df.clim,.(env),function(x) leaf.curve(.endKc = c(.7,3,2),
                                               .lenghPhase = c(15,36,66,90),.dae = x$daysFromStart,.digits=2))

```


> * Then, we can estimate crop evapotranspiration and the aparent photosynthetic radiation

> * from raster files, the cov.by.trial function can be used to collect point-values based on latitude and longitude. As an example, we used a raster file containing soil classification in Brazil.


```{r, eval=FALSE}

df.clim$ET0 =df.clim$ETP*kc$kc # crop evapontranspiration ET0

df.clim$PETP = df.clim$ET0-df.clim$PRECTOT # correct ET0 - PREC

df.clim$aPAR = aPAR(Srad = df.clim$SRAD,IAF=IAF$kc) # aparent photosyntentic radiation

# collecting additional informatiom from raster files (e.g., soil class)
solo=cov.by.trial(cov.raster = soil,reference = df.clim,long = 'LON',latd = 'LAT',trial = 'env')

df.clim = extract_GIS(covraster = soil,Latitude = 'LAT',Longitude = 'LON',env.data = df.clim,env.id = 'env',name.out = 'soil.class')

```

#### **Creating environmental typologies**

> First create the genotype-environment kernel by the kronecker product between K_ET and K_G;

```{r, eval=FALSE}
(id.variables =names(df.clim)[c(2,10:16,21:24,26:29,31:34)] )# transitory effects
(interval = c(0,14,35,65,90,120)) # interval of crop development (days after emergence)
```


> * Enviromic Assembly of environmental covariates (ECs) using statistic criteria

```{r, eval=FALSE}
T_matrix = env_typing(env.data = df.clim, quantiles = c(.25,.50,.75),var.id = id.variables,env.id = 'env',time.window = interval,by.interval = T,format = 'wide') # without cardinals, only statistical criteria

```


> * For some ECs, its possible to include some degree of ecophysiology knowledge by describing cardinal limits.


```{r, eval=FALSE}
cardinals = list(ALT = NULL, T2M = c(0,9,22,32,37,45),T2M_MAX =c(0,9,22,32,37,45),
                 T2M_MIN = c(0,9,22,32,37,45), PRECTOT =c(0,5,10,20,Inf),WS2M = NULL, 
                 RH2M = NULL, T2MDEW = NULL,n = NULL, N = NULL, RTA=NULL, SRAD=NULL,
                 FRUE = c(0,.25,.50,.75,1),T2MR = NULL,VPD=NULL,SPV=NULL,
                 PETP=NULL,ET0=c(0,6,10,15,Inf),aPAR = NULL,soil=NULL)

T_matrix =env_typing(env.data = df.clim, cardinals = cardinals, quantiles = c(.25,.50,.75),var.id = id.variables,env.id = 'env',time.window = interval,by.interval = T,format = 'wide') 


T_matrix = T_matrix[,-which(apply(T_matrix,2,var) == 0)] # you can remove colunmns with no information
panelEnvAssembly(T_matrix) # plot a panel of enviromic assembly

```


#### **Environmental Similarity Kernel**

```{r, eval=FALSE}
Panel=env_kernel(T_matrix)
panelEnvAssembly(Panel$varCov) # ECs relatedness
panelEnvAssembly(Panel$envCov) # Environmental Similarity

K_E = Panel$envCov # for predictive purposes, this is the env. kinship
```


#### **Genomic Prediction with BGGE Package**

>* The genomic prediction models can be run in several R packages. Here we exemplify the use in BGGE. For that, we used the get_kernel function of the EnvRtype package.

```{r, eval=FALSE}

require(BGGE)

Kernels = get_kernel(K_E = NULL,K_G = list(G=K_G),Y = phenoMET,model = 'MDs') # Benchmark GBLUP
Kernels = get_kernel(K_E = list(ET=K_E),K_G = list(G=K_G),Y = phenoMET,model = 'EMDs') # E-GP block diagonal
Kernels = get_kernel(K_E = list(ET=K_E),K_G = list(G=K_G),Y = phenoMET,model = 'RNMM') # E-GP with reaction-norm
Kernels = get_kernel(K_E = list(ET=K_E),K_G = list(G=K_G),Y = phenoMET,model = 'RNMDs') # E-GP block diagonal+reaction norm


fit <- BGGE(y = phenoMET$value,
            K = Kernels, # put here the kernel prepared using get_kernel
            tol = 1e-20,
            XF = model.matrix(~0+env,phenoMET),
            ne = as.vector(table(phenoMET$env)),
            verbose = FALSE)
```

