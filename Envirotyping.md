
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
require(EnvRtype)
require(reshape2)
require(plyr)
require(superheat)

# source data

```




#### **Remote data collection of daily weather**

```{r, eval=FALSE}
env.i      = paste0('Env',1:4)  # id of locations
lat        = c(-22.875,-22.705,-22.875,-22.705) # latitude coordinates
lon        = c(-47.997,-47.637,-47.997,-47.637) # longitude coordinates
plant.date = c("2016-01-26","2016-01-21","2017-01-12","2017-01-10") # planting dates
harv.date  = c('2016-08-01',"2016-07-14","2017-07-25","2017-07-15")  # harvest datte

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
kc = ddply(df.clim,.(env),function(x) 
  leaf.curve(.endKc      = c(.3,1.2,.35),
             .lenghPhase = c(15,36,66,90),
             .dae = x$daysFromStart,.digits=2))

# leaf area index (LAI) across crop development, LAI initial = 0.7; LAI max = 3.0; LAI end = 2
LAI = ddply(df.clim,.(env),function(x) 
      leaf.curve(.endKc      = c(.7,3,2),
                 .lenghPhase = c(15,36,66,90),
                 .dae = x$daysFromStart,.digits=2))


```


> * Then, we can estimate crop evapotranspiration and the aparent photosynthetic radiation

> * from raster files, the cov.by.trial function can be used to collect point-values based on latitude and longitude. As an example, we used a raster file containing soil classification in Brazil.


```{r, eval=FALSE}

df.clim$ET0 =  df.clim$ETP*kc$kc # crop evapontranspiration ET0

df.clim$PETP = df.clim$ET0-df.clim$PRECTOT # P-ET0

df.clim$aPAR = aPAR(Srad = df.clim$SRAD,K = 0.5,fraction = 0.5,LAI = LAI[,3]) # aparent photosyntentic radiation


df.clim = extract_GIS(covraster = soil,Latitude = 'LAT',Longitude = 'LON',env.data = df.clim,env.id = 'env',name.out = 'soil.class')



```

#### **Creating environmental typologies**

> First create the genotype-environment kernel by the kronecker product between K_ET and K_G;

```{r, eval=FALSE}
(id.variables = names(df.clim)[c(11:17,22:25,27:30,32:34,2:3)] )# transitory effects
(interval = c(0,14,35,65,90,120)) # interval of crop development (days after emergence)
```


> * Enviromic Assembly of environmental covariates (ECs) using statistic criteria

```{r, eval=FALSE}
(id.variables = names(df.clim)[c(11:17,22:25,27:30,32:34,2:3)] )# transitory effects
(interval = c(0,14,35,65,90,120)) # interval of crop development (days after emergence)
T_matrix = env_typing(env.data = df.clim, quantiles = c(.25,.50,.75),var.id = id.variables,env.id = 'env',time.window = interval,by.interval = T,format = 'wide') # without cardinals, only statistical

```


> * For some ECs, its possible to include some degree of ecophysiology knowledge by describing cardinal limits.


```{r, eval=FALSE}

(id.variables = names(df.clim)[c(11:17,22:25,27:30,32:34,2:3)] )# transitory effects
(interval = c(0,14,35,65,90,120)) # interval of crop development (days after emergence)

## add cardinals. If you dont know for sure, put = NULL and it will run the statistical approach
cardinals = list(T2M = c(0,9,22,32,37,45),T2M_MAX =c(0,9,22,32,37,45),
                 T2M_MIN = c(0,9,22,32,37,45), PRECTOT =c(0,5,10,20,Inf),WS2M = NULL, 
                 RH2M = NULL,T2MDEW = NULL,n = NULL, N = NULL, RTA=NULL, SRAD=NULL,
                 FRUE = c(0,.25,.50,.75,1), T2M_RANGE = NULL,VPD=NULL,SPV=NULL,
                 PETP=NULL,ET0=c(0,6,10,15,Inf),aPAR = NULL,soil=NULL,alt=NULL)

names(df.clim)
T_matrix =env_typing(env.data = df.clim, cardinals = cardinals, quantiles = c(.25,.50,.75),var.id = id.variables,env.id = 'env',time.window = interval,by.interval = T,format = 'wide') 

T_matrix <- T_matrix[,apply(T_matrix,2,var) > 0]
panel_EnvAssembly(ECs = T_matrix,order.row = T,order.col = T)


```


#### **Environmental Similarity Kernel**


> **T_matrix** of environmental typologies

```{r, eval=FALSE}
  
(id.variables = names(df.clim)[c(11:17,22:25,27:30,32:34,2:3)] )# transitory effects
(interval = c(0,14,35,65,90,120)) # interval of crop development (days after emergence)

cardinals = list(T2M = c(0,9,22,32,37,45),T2M_MAX =c(0,9,22,32,37,45),
                 T2M_MIN = c(0,9,22,32,37,45), PRECTOT =c(0,5,10,20,Inf),WS2M = NULL, 
                 RH2M = NULL,T2MDEW = NULL,n = NULL, N = NULL, RTA=NULL, SRAD=NULL,
                 FRUE = c(0,.25,.50,.75,1), T2M_RANGE = NULL,VPD=NULL,SPV=NULL,
                 PETP=NULL,ET0=c(0,6,10,15,Inf),aPAR = NULL,soil=NULL,alt=NULL)


T_matrix =env_typing(env.data = df.clim, cardinals = cardinals, quantiles = c(.25,.50,.75),var.id = id.variables,env.id = 'env',time.window = interval,by.interval = T,format = 'wide') 

T_matrix <- T_matrix[,apply(T_matrix,2,var) > 0]
panel_EnvAssembly(ECs = T_matrix,order.row = T,order.col = T)
dim(T_matrix)


Panel_T = env_kernel(env.data = T_matrix) # building the environmenal relatedness (linear kernel)
panel_EnvAssembly(ECs = Panel_T$varCov,bottom.size = 3,left.size = 3,order.row = T,order.col = T) # ECs relatedness
panel_EnvAssembly(ECs = Panel_T$envCov,bottom.size = 3,left.size = 3,order.row = T,order.col = T) # Environmental Similarity


K_E = Panel$envCov # for predictive purposes, this is the env. kinship
```

> **W_matrix** of quantitative covariables

```{r, eval=FALSE}
 
W_matrix = W_matrix(env.data = df.clim, statistic = 'quantile',probs = c(.25,.50,.75),var.id = id.variables,env.id = 'env',time.window = interval,by.interval = T) 
W_matrix <- W_matrix[,apply(T_matrix,2,var) > 0]
panel_EnvAssembly(ECs = W_matrix,order.row = T,order.col = T)

Panel_W = env_kernel(env.data = W_matrix) # building the environmenal relatedness (linear kernel)
panel_EnvAssembly(ECs = Panel_W$varCov,bottom.size = 3,left.size = 3,order.row = T,order.col = T) # ECs relatedness
panel_EnvAssembly(ECs = Panel_W$envCov,bottom.size = 3,left.size = 3,order.row = T,order.col = T) # Environmental Similarity

K_E = Panel_T$envCov # using T matrix as example

```

> **Environmental relatedness**

```{r, eval=FALSE}
K_E = Panel_T$envCov # using T matrix as example

Kernels = get_kernel(K_E = NULL,K_G = list(G=K_G),        data = phenoMET,model = 'MDs') # Benchmark GBLUP
Kernels = get_kernel(K_E = list(ET=K_E),K_G = list(G=K_G),data = phenoMET,model = 'EMDs') # E-GP block diagonal (BD)
Kernels = get_kernel(K_E = list(ET=K_E),K_G = list(G=K_G),data = phenoMET,model = 'RNMM') # E-GP with reaction-norm

```

#### **Genomic Prediction with (with Bayesian approach from BGGE package)**

>* The genomic prediction models can be run in several R packages. Here we exemplify the use in BGGE. For that, we used the get_kernel function of the EnvRtype package.

```{r, eval=FALSE}

# example: running E-GP-BD
Kernels = get_kernel(K_E = list(ET=K_E),K_G = list(G=K_G),data = phenoMET,model = 'EMDs') # E-GP with reaction-norm

names(phenoMET)
y   = 'value'
gid = 'gid'
env = 'env'

fixed = model.matrix(~0+env,phenoMET) # you can put it or not

require(BGGE)

fit <- kernel_model(y = y,env = env,gid = gid,data = phenoMET,
                    random = Kernels,fixed = fixed,iterations = 15000,thining = 2000,tol = 10)

fit$yHat   # predicted values

```

