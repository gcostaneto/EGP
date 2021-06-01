# **Enviromic Assembly for purposes of Genomic Prediction (E-GP)**

Author: Germano Costa Neto <germano.cneto@gmail.com> 

Last update: 30th may 2021


# Main

This codes are part of the paper *Enviromic assembly increases accuracy and reduces costs of the genomic prediction for yield plasticity*, (Costa-Neto et al, under review)

**Abstract**: Quantitative genetics states that phenotypic variation is a consequence of genetic and environmental factors and their subsequent interaction. Here, we present an enviromic assembly approach, which includes the use of ecophysiology knowledge in shaping environmental relatedness into whole-genome predictions (GP) for plant breeding (referred to as E-GP). We propose that the quality of an environment is defined by the core of environmental typologies (envirotype) and their frequencies, which describe different zones of plant adaptation. From that, we derive markers of environmental similarity cost-effectively. Combined with the traditional genomic sources (e.g., additive and dominance effects), this approach may better represent the putative phenotypic variation across diverse growing conditions (i.e., phenotypic plasticity). Additionally, we couple a genetic algorithm scheme to design optimized multi-environment field trials (MET), combining enviromic assembly and genomic kinships to provide in-silico realizations of the future genotype-environment combinations that must be phenotyped in the field. As a proof-of-concept, we highlight E-GP applications: (1) managing the lack of phenotypic information in training accurate GP models across diverse environments and (2) guiding an early screening for yield plasticity using optimized phenotyping efforts. Our approach was tested using two non-conventional cross-validation schemes to better visualize the benefits of enviromic assembly in sparse experimental networks. Results on tropical maize show that E-GP outperforms benchmark GP in all scenarios and cases tested. We show that for training accurate GP models, the genotype-environment combinations' representativeness is more critical than the MET size. Furthermore, we discuss theoretical backgrounds underlying how the intrinsic envirotype-phenotype covariances within the phenotypic records of (MET) can impact the accuracy of GP and limits the potentialities of predictive breeding approaches. The E-GP is an efficient approach to better use environmental databases to deliver climate-smart solutions, reduce field costs, and anticipate future scenarios. 

# Data sets

## Multi-Regional Set (Helix Seeds)

> 452 tropical maize single-crosses provided by Helix Sementes®, Brazil. Hybrids were obtained from crosses between 128 inbred lines and were evaluated for grain yield (GY) and plant height (PH). Field trials were carried out using a randomized complete block design with two replicates each, allocated across five sites for GY and three for PH during the growing season of 2015.

> Inbred lines were genotyped via the Affymetrix® Axiom® Maize Genotyping Array (Unterseer et al. 2014) with 660K SNP markers. Quality control for SNPs was made based on call rate (CR), in which all markers with any missing data point were excluded, and minor allele frequency (MAF) procedures, in which markers with a low level of polymorphism (MAF < 0.05) were removed. Hybrid genotypes were scored by an allelic combination of homozygous markers of parental lines. After quality control, 37,625 SNP were used to compare the imputation methods.

Check the [CIMMYT DATA VERSE Repository](https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10887) for the full genotypic and phenotypic data


<div id="p2.2" />

## N-level data set (USP set)

> 906 maize single-crosses obtained from a full dial- lel, according to Griffing’s method 4, divided into two heterotic groups, flint and dent, with 34 and 15 lines, respec- tively. Moreover, each heterotic group has a representative line, frequently used as the tester in our breeding program.

> The experimental scheme used to evaluate the hybrids was an augmented block design (unreplicated trial) consisted of small blocks, each with 16 unique hybrids and two checks. Trials were carried out in Anhembi (22°50′51′′S, 48°01′06′′W, 466 m) and Piracicaba, at São Paulo State, Brazil (22°42′23′′S, 47°38′14′′W, 535 m), during the second growing season of 2016 and 2017, cultivated between January to June. In both sites and years, the hybrids were evaluated under two nitrogen (N) levels, low (LN) with 30 kg N ha−1, and normal (NN) with 100 kg N ha−1.

> The genotyping of the 49 tropical inbred lines was per- formed by Affymetrix® platform, containing about 614,000 SNPs (Unterseer et al. 2014). Then, markers with low call rate (< 95%), minor allele frequency (MAF < 0.05) and heterozygous loci on at least one individual were removed. The missing markers were imputed using the [snpReady](https://github.com/italo-granato/snpReady) R package. Finally, the resulting 146,365 SNPs high-quality polymorphic SNPs were used to build the artificial hybrids genomic matrix, deduced by combining the genotypes from its two parents.

Check the [Mendeley Repository](https://data.mendeley.com/datasets/tpcw383fkm/3) for the full genotypic and phenotypic data

## Data Availability For Running the Codes

> Both data sets can download in R.Data and RDS format direct [here](https://github.com/gcostaneto/KernelMethods/tree/master/Heredity%20Data%20Set) and [Mendeley Repository](https://data.mendeley.com/datasets/cxkzb8mr8b/1).
  
Then, to import and read the file in R you must to use the readRDS() function (more detail below). To download every dataset (USP and Helix), please use the following [file](https://github.com/gcostaneto/KernelMethods/blob/master/Heredity%20Data%20Set/Full%20Data.zip).

# Software

# Pipelines

> **[Envirotyping: data collection, processsing and essembly](https://github.com/gcostaneto/EGP/blob/main/Envirotyping.md)**
> **[Genetic Algorithm for selective phenotyping with genomics and enviromics](https://github.com/gcostaneto/EGP/blob/main/Selective%20Phenotyping.md)**
> **Creating different GxE scenarios** (in prep)
> **Screening yield adaptability from GBLUP-based predictions** (in prep)


