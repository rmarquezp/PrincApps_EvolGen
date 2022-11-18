Visualizing Genetic Structure in Space
============

Visualizing the distribution of genetic variation across geography is a powerful tool in studies of biogeogprahy and historical demography. Today we will be using a few different approaches to visualize patterns of both discrete and continuous population structure on a map, which vary in the degree to which they explicitly consider geography in their estimation of demographic parameters. The first approach will consist of estimating admixture proportions and plotting them directly on the map, while the second will consist of directly estimating rates of migration across the landscape based off of the degree of genetic differentiation and sampling positions of a set of samples. 

## The Data

For today's practical we will be analyzing data obtained from gray wolves across Canada and Alaska by Schweizer et al ([2016](https://onlinelibrary.wiley.com/doi/10.1111/mec.13364)). The authors used a SNP array (i.e. chip) designed for dogs to genotype 111 wolves at roughly 42,000 SNPs.<br><br>
<img src="https://nywolf.org/wp-content/uploads/2022/01/Trumpet_Lighthawk_family_winter_2020-scaled-2-1400x788.jpg" width="600">
<br><br>
Before we start, you may need to install some of the following R packages: `conStruct, rcpp, RcppEigen, raster, rgeos, sp, rworldmap`. Please take a minute to install them if needed.<br><br>

Lets now explore our dataset a little. To begin, plot the collection localities on a map of North America. Use R locally on your computer. For this we will need the file called `wolves.coord`, which contains the coordinate for each one of our samples. 

```R
library(raster)
library(rworldmap)

## Download elevation and coastline/border data
altitude=getData('worldclim', var='alt', res=2.5)
map=getMap(resolution="low")

## crop the data to only include the area of interest: latitudes between 43 and 85 and longitudes between -54 and -175. 
northNA=extent(-175, -54, 43, 85)
altitude_crop=crop(altitude,northNA)
map_crop=crop(map,northNA) ##May throw a warning but seems to work fine regardless

## Load the sampling coordinates
coords=read.table("wolves.coord")

## Plot the elevation layer
plot(altitude_crop, xlim=northNA[1:2], ylim=northNA[3:4])

## Plot the coastline and borders
plot(map_crop, lwd=1, add=T, xlim=c(-175,-54), ylim=c(43,80))

## Plot the samples
points(coords, pch=21, bg="white")
```
From the map we can see that the sampling was pretty consistent across northern North America, which means we will probably not have artifacts associated to uneven sampling. 

<img src="../Images/wolfMap.png" width="600">

## Admixture proportions on a map

A commonly used way to visualize genetic structure on a map is to estimate admixture proportions and plotting them as pie chars at each collection locality. Since the data were obtained using a SNP array (instead of sequencing), we can't use `ngsAdmix` as we did perviously, since we cannot estimate genotype likelihoods. Instead, we will use [Admixture](https://dalexander.github.io/admixture/index.html), which implements the same model, but assumes genotypes are known. Log into `greatlakes` and start an interactive job with 8Gb and one processor. Then load the follwing modules: `Bioinformatics plink clang gcc intel eigen boost`, and define the following variable:<br>

```software_dir=/scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data```
<br>
The data provided by Schweizer et al are in a commonly used format for SNP arrays (and also sequencing data) called `plink` after the program it was designed for. The files are `nacanids_111indiv_unrel_noYNP_42Ksnps.tfam` and `nacanids_111indiv_unrel_noYNP_42Ksnps.tped`. Please upload them to the cluster using `scp` (<b>hint:</b> you can upload multiple files with similar names using the <i>wildcard</i> character `*`, for example `scp nacanids_111indiv_unrel_noYNP_42Ksnps.* marquezr@greatlakes....`). Once your files are on the cluster take a look inside:

```bash
head nacanids_111indiv_unrel_noYNP_42Ksnps.tped 

1 chr1.3398480 0 3398480 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3453895 0 3453895 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3649882 0 3649882 N N 2 2 N N N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N N N 2 2 2 2 2 2 2 2 N N N N N N 2 2 N N N N 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 N N 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3691446 0 3691446 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3713883 0 3713883 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3891091 0 3891091 N N N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 N N 2 2 1 2 2 2 2 2 2 2 2 2 N N 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3897610 0 3897610 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3923133 0 3923133 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3923407 0 3923407 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
1 chr1.3999542 0 3999542 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 N N 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

head nacanids_111indiv_unrel_noYNP_42Ksnps.tfam 

11226.CEL 11226.CEL 0 0 0 -9
11228.CEL 11228.CEL 0 0 0 -9
11232_CLU_NQUE-I 11232_CLU_NQUE-I 0 0 0 -9
11234_CLU_NQUE-I 11234_CLU_NQUE-I 0 0 0 -9
2042_CLU_CANADA-FBG 2042_CLU_CANADA-FBG 0 0 0 -9
2043_CLU_CANADA-FG 2043_CLU_CANADA-FG 0 0 0 -9
2048_CLU_CANADA-FW 2048_CLU_CANADA-FW 0 0 0 -9
2049_CLU_CANADA-FW 2049_CLU_CANADA-FW 0 0 0 -9
2052_CLU_CANADA-FW 2052_CLU_CANADA-FW 0 0 0 -9
2056_CLU_CANADA-FG 2056_CLU_CANADA-FG 0 0 0 -9
```

As you can see, the `.tfam` file contains sample names, while `.tped` contains genotypes. Each line contains the genotypes at one SNP, and there are two entries per sample (222 total). Can you tell how the different alleles and missing genotypes are encoded?<br><br> 

`Admixture` reads a slightly different version of this format, so we first need to convert our files to the appropriate format. We will also take the oportunity to do some filtering. Alleles that ocurr at low frequencies or are genotyped for a small subset of the samples are known to cause problems with admixture proportion estimation, so we will retain only sites genotyped at >75% of individuals, with allele frequencies greater than 0.05. In addition, we will only keep SNPs that are at least 10,000bp away from each other, since `Admixture` and `eems` (which we will use later) assume unlinked sites. 

```bash
plink --tfile nacanids_111indiv_unrel_noYNP_42Ksnps --dog --missing-genotype N --biallelic-only --geno 0.25 --maf 0.05 --bp-space 10000 --make-bed --out wolves_0.25mis_thinned
```
Can you match the filters to the flags passed to `plink`? Out of the ~42,000 sites genotyped, how many were retained after filtering? <br><br>

Now that we have our file in the adequate forat, we can run `Admixture`. It takes two arguments, the fenotype file, and the number of hypothetical populations to which we will be assigning our samples (i.e. $k$). 

```bash
"$software_dir"/admixture_linux-1.3.0//admixture wolves_0.25mis_thinned.bed 2
```

The above code runs `Admixture` for $k=2$. Considering our dataset contains samples covering a very large extent of land and many possible barriers to gene flow (i.e. the Rocky Mountains, Hudson Bay, etc...), we may want to look at a wider range of values fo $k$. We can use a loop to estimate admixture proportions assuming $k=2-8$

```bash 
for i in {2..8}; do ../admixture_linux-1.3.0//admixture wolves_0.25mis_thinned.bed $i; done
```
