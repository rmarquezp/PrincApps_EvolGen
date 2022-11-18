Visualizing Genetic Structure in Space
============

Visualizing the distribution of genetic variation across geography is a powerful tool in studies of biogeogprahy and historical demography. Today we will be using a few different approaches to visualize patterns of both discrete and continuous population structure on a map, which vary in the degree to which they explicitly consider geography in their estimation of demographic parameters. The first approach will consist of estimating admixture proportions and plotting them directly on the map, while the second will consist of directly estimating rates of migration across the landscape based off of the degree of genetic differentiation and sampling positions of a set of samples. 

## The Data

For today's practical we will be analyzing data obtained from gray wolves across Canada and Alaska by Schweizer et al ([2016](https://onlinelibrary.wiley.com/doi/10.1111/mec.13364). The authors used a SNP array (i.e. chip) designed for dogs to genotype 111 wolves at roughly 42,000 SNPs.<br><br>
<img src="https://nywolf.org/wp-content/uploads/2022/01/Trumpet_Lighthawk_family_winter_2020-scaled-2-1400x788.jpg" width="600">
<br><br>
Before we start, you may need to install some of the following R packages: `conStruct, rcpp, RcppEigen, raster, rgeos, sp, rworldmap`. Please take a minute to install them if needed.<br><br>

Lets now explore our dataset a little. To begin, plot the collection localities on a map of North America. For this we will need the file called `wolves.coord`, which contains the coordinate for each one of our samples. 

```R
library(raster)
library(rworldmap)


altitude=getData('worldclim', var='alt', res=2.5)
map=getMap(resolution="low")

northNA=extent(-175, -54, 43, 85)
altitude_crop=crop(altitude,northNA)
map_crop=crop(map,northNA) ##May throw a warning but seems to work fine regardless

coords=read.table("wolves.coord")

plot(altitude_crop, xlim=northNA[1:2], ylim=northNA[3:4])
plot(map_crop, lwd=1, add=T, xlim=c(-175,-54), ylim=c(43,80))
points(coords, pch=21, bg="white")
```
From the map we can see that the sampling was pretty consistent across northern North America, which means we will probably not have artifacts associated to uneven sampling. 

<img src="../Images/wolfMap.png" width="600">
