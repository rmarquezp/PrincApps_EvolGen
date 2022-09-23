Estimating Site Frequency Spectra
===================================

As we discussed in lecture, the site frequency spectrum (SFS; also known as allele frequency spectrum) is a widely-used and very informative statistic to summarize genetic variation within a population. Today's practical will deal with estimating the SFS using different approaches, and using it to calculate other summary statistics, namely the various estimators of $\theta$ that we covered in class.
<bc><bc>

## Known genotypes and ancestral states

Let's begin with the easiest possible case: The file called "W4Matrix.txt" contains data simulated under tha standard neutral coalescent model ($\theta=300$) at 1,000 sites across a hypothetical genome. We've also assumed that these sites are very far from one another along the genome, to the point where they recombine freely. This means each site is part of a locus with a different genealogy. Since our theoretical expectation for the SFS, $\eta_i=\frac{\theta}{i}$, is <i>across all possible genealogies</i>, then these data should closely resemble the expectation. Download this file to your computer. 
<br><br>
Since we know the ancestral states, calculating the SFS for this dataset is quite easy. We just need to count the number of derived alleles at each site, and tabulate our counts. We can do this in R:

```R
# Read in the file. You may need to add the path to your data matrix file. 
genotypes=read.table("W4Matrix.txt")

# Generate counts. 
#Since we score derived alleles as 1 and ancestrals as 0, we can quickly count the derived alleles by adding all the genotypes at each site.
#The colSums() function gives the sum for each column of the matrix
counts=colSums(genotypes)

#Use the table() function to tabulate our results
SFS=table(counts)

## Plot!
barplot(SFS, xlab="Derived Allele count", ylab="Number of Sites", col="black")
```

How does the SFS look? Does it appear to follow the theoretical expectation? We can assess this visaully by plotting them alongside each other. 

```R
# First create a matrix with observed and expected SFS together. Recall Theta=300
ObsExp=matrix(SFS, 300/1:20, ncol=2, byrow=F)

## Plot
barplot(ObsExp, beside=T, col=c("grey15", "grey65"), xlab="Derived Allele count", ylab="Number of Sites", names=1:19)
legend("topright", c("Expected", "Observed"), fill=c("grey15", "grey65"))
```
Now lets estimate $\theta$ using the equations we discussed in class. Fu and Li's estimator is simply the numver of singletons:
$$\hat{\theta}_{FL}=\eta_1$$
```R
Theta_FL=SFS[1]

> ThetaFL
[1] 309
```

Pretty close, now lets calculate Watterson's estimator, which is equal to the number of segregating sites controlled by sample size:

$$\hat{\theta}_{W} = \frac{S}{\sum\frac{1}{i}} = \frac{\sum \eta_i}{\sum\frac{1}{i}}$$ 

Try to write this one out on your own. What answer do you get?

<details>
   <summary> Click here to see the answer</summary>

```R
n=20
Theta_W=sum(SFS)/sum(1/1:(n-1))
  
> Theta_W
[1] 281.8696
```

</details>
  
Finally, lets calculate Tajima's estimator. Try to do this one on your own as well. It is a bit more complicated, since you'll have to use a loop, but it is worth giving it a try. 

$$\hat{\theta}_{\pi}= \frac{\sum d_{i,j}}{n(n-1)\big{/}2}=\frac{\sum i(n-i)\eta_i}{n(n-1)\big{/}2}$$
  
<details>
   <summary> Click here to see the answer</summary>

```R
n=20
dist=c()

for(i in 1:(n-1)){dist[i]=i*(n-i)*SFS[i]}

Theta_pi=sum(dist)/((n*(n-1)/2))
  
> Theta_pi
[1] 277.1421
```

</details>

So overall it seems like our Theta estimators are relatively close to the expected value. 

## Using Real Data
  
Now lets estimate teh SFS from some real data. We will again use the snowshoe hare sequences that we worked with on Week 3. In the interest of time, the reads for 6 snowshoe hares and two individuals from outgroup species (more on this below) have already been mapped back to the snowshoe hare genome. The procedure used for this was very similar to what we did on Week 3, with the only difference being that duplicated reads were marked after mapping, so that the programs that we will use in downstream analyses know to treat duplicates as such. Below is an example of the code used for one file:
  
  ``` bash
  
module load Bioinformatics sratoolkit bwa samtools picard-tools
  
## Download file
  
prefetch SRR11020246; fasterq-dump SRR11020246
  
## Quality-trim the reads
  
skewer -o SRR11020246 -t 4 SRR11020246_1.fastq SRR11020246_2.fastq
  
## Map and sort. Note that we are doing both steps in a single line of code. 
ref=/scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/RefGenomes/LepAme_RefGenome_GCA_004026855.1.fa

bwa mem -t 8 $ref SRR11020246-trimmed-pair1.fastq SRR11020246-trimmed-pair2.fastq | samtools sort -O BAM -@  4 -o SRR11020246.sorted.bam

  
## Add some info to the file prior to duplicate marking 
PicardCommandLine AddOrReplaceReadGroups I=SRR11020246.sorted.bam O=SRR11020246.sorted.RG.bam SORT_ORDER=coordinate RGLB="$ref" RGPU=NONE RGSM=SRR11020246 RGID=SRR11020246 RGPL=Illumina VALIDATION_STRINGENCY=LENIENT


# Mark duplicates
PicardCommandLine MarkDuplicates I=SRR11020246.sorted.RG.bam O=SRR11020246.sorted.deduped.bam M=SRR11020246.dupMetrics.txt AS=true CREATE_INDEX=true MAX_FILE_HANDLES=1000 VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false

# Remove intermediate files
rm SRR11020246.sorted.RG.bam
rm SRR11020246.sorted.bam
```

Before we start, you may have wondered: How will we know which allele is ancestral and which is derived? We can't travel in time and sequence the ancestor of our snowshoe hare population, so we'll need to make an educated guess. A common way to do this, is to use the allele present in a closely related species, often called an <i>outgroup</i> species. of course, there is the possibility that a difference between our outgroup and focal species is due to mutation in the outgroup. To ameliorate this, we can use more than one outgroup. Hypothetically the two outgroups could have experienced the same mutation at the same site in parallel, but this is highly unlikely, so we can safely assume that if two different outgroup species have the same allele, it is most likely the ancestral allele for our focal species. 
<br><br>
Fortunately, there is a pretty good phylogeny of hares available, and exome-enrichment data has been generated for several species in the genus (<i>Lepus</i>). Based on the phylogeny below, generated by Ferreira et al (2020; https://doi.org/10.1093/sysbio/syaa088), I chose the Black-tailed jackrabbit (<i>L. californicus</i>) and the Mountain hare (

  
  
  
