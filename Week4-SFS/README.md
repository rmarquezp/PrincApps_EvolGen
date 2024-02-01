Estimating and Using Site Frequency Spectra
===================================

As we discussed in lecture, the site frequency spectrum (SFS; also known as allele frequency spectrum) is a widely-used and very informative statistic to summarize genetic variation within a population. Today's practical will deal with estimating the SFS using different approaches, and using it to calculate other summary statistics, namely the various estimators of $\theta$ that we covered in class.
<bc><bc>

## Known genotypes and ancestral states

Let's begin with the easiest possible case: The file called "W4Matrix.txt" contains data simulated under tha standard neutral coalescent model $(\theta=5000)$ across a hypothetical chromosome of similar length to the human chromosome 1 (250,000,000bp long). We've also assumed that these sites have a very high rate of recombination between them, so that most (or all) segregating sites do not share a genealogy with other sites. Since our theoretical expectation for the SFS, $\eta_i=\frac{\theta}{i}$, is <i>across all possible genealogies</i>, then these data should closely resemble the expectation. Download this file to your computer. 
<br><br>
Since we know the ancestral states, calculating the SFS for this dataset is quite easy. We just need to count the number of derived alleles at each site, and tabulate our counts. We can do this in R:

```R
# Read in the file. You may need to add the path to your data matrix file. 
genotypes=read.table("W4Matrix.txt")

# Quickly count the number of rows and columns 
nrow(genotypes)
ncol(genotypes)
```
<b>Question 1:</b> How many individuals are in our matrix? How many sites?
<br><br>
On to calculating the SFS. Since we score derived alleles as 1 and ancestrals as 0, we can quickly count the number of derived allele copies at each site by summing over each site's genotypes. We can use the `colSums()` function, which gives the sum for each column of a matrix or data frame. We can then use the `table` function to generate counts for each mutation pattern (i.e. number of derived alleles) across sites. 

```R
counts=colSums(genotypes)
SFS=table(counts)

## Plot!
barplot(SFS, xlab="Derived Allele count", ylab="Number of Sites", col="black")
```

How does the SFS look? Does it appear to follow the theoretical expectation? We can assess this visaully by plotting our observed and expected SFS alongside each other. The first step is generating our expected SFS.

<b>Question 2:</b> How would you obtain the expected SFS? Recall $\theta=5,000$.
<!---
<details>

<summary> Click here to see the answer </summary>
```R
expSFS=5000/1:19
```
</details>
--->
Now lets plot our observed and expected SFS alongside each other. 
```R
expSFS="insert your code here"
ObsExp=matrix(c(expSFS,SFS), nrow=2, byrow=T)

## Plot. Real good match!
barplot(ObsExp, beside=T, col=c("grey15", "grey65"), xlab="Derived Allele count", ylab="Number of Sites", names=1:19)
legend("topright", c("Expected", "Observed"), fill=c("grey15", "grey65"))
```
Now lets estimate $\theta$ using the equations we discussed in class. Fu and Li's estimator is simply the numver of singletons:
$$\hat{\theta}_{FL}=\eta_1$$
```R
Theta_FL=SFS[1]

Theta_FL
[1] 5088
```

Pretty close! Now lets calculate Watterson's estimator, which is equal to the number of segregating sites controlled by sample size:

$$\hat{\theta}_{W} = \frac{S}{\sum\frac{1}{i}} = \frac{\sum \eta_i}{\sum\frac{1}{i}}$$ 

<b>Question 3:</b> Try to write this one out on your own. Provide your code and answer in your report. 

<details>
   <summary> Click here to see the answer (without code)</summary>
   $$\theta_W=5020.38$$

<!---
```R
n=20
Theta_W=sum(SFS)/sum(1/1:(n-1))
  
> Theta_W
[1] 281.8696
```
--->

</details>

Finally, lets calculate Tajima's estimator. Try to do this one on your own as well. It is a bit more complicated, since you'll have to use a loop, but it is worth giving it a try. 
```math
\hat{\theta}_{\pi}= \frac{\sum d_{i,j}}{n(n-1)\big{/}2}=\frac{\sum i(n-i)\eta_i}{n(n-1)\big{/}2}
```
<details>
   <summary> Click here to see the answer (with code)</summary>

```R
n=20
dist=c()

for(i in 1:(n-1)){dist[i]=i*(n-i)*SFS[i]}

Theta_pi=sum(dist)/((n*(n-1)/2))
  
Theta_pi
[1] 4969.742
```

</details>

So overall it seems like our Theta estimators are relatively close to the expected value. It is nice when simulations behave like the model that generated them :-).  

## Using Real Data
  
Now lets estimate teh SFS from some real data. We will again use the snowshoe hare sequences that we worked with on Week 3. In the interest of time, the reads for 6 snowshoe hares and two individuals from outgroup species (more on this below) have already been mapped back to the snowshoe hare genome. The procedure used for this was very similar to what we did on Week 3. Below is an example of the code used for one file, in case you find it useful.
<br><br>
<b>Download, trim, and QC reads</b>
  
  ``` bash
  
#Set variables for the number of threads, memory, and sample name
threads=8
mem=2G
sample=SRR6485258

#Download and dump reads. Note how we gave it more memory and threads so it goes faster
prefetch $sample
fasterq-dump -m $mem -e $threads $sample

## Clean reads
../software/skewer/skewer -o $sample -t $threads -q 30 -l 36 "$sample"_1.fastq "$sample"_2.fastq

#Always good manners to compress our files
gzip "$sample"-trimmed-pair1.fastq & gzip "$sample"-trimmed-pair2.fastq

## fastqc

fastqc --threads $threads "$sample"-trimmed-pair*.fastq.gz

``` 
<b>Map back to the genome</b> (if reads look nice and clean). 

```bash
#!/bin/bash

module load Bioinformatics bwa samtools

threads=8
sample=SRR6485258

#set variable for the reference genome
ref=/scratch/eeb401s002w24_class_root/eeb401s002w24_class/shared_data/ReferenceGenomes/GCF_033115175.1_mLepTim1.pri_genomic.fna

#Map. Note how instead of running each step separately we "pipe" the output of one command into the next using the | symbol.
#This saves time and storage space, since the computer doesn't need to write intermediate files. 

bwa mem -t $threads $ref "$sample"-trimmed-pair1.fastq.gz "$sample"-trimmed-pair2.fastq.gz | samtools fixmate -m - - | samtools sort -@ $threads - | samtools markdup -@ $threads -O bam - "$sample".sort.mkdup.bam

samtools index "$sample".sort.mkdup.bam
samtools flagstat "$sample".sort.mkdup.bam > $sample.flagstat.txt
```
   
All of our samples had >99% read mapping and relatively low (~5-10%) duplication levels. We are now ready to go. Sign into the Greatlakes cluster as we did last time, create a directory named `Week4`, move into it (`cd Week4`) and request a job. This time around we will need some more memory, and we'll again ask for 4 cores. 

```bash
salloc --account eeb401s002w24_class --time 1:30:00 --mem 48G --tasks-per-node 4
```
   
<b> What is the ancestral allele though?</b>
<br><br>
Before we start, you may have wondered: How will we know which allele is ancestral and which is derived? We can't travel in time and sequence the ancestor of our snowshoe hare population, so we'll need to make an educated guess. A common way to do this, is to use the allele present in a closely related species, often called an <i>outgroup</i> species. Of course, there is the possibility that a difference between our outgroup and focal species is due to mutation in the outgroup. To ameliorate this, we can use more than one outgroup. Hypothetically the two outgroups could have experienced the same mutation at the same site in parallel, but this is highly unlikely, so we can safely assume that if two different outgroup species have the same allele, it is most likely the ancestral allele for our focal species. 
<br><br>
Fortunately, there is a pretty good phylogeny of hares available, and exome-enrichment data has been generated for several species in the genus (<i>Lepus</i>). Based on the phylogeny below, generated by Ferreira et al (2020; https://doi.org/10.1093/sysbio/syaa088), I chose the Black-tailed jackrabbit (<i>L. californicus</i>) and the Mountain hare (<i>L. timidus</i>), for which exon enrichment data are also available from a different study. 

<img src="../Images/HareTree.png" width="500">
   
Reads from one individual of these species were mapped to the snowshoe hare genome as detailed above. I then used these read mappings to generate an "ancestral genome sequence", which looks just like the actual genome sequence, but has our best guess for the ancestral allele at each position. The proceudre to generate this is realtively straightworward, but in the itnerest of time we are going to skip it and go straight ahead to SFS estiamtion. 

Lets first load the necessary modules, and create some variables to make our code tidyer. The variables just include the long <i>paths</i> to the folders where many of our files of interest will be found. 

```bash
   module load Bioinformatics angsd R

## Set some directories as variables

ref=/scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/RefGenomes/LepAme_RefGenome_GCA_004026855.1.fa
lists=/scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/lists
bams=/scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/bams
```
Type `ls $bams`. What do you see in there?
<br><br>

OK, time to estimate our SFS frequencies. The first step is to go over all the genotyped sites and calculate the likelihood of all possible SFS at each site.

```bash
angsd -b "$lists"/ingroup.filelist -GL 1 -anc ancestral.fa -ref $ref -sites OutgroupFixedSites.txt -rf OutgroupFixed_chr.txt -dosaf 1 -out L_amer -baq 1 -C 50 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -doMajorMinor 1 -doHWE 1 -minHWEpval 0.001
```

Lets unpack this command a little. <br><br>

The `-b` flag tells Angsd to use the bam files in a list that can be found at "$lists"/ingroup.filelist<br>
   `-GL` specifies the specific equation to ge used for genotypel likelihood calculations. <br>
   `-anc` provides our "ancestral" reference geomne.<br>
   `-ref` provides the reference geomne to which we mapped our reads.<br>
   `-sites` points to a file that contains information on the sites that were identical between both outgroup species, and asks Angsd to onoy look at these sites. <br>
   `-rf` gives a file with the regions that include the "good" sites (see previous flag). This speeds up the process by letting Angsd know where to look for these sites. <br>
   `-dosaf 1` tells Angsd to calculate site allele frequencies and their associated likelihoods.<br>
   `-out` Specifies the name of the output files.<br>
   <br>
   After this there are multiple flags aimed at removing data that isn't in great shape, for example by setting minimum base quality thesholds, using only sites that don't deviate too much from HWE, or removing reads that don't map too well. 

<br><br>
   The command will take about 30 minutes to run. Since you may not want to wait this long, feel free to cancel the run (`ctrl+C`) and copy the output of this step into your directory: `cp /scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/L_amer.saf.* .`
   
   Now that we have our counts, we can try to find the most likely SFS. Angsd comes with an auxiliary program that maximizes the likelihood of the SFS, called realSFS. We can ask it to run for at most 200 iterations while trying to find the best SFS. 
   
   ```bash 
   realSFS -maxIter 200 L_amer.saf.idx > L_amer.ml.sfs.txt
   ```
Download the file called L_amer.ml.sfs.txt to the desktop as we did last time. You can read it into R using `sfs<-scan("L_amer.ml.sfs.txt"), and plot it as we did above. Note that this SFS has 13 entries, as it includes values for fixed sites both at the ancestral (first entry) and derived (last entry) allele. You may want to plot values 2-12. <br>
<br>
How does the SFS look? Does it seem real?
<br><br>
   <b>What if we can't confidently infer the ancestral allele</b><br><br>
   
   I'm sure you can imagine several situations in which infering ancestral bases confidently may be difficult or problematic. In those cases, we can use a version of the SFS often called the <i>folded</i> or <i>minor</i> allele frequency spectrum. It is still a frequency spectrum, but uses the frequency of the minor (ie least comon) allele, regardless of which one is ancestral. Therefore, for n samples it goes from 1 to (n/2), instead of 1 to n. We can obtain this spectrum from the unfolded (ie. regular) SFS by just adding $\eta_i+\eta_{n-i}$. For example, for 10 samples, sites with 1 derived alleles get added to sites with 9 derived alleles, since in both cases the least common allele is at frequency 0.1. <br><br>
   
   We can easily produce a folded spectrum on reafSFS by just adding the `-fold 1` flag. 
   
   ```bash
   
  realSFS -maxIter 200 -fold 1 L_amer.saf.idx > L_amer.ml.folded.sfs.txt
```
   Does this look better?
