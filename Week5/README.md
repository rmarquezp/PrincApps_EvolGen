Demographic Inference using Coalescence
==========

The evolutionary history of a lineage is often marked by changes in population size. Therefore, estimating historical demographic changes can be a highly informative approach for evolutionary biologists. Coalescent theory provides a powerful framework for this, since the rate of coalescence depends on the (effective) population size. Today we will be impleenting two approaches for demographic inference under the coalescent model: The "skyline" plot and demographic parameter estimation from the site frequency spectrum (SFS). For today's lab you will some software installed: <a href="https://www.beast2.org/" >BEAST2</a>, a very powerful suite of programs for the inference and analysis of evolutionary trees, <a href="https://github.com/beast-dev/tracer/releases/tag/v1.7.2" >Tracer</a>, used to visualize the results of Bayesian MCMC analysis, and <a href="http://cmpg.unibe.ch/software/fastsimcoal27/">fastsimcoal2<a/>. If you haven't please install the first two (locally) on your cumputer. The third one is already on greatlakes. 

## Skyline plot

In cases where we can confidently estimate the genealogy of a locus <i>that does not recombine</i>, we can use the coalescence times in this genealogy to estimate the population size at different intervals of our focal lineage's history. In many cases finding non-recombining loci longer than a few base-pairs is very difficult, since recombination is pretty ubiquitous across genomes, and estimating the genealogy of a short sequence is challenging, since it provides very little data to work with. However some types of genetic material, such as mitochondrial/chloroplast DNA or viral DNA/RNA sequences, do not recombine, and have relatively high mutation rates, so we can easily and confidently estimate their genealogies. For today's lab, we will be estimating recent changes in the population size of the SARS-Cov-2 virus, which causes covid 19, in the state of Michigan. 
<br><br>
The National Center of Biotechnology Information (NCBI) maintains a database with a large number of the SARS-Cov-2 sequenced by genetic surveillance efforts. Today we will be analyzing the 911 sequences in this database that come from samples collected from patients in Michigan this month. The file `MI_Covid.aln.fasta` contains an alignment of all of these sequences. Ideally, we would analyze them all together. However, in the interest of time, we will be analyzing 100 randomly chose sequences out of this alignment. Since each of us will be analyzing a different random subset of 100 sequences, we'll be able to evaluate how much our results vary by our smaller sample size. 
  <br><br>
  Lets beging by sub-sampling our alignment. For this we will use the `ape` package in R. 
  ```R
  library(ape)

#Read in alignment. Remember to adjust your path!
aln=read.dna("MI_Covid.aln.fasta", format="fasta")

#Create a new alignment that only includes 100 random rows of the complete alignment
random=aln[sample(1:nrow(aln), 100, rep=F),]

#Write out the new file in nexus format. Remember to adjust your path!
write.nexus.data(random, file="MI_Covid.aln.100random.nex", format="DNA", interleaved=F)
```
  <b>Setting up the analysis in BEAUTi</b><br>
  
Now we'll set up a configuration file for BEAST to run our skyline plot analysis. This can be done using the program BEAUTi, which is distributed with BEAST. When you open BEAUTi, you'll find multiple different tabs at the top of the window. In the "Partitions" tab we will load our alignment. To do so go to `File -> Import Alignment`, and select your file with 100 random sequences. Once the file loads, you should see a new partition appear, with 100 taxa and 29,863 nucleotide positions. 
<br><br>
The "Tip Dates" tab can be used to give BEAST the dates at which out sequences were sampled, which can help callibrate the branch lengths from units of genetic distance to time. Since all of our sequences were collected this month, we will use an experimentally-determinaed mutation rate instead (see below), so we can skip the "Tip Dates" tab . 
<br><br>
In the "Site Model" tab we will specify our mutation model (also known as a model of DNA evolution). Since we are dealing with a relatively short DNA sequence with a high mutation rate, we will use a finitely-many alleles model. This is commonly done for tree-inference, since we need sequences that evolve fast enough to generate enough information that allows us to infer the tree, so back-mutation must be taken into account. We will use a model commonly known as "HKY", which uses the following rate matrix:
  
$$Q = \begin{bmatrix}
\cdot&  \pi_{C}& k\pi_{G}  & \pi_{T}  \\
\pi_{A}&  \cdot& \pi_{G}  & \pi_{T}  \\
k\pi_{A}&  \pi_{GC}&\cdot  & \pi_{T}  \\
\pi_{A}&  k\pi_{C}& \pi_{G}  & \cdot  \\
\end{bmatrix}
$$
