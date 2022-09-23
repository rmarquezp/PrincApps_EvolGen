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
  
  ```R
  Theta_W=sum(SFS)/sum(1/1:19)
  
> Theta_W
[1] 281.8696
```

