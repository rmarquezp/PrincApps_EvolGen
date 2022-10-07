Calculating Linkage Disequilibrium
==============================
Linkage disequilibrium is a very useful source of information in population genetics, which can be used to make inferences about demography, selection, and genome architecture, among others. Today, we will be exploring two ways of visualizing linkage disequilibrium in the genome: its decay with physical distance in the genome, and the length and distribution of runs of homozygocity (ROHs) along the genome. 
## Study System
In today's practical, we will focus on populations of two closely-related species of <i>Drosophila</i> flies that inhavit the Seychelles archipelago in West Africa, <i> D. sechellia</i> and <i> D. simulans</i>. The former is endemic to some of the islands in the Seychelles, while the latter is a human commensal with a worldwide distribution. We will be comparing the degrees of LD in populations of both species collected at a single locality in the Seychelles where they coexist. A big thanks to Daniel Matute and Adam Stuckert for making these data available prior to publication. 
## The data
Today we will be using data in the form of called genotypes. That is, reads have already been mapped to the genome, genotype likelihoods have been calculated, and a particular genotype at each site for each individual has been determined based on the GLs. We will therefore be assuming that these genotypes are true, and will not incorporate the unceratainty associated to this process. The data are contained in a file in variant call format, also known as a VCF. This format is basically a very large table where, sites are rows and individuals are columns. For each site, there is inforamtion on genotype, genotype likelihood, depth, and other associated quality measures for each individual. 


Therefore, before we begin, it is a good idea to make sure we have good enough depth
