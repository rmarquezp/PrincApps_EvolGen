Handling Illumina data on the Greatlakes HPC Cluster
====================================================

This week we will begin working with data. As we discussed in lecture, one of the most comonly used genotyping strategies in current population genetic studies consists on using the Illumina technology to sequence short fragments of DNA, which we will call <b>reads</b>. These short reads are them aligned against a previously generated reference genome, which allows us to 1. find their location in the genome, and 2. obtain genotypic information for our individual(s) of interest. 
<br><br>In this practical we will learn how to:<br>
* Interact with a remote high performance computing (HPC) cluster through the terminal.
* Download data from the NCBI's [Short Read Archive](http://www.ncbi.nlm.nih.gov/sra) (SRA), a public repository of sequence data.
* Conduct quality-control analyses on Illumina data.
* Align these data against a reference genome.
* Asses the quality of mapped data.
<br><br>
## Study System

For this and several other practicals we will be using data from a sample of snowshow hares (<i>Lepus americanus</i>), which was collected by Jones et al. ([2018](https://doi.org/10.1126/science.aar5273)). Briefly, the authors collected tissue samples, and generated sequence data using a technique called [exome enrichment](https://en.wikipedia.org/wiki/Exome_sequencing), in which only the coding sequences in the genome are sequenced. This is often done to reduce costw and computational effort when sequencing the entire genome is not necessary. 

<img src="https://www.nrcm.org/wp-content/uploads/2021/12/snowshoe-hare2-bcomeau.jpg" width="600">

## Interacting with the Greatlakes Cluster
We usually don't interact with computer clusters such as UM's Greatlakes the way we do with personal computers. Instead of using the mouse and keyboard to give the computer instructions and receiving results on a screen, all directly connected with the computer, we will be remotely logging into the cluster and giinteracting with it through the <b>command line</b>. That means we will be typing commands into a prompt, and receiving results as text on that same window, similar to how we have interacted with R in previous sessions. To start using the command line open the Terminal (on Mac this is found under Applications -> Utilities). A window like this should open:

<img src="https://cdn2.macpaw.com/images/content/Screen%20Shot%202021-09-03%20at%2014.32.58_1630671309.png" width=500>
<br>
This window is where we will type commands and receive outputs from Greatlakes. The first step to do so is logging into Greatlakes. To do so we must use a command called "ssh" and our UM credentials.

```bash
ssh uniqname@greatlakes.arc-ts.umich.edu
```
You will then be asked for your password (this is the same password you use to log into other UM services, such as your email). Type it it and hit enter. <b>Don't worry if the cursor doesn't move as you type, this is normal</b>. If you entered your password correctly, you should connect to Greatlakes. Before letting you in, the cluster will ask for Duo two-factor authentication, giving you multiple options to do so. 

```bash
Enter a passcode or select one of the following options:

 1. Duo Push to XXX-XXX-NNNN
 2. Phone call to XXX-XXX-NNNN
 3. SMS passcodes to XXX-XXX-NNNN
 ```
 
Choose your preferred option and type its number (or enter a passcode), and hit enter. After authenticating you should see a welcome screen, and a command prompt ready for you to type. On the terminal we can do 

## Obtaining Fastq files
The first step in most bioinformatic pipelines is transferring the data to our work environment. If you have generated these data yourself this may involve transferring it from the sequencing facility's computer to yours. If you are using data available in 
