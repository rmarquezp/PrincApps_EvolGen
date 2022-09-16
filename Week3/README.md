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
 
Choose your preferred option and type its number (or enter a passcode), and hit enter. After authenticating you should see a welcome screen, and a command prompt ready for you to type. On the terminal we can run analyses, but also perform routine tasks such as moving between directories or creating new folders.files. Before we start analyzing data, lets create some directories for our work. We can use the "mkdir" command for this. 
```bash
 mkdir Week3
 ```
 This will create a directory called "Week3". To make sure this worked, we can use the "ls" command to list all the files and directories at our current location
 ```bash
 ls
 
 Week3
```
Running this command confirms we've succesfully created our folder. Lets now go to this folder, where we can run our commands. 

```bash
cd Week3
```
Until now, we've been working on the <i>head node</i> of the cluster. This is a computer meant for logging in and running menial tasks, such as moving/copying files and creating new directories. To run computationally intensive tasks, 
we use <i>compute nodes</i>, which are more powerful, and exclusively allocatod for this purpose. To gain access to a compute node, we can use the "srun" command:
```bash 
srun --account eeb401s002f22_class --time 1:30:00 --mem 8G --tasks-per-node 1 --pty bash
```
This command asks for acces to a compute node with 8Gb RAM and one processor for 1.5 hours. The resources used will come from our class allocation (`--account eeb401s002f22_class`). A few momments after typing this you should get a message saying the requested resources have been allocated, together with a command prompt in which you can type. This is where we will work today. Before we start, we need to load some <i>modules</i>, which contain the programs that we will use. This is analogous to loading packages in R. 
```bash
module load Bioinformatics bwa sratoolkit samtools fastqc
```

## Obtaining Fastq files
The first step in most bioinformatic pipelines is transferring the data to our work environment. If you have generated these data yourself this may involve transferring it from the sequencing facility's computer to yours. If you are using data available in the SRA, we can use a set of programs specifically created by the NIH to interact with the SRA, called the SRA toolkit. Each file uploaded to the SRA is assigned a unique ID number. Below are the ID numbers for 12 read files from our snowshoe hare system. 

```bash
1. SRS6102628
2. SRS6102626
3. SRS6102627
4. SRS6102625
5. SRS6102623
6. SRS6102624
7. SRS6102622
8. SRS6102621
9. SRS6102620
10. SRS6102619
11. SRS6102618
12. SRS6102617
```
Pick one of these files so that each person in the class analyzes a different file. Nowlets download it using the SRA Toolkit. We will first download each file and then extract it in <i>fastq</i> format. This format is a commonly used file format for sequence data, which includes base calls and associated quality values (ie. error rates) for each base. 

``` bash
#Download data, replace the file ID for your chosen file. 
prefetch fileID

#Extract it as fastq
fasterq-dump fileID
```
These sequences were generated using a technology called<i>paired-ejnd</i> sequencing where the beginind and end of each DNA fragment is sequenced. Therefore, after running the commands above, you should find two fastq files, one containing the beginings, and another the eands of the sequenced DNA fragments. Use the "ls" command to verify they are there. Lets now take a look inside these files 
```bash
less fileID.fastq
```
Note how each sequence is represented by four lines:
```bash
@SRR11020240.1 1 length=101
ATCGATGATTAAATCACCCTAATTTGCATTGTCTGAGCTAATCACCGATGATTTATTTACCTGCTATGTTTACACGAAGTGGATAGCAACGATTAAGTTTA
+SRR11020240.1 1 length=101
AAFFFJJJFJJJJJJJJJJFJJJJFJJJJJJJJJJJJJJJJJJJJJJJJAJ-FJJJJJJJFJJJFFJF7FFJJJJFJFFJJJJFJFJJJFJJA7JJFF<FF
```
The first and third lines include information avout the read, in this case the file name, read number, and read length. The second line contains the sequence, and the fourth line contains the quality for each base. To save space, numerical values are stored as letters istead of numbers. YOu can read more about the fastq file [here](https://en.wikipedia.org/wiki/FASTQ_format). 

## Assessing read quality and trimming reads

As we covered ijn class, Illumina sequences have some probability of including technical artifacts, mostly in the form of sequencing errors and contamination introduced by the library preparation process. Therefore, it is advisable to evaluate the quality of our reads and to process them in order to remove low-quality and contaminated segments of reads. 
<br><br>
To assess read quality we can use aprogram called [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
```bash
fastqc fileID_1.fastq fileID_2.fastq
```
This will produce two files called fileID_1_fastqc.html and fileID_2_fastqc.html. These files can be visualized in an internet browser. To do so they need to be downloaded to your local machine. Open a new terminal window by clicking on "Shell -> New Window" on the top menu bar. In the new terminal window type:
```bash
scp "marquezr@greatlakes.arc-ts.umich.edu:Week3/*.html" ~/Desktop
```
This command copies all files with the extension ".html" in the folder Week3 to your Desktop. Wait here while everyone catches up, and we will inspect this file as a group. 
<br><br>
Overall our reads seem to be OK, but there seems to be some adapter contamination. There are many many programs to deal with this sort of artifact. One that I particularly like is called [skewer](https://github.com/relipmoc/skewer). Run it as follows:

```bash
/home/marquezr/software/skewer/skewer -o fileID fileID_1.fastq fileID_2.fastq
```
This should create two new fastq files. Use the "ls" command to see their names, and run fastqc on tehse new files. This should generate two new html files. Download them to the desktop and check them out. Was the adapter contamination removed?

## Read Mapping

Now that we have high-quality, uncontaminated read files, we can map them back to the reference genome, using the program "bwa". The reference genome file is located in `/scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/RefGenomes/LepAme_RefGenome_GCA_004026855.1.fa`.

```bash
#Index the genome file. This only needs to be done once. It takes a few hours, so it has been done for you. 

# bwa index /scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/RefGenomes/LepAme_RefGenome_GCA_004026855.1.fa

## Now map the reads back to the genome

bwa mem -o fileID.sam /scratch/eeb401s002f22_class_root/eeb401s002f22_class/shared_data/RefGenomes/LepAme_RefGenome_GCA_004026855.1.fa fileID-trimmed-pair1.fastq fileID-trimmed-pair2.fastq
```
The mapped reads in our output file are organized in the order in which there were sequenced. We can save a considerable amount of space by 1. organizing them according to the regions of the genome that they map back to, and 2. compressing the file. We can achieve this using the program "samtools":

```bash
samtools sort -O BAM -o fileID.sorted.bam SRR11020240.sam

#We can remove the unsorted filel

rm -r fileID.sorted.bam
```
Now lets have a look at what is inside our file:

```bash
samtools view
```

NEXT DO SAMTOOLS FLAGSTAT AND SAMTOOLS DEPTH, PLOT DEPTH HISTOGRAM. IF THERE IS TIME DOWNLOAD TO IGV. 
