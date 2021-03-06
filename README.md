<h1 align="center"> DNA Methylation Data Analysis </h1>


<img src = "https://user-images.githubusercontent.com/88287648/131827248-562cf188-40c2-48a6-82c9-75e683cef24d.jpeg">


##  π The Goal of Our Analysis

Oncogenic transformation of normal cells often involves epigenetic alterations, such as histone modification and DNA methylation. For this analysis, we obtained a subset of the whole genome bisulfite sequencing (WGBS) data to analyse the pattern of DNA methylomes of normal breast and cancerous breast cells. 
This workflow is a complete reproduction of the workflow mentioned in the Galaxy training material on Methyl-Seq Analysis.

## π Our Workflow
 <p align = "center"> 
  <kbd> <img src ="https://user-images.githubusercontent.com/88287648/131124866-4f63e72a-6524-418f-8cf6-5786e00de8ac.png" </kbd>
</p>
  
## π Dataset and its source

The datasets we used are mainly a subset of the data from the publication "Hierarchical Clustering     of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes" (10.1371/journal.pone.0118453). 

## **Lets Move to the Analysis π**

### Get Data

````python
  wget https://zenodo.org/record/557099/files/subset_1.fastq?download=1
  ````
````python
  wget https://zenodo.org/record/557099/files/subset_2.fastq?download=1
  ````

### Quality Control

- Download the FastQC module
Note: FASTQC requires java and javac installed for implementation and you need to run the fastqc file from the folder (using the relative/absolute links to the sequence reads)<bR> 
```python
  sudo apt install default-jre
  ```
```python
  sudo apt install default-jdk
  ```
  Make the βfastqcβ an executable file<bR>
```python 
  chmod 755 fastqc
  ```
- Run the fastqc on all sequenced reads from its folder<bR> 
```python 
  fastqc subset_1.fastq?download=1  subset_2.fastq?download=1
  ```
  
If you open the html files you will see the result as follows-

<p align="center"> <kbd> <img src="https://user-images.githubusercontent.com/88287648/131125660-6a6a6eed-eec0-44b0-bd74-9bd76bff8062.png" width = "800" height = "600"> </kbd> </p> 
    
This tells that the per base sequence quality was quite good for the datasets. The mean quality of reads (indicated by the blow lines) showed that even the lowest score is above 28.


<p align="center">  <kbd> <img src="https://user-images.githubusercontent.com/88287648/131812509-8debe5ec-c545-48aa-ba8f-d97688c3c93f.png" width = "800" height = "600"></kbd> </p> 

The βper base sequence contentβ segment from the fastq reports shows a drop of βCβ and a rise of βTβ bases. Itβs because- Every C-meth stays a C and every normal C becomes a T during the bisulfite conversion.

### Mapping
Mapping of bisulfite-sequencing reads needs different aligners from the normal NGS sequencing read aligners. As in a BS-seq reads, all the Cβs are C-methβs and a T can be a T or a C, the mapper for methylation data needs to find out what is what. Bismark is such a suitable aligner for BS-seq reads.
Due to the memory shortage of our device, we will use only Chromosome_1 as our reference genome for aligning reads.
  
Download reference genome  

```python 
  wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
  ```  
Download Perl (a prerequisite for Bismark) 
```python 
  conda install -c bioconda perl-lwp-simple
  ```  
Download Bismark 
```python 
  conda install -c bioconda bismark
  ```  
Create a folder for moving the chr1 reference genome into that 
```python 
  mkdir genome_folder
  ```  
 Now move the chr1 reference genome file 
```python 
  mv Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz genome_folder/
  ```  
 Create indexing of chromosome 1 using bismark 
```python 
  bismark_genome_preparation --verbose genome_folder/
  ```  	
 Time for mapping 
 ```python 
  Bismark -genome genome_folder/ -1 subset_1.fastq?download=1 -2 subset_2.fastq?download=1
  ```  
  ```python 
  bismark_methylation_extractor -p -- no_overlap file_1.fastq_bismark_pe.bam
  ```  		
Bismark generates a .bam file as its final output. The alignment stats from the Bismark tool is as follows. As we used only Chr1 as our reference genome, it is showing the largest part of the reads for no alignment.

<h3 align = "center "> 
<kbd><img src="https://user-images.githubusercontent.com/88287648/131815589-b059bd4e-231f-48ca-831c-1fab893ea69c.png" width = "400" height = "500"></kbd>

<kbd><img src = "https://user-images.githubusercontent.com/88287648/131817223-da72afcf-f061-49af-9249-a77b986a8fe8.png" width = "700" height = "500"></kbd>

### Generating a sorted bam file with Samtools

MethylDackel needs a sorted bam file for its execution. For that, Samtools is used to generate a sorted bam file from the output bam file of aligner.

Install SamTools  
```python 
  conda install -c bioconda samtools
  ```  
Sort the output file from alignment with Samtools 
```python 
  samtools sort subset_1_bismark_bt2_pe.bam -o bam_sorted_by_samtools.bam
  ```  
  
### Methylation Bias & Metric Extraction
  
Methylation bias plot is helpful for looking at the distribution of methylation and searching for any possible bias.

Install MethylDackel   
```python 
  conda install -c bioconda methyldackel
  ```  
Now Plot the Methylation Bias 
```python 
  MethylDackel mbias Chr1_ref.fa.gz bam_sorted_by_samtools.bam methylation_bias_by_methylDackel
  ```  
This plot shows distribution of the CpG methylation here is a bit biased along the genome which is not expected. If we were to trim the reads, we would include roughly the positions 0 to 134, for the both strands just to exclude the most biased part from the 3' end. The difference of DNA methylation pattern between the two datasets are clearly visible.

<p align="center"> <kbd> <img src="https://user-images.githubusercontent.com/88287648/131824142-c91f568e-9c40-4071-b4dc-ce6421a72c7b.png"> </kbd>
  
To extract the methylation on the resulting BAM file of the alignment step 
```python 
  MethylDackel extract --mergeContext Chr1_ref.fa bam_sorted_by_samtools.bam
  ```  
### Visualization
  
BedGraph-to-bigWig was used to convert the Bedgraph file containing methylation level to a bigwig file

Install bedGraphtoBigWig 
```python 
  conda install -c bioconda ucsc-bedgraphtobigwig
  ```  
Now run the following commands-
```python 
  awk 'NR!=1' bam_sorted_by_samtools_CpG.bedGraph > out_from_bigwig.deheader.bedGraph
  ```  
```python 
  sort -k1,1 -k2,2n out_from_bigwig.deheader.bedGraph
  ```  
```python 
  sort -k1,1 -k2,2n out_from_bigwig.deheader.bedGraph > sorted.bedGraph
  ```  
Create a file named Chr1.chrom.sizes and write the size of Chromosome 1 in that in a tab delimited style
```python 
  nano Chr.chrom.sizes
  ```  
Copy chr1	249250621  to the file, save and exit. And then run-
  
```python 
  awk '{print $1,$2,$3,$4}' sorted.bedGraph > final_sorted_with4.bedGraph
  ```  
```python 
  bedGraphToBigWig final_sorted_with4.bedGraph Chr1.chrom.sizes final_output_from_bigwig.bw
  ```  
Deeptools contains several helpful tools for helpful visualization, including ComputeMatrix and plotProfile. With CpG islands data as regions to plot, ComputeMatrix helps to convert the bigwig file into a matrix file.

```python 
  conda install -c bioconda deeptools
  ```  
```python 
  computeMatrix scale-regions -S final_output_from_bigwig.bw -R CpGIslands.bed -b 1000 -out output_from_ComputeMatrix
  ```  
Finally, the methylation level of our data around the transcription start site (TSS) was plotted with plotProfile.
```python 
  plotProfile --matrixFile output_from_ComputeMatrix --outFileName output_from_plotProfile
  ```  
<p align="center"><kbd> <img src="https://user-images.githubusercontent.com/88287648/131825989-3e7ce08c-40b0-4a15-955a-391b842af718.png"></kbd>
  
This plot shows the methylation level around all Transcription Start Sites of chrosome_1. When located at gene promoters, DNA methylation is usually a repressive mark.

### References 
1. Lin, I.-H., D.-T. Chen, Y.-F. Chang, Y.-L. Lee, C.-H. Su et al., 2015 Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes (O. El-Maarri, Ed.). PLOS ONE 10: e0118453.
2. RamΓ­rez, Fidel, Devon P. Ryan, BjΓΆrn GrΓΌning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike DΓΌndar, and Thomas Manke. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.     
3. Joachim Wolff, Devon Ryan, 2019 DNA Methylation data analysis (Galaxy Training Materials). https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html Online; accessed Mon Aug 23 2021
4. Batut et al., 2018 Community-Driven Data Analysis Training for Biology Cell Systems 10.1016/j.cels.2018.05.012
5. Krueger, Felix, and Simon R. Andrews. βBismark: A Flexible Aligner and Methylation Caller for Bisulfite-Seq Applications.β Bioinformatics, vol. 27, no. 11, 2011, pp. 1571β1572., doi:10.1093/bioinformatics/btr167
6. Devon Ryan, 2021 dpryan79/MethylDackel: A (MOSTLY) UNIVERSAL METHYLATION extractor For BS-seq experiments, https://github.com/dpryan79/MethylDackel
  
### Acknowledgements
We acknowledge the HackBio Team for the platform created to bring us all together learning and working as a team progressively on a bioinformatics track.
