# LBD_CWOW
Lewy Body Dementia Center With Out Walls (LBD CWOW) RNAseq processing and analysis.
Anterior cingulate cortex tissue samples from the Mayo Clinic brain bank were collected for 619 individuals. 
Raw fastq files are available on SRA, PRJNA1023207.

| Disease                   | Count   |
| ------------------------- |:-------:|
| Control                   | 86      |
| Pathological  aging (PA)  | 39      |
| Alzheimer’s disease (AD)  | 54      |
| Lewy body dementia (LBD)  | 440     |

This git repo contains scripts for the following:
-   Metadata analysis of the 619 individuals. 
-   Processing of bulk-tissue RNA-sequencing data.
-   Gene-level differential expression among pairwise groups and within each genetic sex.
-   Weighted gene co-expression network analysis (WGCNA).  
-   Generation of manuscript figures from Olney et al. 202x publication 
-   Generation of shiny app for exploration of the results presented in Olney et al. 202x publication, view app [here](https://fryerlab.shinyapps.io/LBD_CWOW/)


## Create conda environment
The necessary software for bulk RNAseq data processing is contained in: `LBD.yml`.

To create the environment:
```
conda env create -n LBD --file LBD.yml

# To activate this environment, use
#
#     $ conda activate LBD
#
# To deactivate an active environment, use
#
#     $ conda deactivate

```
After the conda LBD environment has been created, you will need to additionally install GATK.
This step must be done manually and not through conda; see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4). The link will explain how to download GATK4, then you will need to add an alias to your bash profile. 
```
alias gatk='/path/to/gatk-package/gatk'
```
## Download fastq files and human reference genome 
### Download the raw fastq files from the SRA PRJNA1023207
The raw fastq files may be obtained from SRA PRJNA1023207. There are 619 individuals. Samples were sequenced to ~50 million (M) 2 × 100 bp paired-end reads across two lanes. Total storage requirements for the raw gzip sequences exceed 3TB. Information on how to download from SRA may be found [here](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/). 

Download or move the raw sequences to the raw_fastq folder. 

### Download the human GRCh38 reference genome and create sex chromosome complement reference genomes. 
Reads were aligned to the default gencode human GRCh38 reference genome. Samples were sex checked, and then samples were re-sequenced to a reference genome informed on the sex chromosome complement (SCC) of the sample. See [Olney et al. 2020 Biol Sex Differ](https://bsd.biomedcentral.com/articles/10.1186/s13293-020-00312-9) for more details about the SCC approach. 

Download the gencode GRCh38 reference genome and gene annotation file. First, change the working directory to the references folder OR once the download is complete, move the files to the references folder. 
```
cd references
```
```
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
```
Change the working directory to the scripts/snakemake folder. 
```
cd scripts/snakemake/
```
The build_gencode.v38.Snakefile script located in `scripts/snakemake/` will create three STAR reference indexes: default, YPARSmask, and Ymasked. Be sure to have the LBD conda environment activated before running snakemake. The human_ref_config.json has already been created and is located in the scripts/snakemake/ folder. 
```
snakemake -s build_gencode.v38.Snakefile
```


## Run Snakemake for bulk RNAseq data processing 
The following scripts are located in `scripts/snakemake/`
1. The `get_read_group_info.sh` script will list the fastq files that are located in the raw_fastq folder, and then collect header information that will be used to create the config file in the next step. 
```
sh get_read_group_info.sh
```

2. The `create_RNA_config.py` python script will create a config file of the sample IDs that are outputted from running the step above. 
```
python create_RNA_config.py
```

3. Finally, the `RNA.alignment.Snakefile` will merge samples that were sequenced on multiple lanes, generate fastqc reports, align to the reference genome, and collect RNA metrics. It is highly recommended to run snakemake on an HPC cluster. 
```
snakemake -s RNA.alignment.Snakefile -j 100 --rerun-incomplete 
```
The above command will submit 100 jobs in parallel and re-run any incomplete jobs. The outputs will include fastqc reports, BAM files, gene-level counts data, and RNA metrics. Details regarding each processing step are outlined within the `RNA.alignment.Snakefile`. 

To have the sequences align to a reference genome informed on the sex chromosome complement, create a male and female sample ID list within the `RNA.config.json` file. See the premade `RNA.config.json` file located in `scripts/snakemake/` for the list of male and female sample IDs. The individuals have already been sex checked, see Olney et al. 202x for more details. 

## Bulk RNAseq gene differential expression analysis 
Counts data obtained from `RNA.alignment.Snakefile` along with sample information is now ready to be read into R for further examination. The R scripts will source two  R files that contain libraries and visualization variables. 

Source files:
a. `file_paths_and_colors.R` This script will be sourced in the  R scripts listed below to load the necessary libraries, read in the metadata table, and define colors and shapes for each disease type and genetic sex.
b. `gtf_path.R` This script will be sourced in only some of the R scripts listed above. The script will load in the human gene annotation file that is needed for some of the downstream analyses. 

Overview of R scripts:
1. `01_sex_check.Rmd` evaluates the counts data obtained from aligning the reads to the default reference genome and determines if the reported sex matches the inferred sex for each individual. This script is optional and may be skipped as all samples have already been sex checked and reported sex matched inferred sex. 
2. `02_metadata_and_alignment_mertic_assessment.Rmd` evaluates and plots metadata and alignment metric information among the disease types and between the sexes. The alignment metrics have already been added to the `metadata.tsv` file and thus this script is optional and may be skipped. 
3. `03_create_dge_object.Rmd` will collect counts data for each sample and create an DGEList object which holds the dataset to be analyzed by edgeR and the subsequent calculations performed on the dataset. Specifically, the DGEList object contains counts, library size, norm.factors, group and metadata information, and gene annotation information. The 03a script will also remove MT genes, retain only protein coding genes, and filter to remove lowly expressed genes. 
4. `04_assess_variance.Rmd` will run variance partitioner to quantify sources of variation. Additionally, the script will perform a forward stepwise regression Bayesian information criterion (BIC) to determine the best model and then perform differential expression analysis among pairwise groups. The optional `04b_assess_variance_sex.Rmd` script assesses sources of variation in expression data and performs differential expression analysis within each genetic sex. 
5. `05_DEGs.Rmd` will generate volcano plots and differentially expressed gene (DEG) tables for all pairwise comparisons. The optional `05_DEGs_sex.Rmd` will generate volcano plots and DEGs tables for each sex separately and examine sex differential expression patterns. 
6. `06_WGCNA.Rmd`


## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | olney.kimberly@mayo.edu |
| John Fryer, PhD | fryer.john@mayo.edu |


