# LBD_CWOW
Lewy Body Dementia Center With Out Walls (LBD CWOW) RNAseq processing and analysis.
Anterior cingulate cortex tissue samples from the Mayo Clinic brain bank were collected for 619 individuals. 
Raw fastq files are available on SRA, PRJNA1023207

| Disease                   | Count   |
| ------------------------- |:-------:|
| Control                   | 86      |
| Pathological  aging (PA)  | 39      |
| Alzheimer’s disease (AD)  | 54      |
| Lewy body dementia (LBD)  | 440     |

This git repo contains scripts for the following:
-   Metadata analysis
-   Processing of bulk-tissue RNA-sequencing data
-   Analysis of bulk-tissue RNA-sequencing data 
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
After the conda LBD environment has been created, you will need to additionally install GATK
This step must be done manually and not through conda; see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
The above link will explain how to download GATK4, then you will need to add an alias to your bash profile:
```
alias gatk='/path/to/gatk-package/gatk'
```
## Download fastq files and human reference genome 
### Download the raw fastq files from the SRA PRJNA1023207

### Download the human GRCh38 reference genome and create sex chromosome complement reference genomes. 

## Run Snakemake for bulk RNAseq  data processing 
1. get_read_group_info.sh script will list the fastq files that are located in the raw_fastq folder and collect header information that will be used to create the config file in the next step. 
```
sh scripts/snakemake/get_read_group_info.sh
```

2. create_RNA_config.py will create a config file of the sample IDs that are outputted from running the step above. 
```
python create_RNA_config.py
```

3. Run Snakemake which will merge samples that were sequenced on multiple lanes, generate fastqc reports, align to the reference genome, and collect RNA metrics. It is highly recommended to run snakemake in script on an HPC cluster.
```
snakemake -s RNA.alignment.Snakefile -j 100 --rerun-incomplete 
```
The above command will submit 100 jobs in parallel and re-run any incomplete jobs. 



## Bulk RNAseq gene differential expression analysis 

## Manuscript figures 
Scripts to make the following figures can be found under `scripts/R/manuscript_figures/`

Main figures:
1. **Figure 1:**
2. **Figure 2:** 
3. **Figure 3:** 
4. **Figure 4:** 
5. **Figure 5:** 

Supplemental figures:
1. **Fig S1:** 

## Contact information 
## Contact information

| Contact | Email |
| --- | --- | --- |
| Kimberly Olney, PhD | olney.kimberly@mayo.edu |
| John Fryer, PhD | fryer.john@mayo.edu |


