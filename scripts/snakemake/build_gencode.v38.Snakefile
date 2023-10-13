import os

configfile: "human_ref_config.json"

#Tools
star_path = "STAR"

rule all:
	input:
		expand("{gencode_ref_Ymasked}_XY.fa", gencode_ref_Ymasked = config["gencode_Ymasked"]),
		expand("{gencode_ref_YPARsmasked}_XX.fa", gencode_ref_YPARsmasked = config["gencode_YPARsmasked"]),
		expand("{gencode_ref_fa}_def_STAR", gencode_ref_fa = config["gencode_fa"]),
		expand("{gencode_ref_Ymasked}_XX_STAR", gencode_ref_Ymasked = config["gencode_Ymasked"]),
		expand("{gencode_ref_YPARsmasked}_XY_STAR", gencode_ref_YPARsmasked = config["gencode_YPARsmasked"])
      
#---------------------
# Reference genome and annotation were downloaded prior to running snakemake. 
# https://www.gencodegenes.org/human/
# Gencode GRCh38.p13 fasta
# 	wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
# Gencode GRCh38.p13 fasta gtf
# 	wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# Gencode GRCh38.p13 transcriptome cdna
#       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz 
# prep reference genome. 
# Hard mask Y chromosome for aligning XX samples. 
# Create index and dictionary. 

rule Ymask_ref:
    output:
        ref_out = "{gencode_ref_Ymasked}_XY.fa"
    shell:
    	"""
    	python maskChrY.py
    	"""

rule YPARsmask_ref:
    output:
        ref_out = "{gencode_ref_YPARsmasked}_XX.fa"
    shell:
    	"""
    	sh maskPARY.sh
    	"""
    	
rule star_def_genome_index:
	input:
		fa_def = "{gencode_ref_fa}.fa",
		gtf_star_def = (config["gencode_gtf"]+".gtf")
	output:
		star_index_def = "{gencode_ref_fa}_def_STAR",
	params:
		star = star_path,
		fa_def_dir = "{gencode_ref_fa}",
	shell:
		"""
		mkdir {params.fa_def_dir}_def_STAR;
		{params.star} --runThreadN 8 --runMode genomeGenerate --genomeDir {output.star_index_def} --genomeFastaFiles {input.fa_def} --sjdbGTFfile {input.gtf_star_def}
		"""
rule star_Ymask_genome_index:
	input:
		fa_Ymask = "{gencode_ref_Ymasked}_XX.fa",
		gtf_star = (config["gencode_gtf"]+".gtf")
	output:
		star_index_Ymask = "{gencode_ref_Ymasked}_XX_STAR",
	params:
		star = star_path,
		fa_Ymask_dir = "{gencode_ref_Ymasked}",
	shell:
		"""
		mkdir {params.fa_Ymask_dir}_XX_STAR;
		{params.star} --runThreadN 8 --runMode genomeGenerate --genomeDir {output.star_index_Ymask} --genomeFastaFiles {input.fa_Ymask} --sjdbGTFfile {input.gtf_star}
		"""

rule star_YPARsmask_genome_index:
	input:
		fa_YPARsmask = "{gencode_ref_YPARsmasked}_XY.fa",
		gtf_star = (config["gencode_gtf"]+".gtf")
	output:
		star_index_YPARsmask = "{gencode_ref_YPARsmasked}_XY_STAR",
	params:
		star = star_path,
		fa_YPARsmask_dir = "{gencode_ref_YPARsmasked}",
	shell:
		"""
		mkdir {params.fa_YPARsmask_dir}_XY_STAR;
		{params.star} --runThreadN 8 --runMode genomeGenerate --genomeDir {output.star_index_YPARsmask} --genomeFastaFiles {input.fa_YPARsmask} --sjdbGTFfile {input.gtf_star}
		"""
		
# KEY
#--sjdbOverhang ReadLength-1, default value of 100
#--runThreadN # of threads
#--runMode build genome indices (default is aligns reads)
#--genomeFastaFiles /path/to/genome/fasta
#--sjdbGTFfile /path/to/annotations.gtf
#--genomeDir /path/to/genomeDir specifies where the genome indices are stored. Note needs to be made prior to this step 
#---------------------
