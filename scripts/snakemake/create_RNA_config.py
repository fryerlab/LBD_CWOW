#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

allSamples = list()
read = ["R1", "R2"]
numSamples = 0

with open('sample_read_group_info.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace(".", "_")
        split = line.split()
        sampleAttributes = split[0].split('_') # NA00-190.FCH5MYFDMXY_L1_R1_ITCTCTACT-CGCGGTTC.fastq.gz 
        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2]
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "rawReads" : "/LBD_CWOW/raw_fastq/",
    "rawQC" : "../../rawQC/",
    "trimmedReads" : "../../trimmedReads/",
    "trimmedQC" : "../../trimmedQC/",
    "starAligned" : "../../starAligned/",
    "starAligned_SCC" : "../../startAligned_SCC/",
    "bamstats" : "../../bamstats/",

    "Comment_Reference" : "This section specifies the location of the Sus scrofa, Ensembl reference genome",
    "GRCh38.fa" : "/LBD_CWOW/references/GRCh38.primary_assembly.genome.fa",
    "GRCh38.star" : "/LBD_CWOW/references/GRCh38_def",
    "GRCh38.gtf" : "/LBD_CWOW/references/gencode.v38.annotation",
    
    "GRCh38.Ymasked.fa" : "/LBD_CWOW/references/GRCh38_Ymasked_XX.fa",
    "GRCh38.YPARs_masked.fa" : "/LBD_CWOW/references/GRCh38_YPARsmasked_XY.fa",

    "GRCh38.Ymasked.star" : "/LBD_CWOW/references/GRCh38_Ymasked_XX",
    "GRCh38.YPARs_masked.star" : "/LBD_CWOW/references/GRCh38_YPARsmasked_XY",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
    "male_names": [],
    "female_names": [],
    "read": {1},
'''
outfile.write(header.format(allSamples, read))

# config formatting
counter = 0
with open('sample_read_group_info.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace(".fastq.gz", "")
        sampleName1 = base
        sampleName2 = sampleName1.replace("R1","R2")
        base = base.replace("_R1_", "")
        sampleInfo = split[1]

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace(".", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # NA00-190.FCH5MYFDMXY_L1_R1_ITCTCTACT-CGCGGTTC.fastq.gz
        # uniqueNum-number_sequencer_lane_read.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2]
        shortName1 = stemName + '_R1'
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        ID = stemName  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = stemName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "fq_path": "/LBD_CWOW/raw_fastq/",
        "fq1": "{1}",
        "fq2": "{2}",
        "shortName1": "{3}",
        "shortName2": "{4}",
        "ID": "{5}",
        "SM": "{5}",
        "PU": "{6}",
        "LB": "{7}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, shortName1, shortName2, stemName, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
