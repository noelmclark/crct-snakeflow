### start with some stuff to get the Snakemake version
from snakemake import __version__ as snakemake_version
smk_major_version = int(snakemake_version[0])


### import modules as needed. Also import the snakemake
## internal command to load a config file into a dict
import pandas as pd

if smk_major_version >= 8:
  from snakemake.common.configfile import _load_configfile
else:
  from snakemake.io import _load_configfile





### Get a dict named config from config/config.yaml
configfile: "config/config.yaml"





### Get the sample info table read into a pandas data frame
sample_table=pd.read_table(config["sample_info"], dtype="str").set_index(
    ["sample", "unit"], drop=False
)





### Transfer values from the yaml and tabular config to
### our familiar lists, SAMPLES and CHROMOS
# Populate our SAMPLES list from the sample_table using a little
# pandas syntax
SAMPLES=sample_table["sample"].unique().tolist()

# Define chromosomes and scaffold groups from the values in the config file
chromosomes = pd.read_table(config["chromosomes"]).set_index("chrom", drop=False)
scaffold_groups = pd.read_table(config["scaffold_groups"]).set_index("id", drop=False)

# ensure that column 1 of the scaffold group file is "id" and
# column 2 is "chrom".  This is essential because Eric used those
# column positions in some awk to pull things out.
scaff_cols = list(scaffold_groups.columns)
if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom': 
    raise Exception("Column order is important in the scaffold_groups file.  The first column must be 'id' and the second column must be 'chrom'.")





### Input Functins that use the tabular sample_info

# define a function to get the fastq path from the sample_table. This
# returns it as a dict, so we need to unpack it in the rule
# This function groups fastqs by sample and unit so we can account for cases
# where there are more than 2 fastq files per sample. 
# Also the if function deals with if we only have one fastq file for a sample. 
def get_fastqs(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = sample_table.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

# define a function for getting the read group information
# from the sample table for each particular sample (according
# to the wildcard value)
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{sample_id}_{library}_{flowcell}_{lane}\tSM:{sample_id}\tPL:{platform}\tLB:{library}\tPU:{flowcell}.{lane}.{library}'".format(
        sample=sample_table.loc[(wildcards.sample, wildcards.unit), "sample"],
        sample_id=sample_table.loc[(wildcards.sample, wildcards.unit), "sample_id"],
        platform=sample_table.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=sample_table.loc[(wildcards.sample, wildcards.unit), "library"],
        flowcell=sample_table.loc[(wildcards.sample, wildcards.unit), "flowcell"],
        lane=sample_table.loc[(wildcards.sample, wildcards.unit), "lane"],
    )

#define function to get all the bam files for different units of same sample
def get_all_bams_of_common_sample(wildcards):
    s=sample_table.loc[(sample_table["sample"] == wildcards.sample)]
    return(expand("results/mapped/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist()
    ))





### Deal with the indel_grps if present (i.e. groupings of the samples
### into different species so that indel realignment is done species by species).
### At the same time we do this, we are also going to define the lists of output bams

# this is the default, but we update it if the indel_grps is defined in the config
realigned_bams_output_list=expand("results/indel_realigned/__ALL/{samp}.bam", samp=sample_list)
indel_grps_list=["weirdwhacko", "__ALL"]  # this is for the wildcard constraints. For some reason it fails if it is just "__ALL"
if "indel_grps" in config and config["indel_grps"] != "":
    indel_grps = pd.read_table(config["indel_grps"], dtype=str).set_index("sample", drop=False)
    validate(indel_grps, schema="../schemas/indel_grps.schema.yaml")

    # check to make sure that every sample is accounted for
    if(set(sample_list) != set(indel_grps["sample"].unique())):
      raise Exception("Every sample in units must be represented in the indel_grps file")

    # update the realigned_bams_output_list
    sm=indel_grps['sample'].tolist()
    ig=indel_grps['indel_grp'].tolist()
    realigned_bams_output_list=["results/indel_realigned/{ig}/{samp}.bam".format(ig=ig[i], samp=sm[i]) for i in range(len(ig))]
    indel_grps_list=ig


# down here, if we choose not to do indel realignment in the config
# then we just set realigned_bams_output_list to an empty list
if "do_indel_realignment" in config and config["do_indel_realignment"] == False:
    realigned_bams_output_list = []





### loose ends from eric's common.smk
def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL"
    )


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}