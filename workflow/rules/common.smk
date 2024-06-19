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



### Load MAF cutoffs from the config file
# these are our MAF cutoffs to prepare.  We take the unique values
# of the the maf_cutoffs from the config.
mafs = list(dict.fromkeys([str(x) for x in config["maf_cutoffs"]]))



### Get the sample info table read into a pandas data frame
sample_table=pd.read_table(config["sample_info"], dtype="str").set_index(
    ["sample", "unit"], drop=False
)

# rather than have a separate samples.tsv, we can just get a list of
# the samples from sample_table
sample_list = sample_table["sample"].unique().tolist()

# this is handy for getting sample info from the bcftools summaries
unique_sample_ids = list(sample_table["sample_id"].unique())

# Define chromosomes and scaffold groups from the values in the config file
chromosomes = pd.read_table(config["chromosomes"]).set_index("chrom", drop=False)
scaffold_groups = pd.read_table(config["scaffold_groups"]).set_index("id", drop=False)

# ensure that column 1 of the scaffold group file is "id" and
# column 2 is "chrom".  This is essential because Eric used those
# column positions in some awk to pull things out.
scaff_cols = list(scaffold_groups.columns)
if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom': 
    raise Exception("Column order is important in the scaffold_groups file.  The first column must be 'id' and the second column must be 'chrom'.")

# get a list of just the unique values of the scaffold_group and of the chromosomes
unique_scaff_groups = list(scaffold_groups.id.unique())
unique_chromosomes = list(chromosomes.chrom.unique())  # don't need to unique it, but I do anyway
sg_or_chrom = list(unique_scaff_groups + unique_chromosomes)

#get list of desired max lengths of the binsizes for the scatter_intervals_file R script from the config
binsize = list(dict.fromkeys([str(x) for x in config["binsize"]]))

# finally, get all the scatter groups, indexed two different ways.
scatter_wc_constraint="scat_0[0-9]*"  # this is just here for the case where `scatter_intervals_file: ""`
if config["scatter_intervals_file"] != "":
    scatter_groups = pd.read_table(config["scatter_intervals_file"]).set_index("id", drop=False)
    scatter_cols = list(scatter_groups.columns)
    if scatter_cols[0] != 'id' or scatter_cols[1] != 'scatter_idx' or scatter_cols[2] != 'chrom' or scatter_cols[3] != 'start' or scatter_cols[4] != 'end' or scatter_cols[5] != 'scatter_length':
        raise Exception("Column order is important in the scaffold_groups file.  The columns must be in order: id, scatter_idx, chrom, start, end, scatter_length.")
    unique_scats = list(scatter_groups.scatter_idx.unique())
    scatter_wc_constraint="|".join(unique_scats)
    # get a pandas frame of unique values of scaff_group and scat. This is for force-calling VCFs
    unique_scatters_table=scatter_groups[['id', 'scatter_idx']].drop_duplicates()


#scat_ids=scatter_groups.loc[(scatter_groups["id"] == sg_or_chrom), (scatter_groups["scatter_idx"] == unique_scats)].unique()
#scatter_groups.loc["id", "scatter_idx"].unique()


##### Wildcard constraints #####
# copied from Eric
# we have to deal with the cases where scaff groups or chroms might be empty
if not unique_scaff_groups:
    sg_constraint = "CrazyScaffGroupName"  # this is sort of a hack.  Just put something that won't be a chromosome
else:
    sg_constraint="|".join(unique_scaff_groups)

if not unique_chromosomes:
    chrom_constraint="CrazyChromName"
else:
    chrom_constraint="|".join(unique_chromosomes)

wildcard_constraints:
    sample="|".join(sample_list),
    unit="|".join(sample_table["unit"]),
    chromo=chrom_constraint,
    scaff_group=sg_constraint,
    sg_or_chrom="|".join(unique_scaff_groups + unique_chromosomes),
    filter_condition="ALL|PASS|FAIL",
    maf="|".join(mafs),
    binsize="|".join(binsize),
    scatter=scatter_wc_constraint,
    #igrp="|".join(indel_grps_list)



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
    return(expand("results/mapping/mapped/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist()
    ))

# given a chromosome or a scaff_group, get all the scattered vcf.gz of vcf.gz.tbi files
#def get_scattered_vcfs(wildcards, ext):
#    scat_ids=scatter_groups.loc[(scatter_groups["id"] == wildcards.sg_or_chrom), (scatter_groups["scatter_idx"] == wildcards.unique_scats)].unique()
#    return(expand("results/calling/vcf_sections/{sgc}/{scat}.vcf.gz{e}", 
#        sgc=sg_or_chrom, 
#        scat=scat_ids, 
#        e=ext
#    ))

def get_scattered_vcfs(wildcards, ext):
    scat_ids=scatter_groups.loc[(scatter_groups["id"] == wildcards.sg_or_chrom), "scatter_idx"].unique()
    return expand("results/calling/vcf_sections/{{sg_or_chrom}}/{scat}.vcf.gz{e}", scat=scat_ids, e=ext)


# here, we deal with indel realignment by species (or "igrp")
def get_igrp_sample_names(wildcards):
    if "indel_grps" not in config or config["indel_grps"] == "":
        if wildcards.igrp != "__ALL":
            raise Exception("Requesting igrp that is not \"__ALL\", but indel_grps not defined in config.")
        else:
            ret="JustStuffNotNeeded: IDs gotten through Unix in the rule if __ALL"
    else:
        if wildcards.igrp == "__ALL":
          raise Exception("You can't request igrp \"__ALL\" when indel_grps is defined in the config. ")
        else:
          ret=indel_grps.loc[indel_grps['indel_grp'] == wildcards.igrp]["sample_id"].unique().tolist()
    return ret

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