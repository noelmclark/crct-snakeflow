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

# sample table for replacing read group info indexed by sample and unit 
# (this is basically the same as the table above but I wanted two with different names to keep track)
replace_RG_table=pd.read_table(config["sample_info"], dtype="str").set_index(
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
    post_id="ALL",
    psmc_id="lineage",
    subsamp="all|crct-blue|crct-green|crct-both|srm|outgroups",



##### Input Functins #####

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
    return r"-R '@RG\tID:{sample}_{sample_id}_{library}_{flowcell}_{lane}\tSM:{sample}\tPL:{platform}\tLB:{library}\tPU:{flowcell}.{lane}.{library}'".format(
        sample=sample_table.loc[(wildcards.sample, wildcards.unit), "sample"],
        sample_id=sample_table.loc[(wildcards.sample, wildcards.unit), "sample_id"],
        platform=sample_table.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=sample_table.loc[(wildcards.sample, wildcards.unit), "library"],
        flowcell=sample_table.loc[(wildcards.sample, wildcards.unit), "flowcell"],
        lane=sample_table.loc[(wildcards.sample, wildcards.unit), "lane"],
    )

# define a function for getting the read group information 
# from the sample table for each particular sample (according
# to the wildcard value) for use in the AddOrReplaceReadGroups rule 
# except this is not working properly so I am going to put it in the shell of the rule
def replace_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-ID {sample}_{sample_id}_{library}_{flowcell}_{lane} -SM {sample} -PL {platform} -LB {library} -PU {flowcell}.{lane}.{library}".format(
        sample=replace_RG_table.loc[(wildcards.sample, wildcards.unit), "sample"],
        sample_id=replace_RG_table.loc[(wildcards.sample, wildcards.unit), "sample_id"],
        platform=replace_RG_table.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=replace_RG_table.loc[(wildcards.sample, wildcards.unit), "library"],
        flowcell=replace_RG_table.loc[(wildcards.sample, wildcards.unit), "flowcell"],
        lane=replace_RG_table.loc[(wildcards.sample, wildcards.unit), "lane"],
    )

#define function to get all the bam files for different units of same sample
def get_all_bams_of_common_sample(wildcards):
    s=sample_table.loc[(sample_table["sample"] == wildcards.sample)]
    return(expand("results/mapping/mapped/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist()
    ))

#define function to ask for all the RG fixed bams in rule all for an intermediate step
#if the RGs are set properly from the beginning, you don't need this
def get_RG_fixed_bams_of_common_sample(wildcards):
    s=sample_table.loc[(sample_table["sample"] == wildcards.sample)]
    return(expand("results/angsd_bams/RG-fixed/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist()
    ))

# given a chromosome or a scaff_group, get all the scattered vcf.gz of vcf.gz.tbi files
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


# for setting options for genomics DB import
def chromo_import_gdb_opts(wildcards):
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --intervals {chr} ".format(chr = wildcards.chromo))

def scaff_group_import_gdb_opts(wildcards):
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --intervals results/calling/interval_lists/{sg}.list --merge-contigs-into-num-partitions 1 ".format(sg = wildcards.scaff_group))

def get_bcftools_stats_filter_option(wildcards):
    if(wildcards.filter_condition == "ALL"):
        return(" ")
    elif(wildcards.filter_condition == "PASS"):
        return( " -i 'FILTER=\"PASS\"' ")
    elif(wildcards.filter_condition == "FAIL"):
        return( " -i 'FILTER!=\"PASS\"' ")
    else:
        raise Exception("Wildcard filter_condition must be ALL, PASS, or FAIL.")

def get_sgc_list(wildcards):
    if wildcards.sample == wildcards.sample:
        ret=sg_or_chrom
    else: 
        raise Exception("That didn't work.")
    return ret

def get_comma_sep_sample_names(sample_list):
    sample_list=sample_table["sample"].unique().tolist()
    return ','.join(sample_list)


##### PSMC grouping things #####
### Get a PSMC sample info table read into a pandas data frame
psmc_table=pd.read_table(config["psmc_info"], dtype="str").set_index(
    ["sample"], drop=False
)

psmc_id="lineage"
subsamp=["all", "crct-blue", "crct-green", "crct-both", "srm", "outgroups"]
psmcbootround = list(range(1, 101)) #specify wildcard for psmc bootstrap rounds, that is just a list of 1-100

lineage_list = psmc_table["lineage"].unique().tolist()
population = psmc_table["population"].unique().tolist()

def get_psmc_subsamps(wildcards):
    if(wildcards.subsamp == "all"):
        return(expand("results/psmc/run-psmc/{s}.psmc", s=sample_list))
    elif(wildcards.subsamp == "crct-blue"):
        b=psmc_table.loc[(psmc_table["lineage"] == "BL")]
        return(expand("results/psmc/run-psmc/{s}.psmc", s=b["sample"].tolist()))
    elif(wildcards.subsamp == "crct-green"):
        g=psmc_table.loc[(psmc_table["lineage"] == "GR")]
        return(expand("results/psmc/run-psmc/{s}.psmc", s=g["sample"].tolist()))
    elif(wildcards.subsamp == "crct-both"):
        bgvalues=["BL", "GR"]
        bg=psmc_table.loc[(psmc_table["lineage"].isin(bgvalues))]
        return(expand("results/psmc/run-psmc/{s}.psmc", s=bg["sample"].tolist()))
    elif(wildcards.subsamp == "srm"):
        srmvalues=["BL", "GR", "GB", "RG", "SJ"]
        srm=psmc_table.loc[(psmc_table["lineage"].isin(srmvalues))]
        return(expand("results/psmc/run-psmc/{s}.psmc", s=srm["sample"].tolist()))
    elif(wildcards.subsamp == "outgroups"):
        bgvalues=["BL", "GR"]
        out=psmc_table.loc[(~psmc_table["lineage"].isin(bgvalues))]
        return(expand("results/psmc/run-psmc/{s}.psmc", s=out["sample"].tolist()))
    else:
        raise Exception("Wildcard subsamp must be all, crct-blue, crct-green, crct-both, srm, or outgroups.")

def get_comma_sep_subsamp_names(wildcards):
    if(wildcards.subsamp == "all"):
        return','.join(sample_list)
    elif(wildcards.subsamp == "crct-blue"):
        b=psmc_table.loc[(psmc_table["lineage"] == "BL")]
        return','.join(b["sample"].tolist())
    elif(wildcards.subsamp == "crct-green"):
        g=psmc_table.loc[(psmc_table["lineage"] == "GR")]
        return','.join(g["sample"].tolist())
    elif(wildcards.subsamp == "crct-both"):
        bgvalues=["BL", "GR"]
        bg=psmc_table.loc[(psmc_table["lineage"].isin(bgvalues))]
        return','.join(bg["sample"].tolist())
    elif(wildcards.subsamp == "srm"):
        srmvalues=["BL", "GR", "GB", "RG", "SJ"]
        srm=psmc_table.loc[(psmc_table["lineage"].isin(srmvalues))]
        return','.join(srm["sample"].tolist())
    elif(wildcards.subsamp == "outgroups"):
        bgvalues=["BL", "GR"]
        out=psmc_table.loc[(~psmc_table["lineage"].isin(bgvalues))]
        return','.join(out["sample"].tolist())
    else:
        raise Exception("Wildcard subsamp must be all, crct-blue, crct-green, crct-both, srm, or outgroups.")

def get_comma_sep_pop_names(wildcards):
    if(wildcards.population == "navajo"):
        nav=psmc_table.loc[(psmc_table["population"] == "navajo")]
        return','.join(nav["sample"].tolist())
    elif(wildcards.population == "williamson"):
        wil=psmc_table.loc[(psmc_table["population"] == "williamson")]
        return','.join(wil["sample"].tolist())
    elif(wildcards.population == "steelman"):
        ste=psmc_table.loc[(psmc_table["population"] == "steelman")]
        return','.join(ste["sample"].tolist())
    elif(wildcards.population == "nanita"):
        nan=psmc_table.loc[(psmc_table["population"] == "nanita")]
        return','.join(nan["sample"].tolist())
    elif(wildcards.population == "como"):
        com=psmc_table.loc[(psmc_table["population"] == "como")]
        return','.join(com["sample"].tolist())
    elif(wildcards.population == "s_hayden"):
        sha=psmc_table.loc[(psmc_table["population"] == "s_hayden")]
        return','.join(sha["sample"].tolist())
    elif(wildcards.population == "severy"):
        sev=psmc_table.loc[(psmc_table["population"] == "severy")]
        return','.join(sev["sample"].tolist())
    elif(wildcards.population == "abrams"):
        abr=psmc_table.loc[(psmc_table["population"] == "abrams")]
        return','.join(abr["sample"].tolist())
    elif(wildcards.population == "e_fk_piedra"):
        efk=psmc_table.loc[(psmc_table["population"] == "e_fk_piedra")]
        return','.join(efk["sample"].tolist())
    elif(wildcards.population == "w_fk_boulder"):
        wfk=psmc_table.loc[(psmc_table["population"] == "w_fk_boulder")]
        return','.join(wfk["sample"].tolist())
    elif(wildcards.population == "kelso"):
        kel=psmc_table.loc[(psmc_table["population"] == "kelso")]
        return','.join(kel["sample"].tolist())
    elif(wildcards.population == "roan"):
        roa=psmc_table.loc[(psmc_table["population"] == "roan")]
        return','.join(roa["sample"].tolist())
    elif(wildcards.population == "dry_gulch"):
        dry=psmc_table.loc[(psmc_table["population"] == "dry_gulch")]
        return','.join(dry["sample"].tolist())
    elif(wildcards.population == "s_twin"):
        stw=psmc_table.loc[(psmc_table["population"] == "s_twin")]
        return','.join(stw["sample"].tolist())
    elif(wildcards.population == "hunter"):
        hun=psmc_table.loc[(psmc_table["population"] == "hunter")]
        return','.join(hun["sample"].tolist())
    elif(wildcards.population == "w_antelope"):
        wan=psmc_table.loc[(psmc_table["population"] == "w_antelope")]
        return','.join(wan["sample"].tolist())
    elif(wildcards.population == "yellowstone"):
        yel=psmc_table.loc[(psmc_table["population"] == "yellowstone")]
        return','.join(yel["sample"].tolist())
    elif(wildcards.population == "san_juan"):
        san=psmc_table.loc[(psmc_table["population"] == "san_juan")]
        return','.join(san["sample"].tolist())
    elif(wildcards.population == "bonneville"):
        bon=psmc_table.loc[(psmc_table["population"] == "bonneville")]
        return','.join(bon["sample"].tolist())
    elif(wildcards.population == "greenback"):
        gre=psmc_table.loc[(psmc_table["population"] == "greenback")]
        return','.join(gre["sample"].tolist())
    elif(wildcards.population == "lahontan"):
        lah=psmc_table.loc[(psmc_table["population"] == "lahontan")]
        return','.join(lah["sample"].tolist())
    elif(wildcards.population == "coastal"):
        coa=psmc_table.loc[(psmc_table["population"] == "coastal")]
        return','.join(coa["sample"].tolist())
    elif(wildcards.population == "westslope"):
        wes=psmc_table.loc[(psmc_table["population"] == "westslope")]
        return','.join(wes["sample"].tolist())
    elif(wildcards.population == "rio_grande"):
        rio=psmc_table.loc[(psmc_table["population"] == "rio_grande")]
        return','.join(rio["sample"].tolist())
    else:
        raise Exception("Wildcard population must be one of those listed in psmc-info.tsv.")

def get_psmc_boot_subsamps(wildcards):
    if(wildcards.subsamp == "all"):
        return(expand("results/psmc/bootstrap/run-psmc/{s}/{s}-100bootstrap.psmc", s=sample_list))
    elif(wildcards.subsamp == "crct-blue"):
        b=psmc_table.loc[(psmc_table["lineage"] == "BL")]
        return(expand("results/psmc/bootstrap/run-psmc/{s}/{s}-100bootstrap.psmc", s=b["sample"].tolist()))
    elif(wildcards.subsamp == "crct-green"):
        g=psmc_table.loc[(psmc_table["lineage"] == "GR")]
        return(expand("results/psmc/bootstrap/run-psmc/{s}/{s}-100bootstrap.psmc", s=g["sample"].tolist()))
    elif(wildcards.subsamp == "crct-both"):
        bgvalues=["BL", "GR"]
        bg=psmc_table.loc[(psmc_table["lineage"].isin(bgvalues))]
        return(expand("results/psmc/bootstrap/run-psmc/{s}/{s}-100bootstrap.psmc", s=bg["sample"].tolist()))
    elif(wildcards.subsamp == "srm"):
        srmvalues=["BL", "GR", "GB", "RG", "SJ"]
        srm=psmc_table.loc[(psmc_table["lineage"].isin(srmvalues))]
        return(expand("results/psmc/bootstrap/run-psmc/{s}/{s}-100bootstrap.psmc", s=srm["sample"].tolist()))
    elif(wildcards.subsamp == "outgroups"):
        bgvalues=["BL", "GR"]
        out=psmc_table.loc[(~psmc_table["lineage"].isin(bgvalues))]
        return(expand("results/psmc/bootstrap/run-psmc/{s}/{s}-100bootstrap.psmc", s=out["sample"].tolist()))
    else:
        raise Exception("Wildcard subsamp must be all, crct-blue, crct-green, crct-both, srm, or outgroups.")




##### hPSMC grouping things #####

bams=pd.read_table(config["hpsmc-test"]["bams"], dtype=str).set_index(
    ["sample"], drop=False
)

pwcomps=pd.read_table(config["hpsmc-test"]["pwcomps"], dtype=str).set_index(
    ["pop1", "pop2"], drop=False
)

hpsmcchroms=pd.read_table(config["hpsmc-test"]["chroms"], dtype=str).set_index(
    ["chrom"], drop=False
)

hpsmcpops=list(set(pwcomps.pop1.tolist() + pwcomps.pop2.tolist()))
pop1=pwcomps.pop1.tolist()
pop2=pwcomps.pop2.tolist()

hpsmcchroms=hpsmcchroms.chrom.tolist()

def get_hpsmc_bams_in_pop(wc):
  b=bams.loc[(bams["group"] == wc.hpsmcpops)]
  return b.path.tolist()