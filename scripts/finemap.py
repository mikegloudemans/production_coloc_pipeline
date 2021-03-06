# Author: Mike Gloudemans
#
# Run FINEMAP algorithm for colocalization.
# Store full results in the output directory for this run.
# Return the CLPP score.
#

import sys
import subprocess
from scipy import stats
from shutil import copyfile
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
import pandas as pd
import numpy as np
import math

def run_finemap(locus, window=500000):

    pf = prep_finemap(locus, window)
    if pf == "Fail":
        return "Fail"
    return launch_finemap(locus, window, pf)


# Separates the preparation of data for finemap from the actual 
# launching of finemap. This is done to avoid duplicating code,
# since eCAVIAR and FINEMAP use very similar setups.
def prep_finemap(locus, window):

    locus.conditional_level = 0   # Currently not used at all, but we may re-add this functionality later.
    combined = locus.data.copy()

    # Make required directories for FINEMAP eCAVIAR analysis
    subprocess.call("mkdir -p {0}/plink/{1}/{2}_{3}/{4}".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p {0}/ecaviar/{1}/{2}_{3}/{4}".format(locus.tmpdir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix), shell=True)
    subprocess.call("mkdir -p {0}/finemap".format(locus.basedir), shell=True)
    
    # Get VCF file paths
    gwas_vcf = locus.vcf_file.format(locus.chrom)
    eqtl_vcf = locus.vcf_file.format(locus.chrom)
    
    # Two different cases depending on whether GWAS and eQTL
    # are using same reference genome.
    if eqtl_vcf == gwas_vcf:
        # Get and filter the single VCF.
        vcf, combined = load_and_filter_variants(eqtl_vcf, locus, combined, window, ["eqtl", "gwas"])
        assert vcf.shape[0] == combined.shape[0]

        # Run PLINK on just one VCF.
        removal_list = compute_ld(vcf, locus, "eqtl")
        if removal_list == "Fail":
            return "Fail"
        subprocess.check_call("cp {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True)

        # Remove indices that produced NaNs in the LD computations
        removal_list = list(set(removal_list))
        combined = combined.drop(combined.index[removal_list])

    else:
        # Get and filter both VCFs.
        evcf, combined = load_and_filter_variants(eqtl_vcf, locus, combined, window, ["eqtl"])
        gvcf, combined = load_and_filter_variants(gwas_vcf, locus, combined, window, ["gwas"])

        # Subset to overlapping SNPs
        evcf, gvcf, combined = intersect_reference_vcfs(evcf, gvcf, combined)
        assert evcf.shape[0] == combined.shape[0]
        assert gvcf.shape[0] == combined.shape[0]

        # Run PLINK on both VCFs.
        while True:
            removal_list = compute_ld(evcf, locus, "eqtl")
            if removal_list == "Fail":
                return "Fail"
            
            extension_list = compute_ld(gvcf, locus, "gwas")
            if extension_list is "Fail":
                return "Fail"
 
            removal_list.extend(extension_list)

            # Continue until no more NaNs
            if len(removal_list) == 0:
                break

            removal_list = list(set(removal_list))

            combined = combined.drop(combined.index[removal_list])
            evcf = evcf.drop(evcf.index[removal_list])
            gvcf = gvcf.drop(gvcf.index[removal_list])

    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "w") as w:
        snps = combined[['snp_pos', 'ZSCORE_eqtl']]
        snps.to_csv(w, index=False, header=False, sep=" ")

    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "w") as w:
        snps = combined[['snp_pos', 'ZSCORE_gwas']]	
        snps.to_csv(w, index=False, header=False, sep=" ")
    
    return (min(combined["pvalue_gwas"]), min(combined["pvalue_eqtl"]))

# This function contains the code that's specific to FINEMAP,
# not shared with eCAVIAR.
def launch_finemap(locus, window, top_hits):

    # Load sample sizes
    eqtl_n = locus.N
    gwas_n = locus.N

    # Write config file for finemap
    subprocess.check_call('echo "z;ld;snp;config;n-ind" > {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True)
    subprocess.check_call('echo "{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.z;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_gwas.fixed.ld;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.config;{6}" >> {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, gwas_n, locus.tmpdir), shell=True)
    subprocess.check_call('echo "{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.z;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_eqtl.fixed.ld;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp;{7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.config;{6}" >> {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, eqtl_n, locus.tmpdir), shell=True)
    
    # Run FINEMAP
    subprocess.check_call('../bin/finemap --sss --in-files {6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_finemap.in --n-causal-max 1 --n-iterations 1000000 --n-convergence 50000 > /dev/null'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), shell=True)

    # Parse FINEMAP results to compute CLPP score
    gwas_probs = []
    eqtl_probs = []
    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
            f.readline()
            for line in f:
                    data = line.strip().split()
                    gwas_probs.append((int(data[1]), float(data[2])))
    with open("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir)) as f:
            f.readline()
            for line in f:
                    data = line.strip().split()
                    eqtl_probs.append((int(data[1]), float(data[2])))
    gwas_probs = sorted(gwas_probs)
    eqtl_probs = sorted(eqtl_probs)

    assert len(gwas_probs) == len(eqtl_probs)
    for i in range(len(gwas_probs)):
            assert gwas_probs[i][0] == eqtl_probs[i][0]

    finemap_clpp = sum([gwas_probs[i][1] * eqtl_probs[i][1] for i in range(len(gwas_probs))])

    # Save full FINEMAP results for the locus
    copyfile("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.gwas.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "{0}/finemap/{1}_{2}_{3}_{4}_{5}_{6}_finemap_gwas.snp".format(locus.basedir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level))
    copyfile("{6}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.finemap.eqtl.snp".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.tmpdir), "{0}/finemap/{1}_{2}_{3}_{4}_{5}_{6}_finemap_eqtl.snp".format(locus.basedir, locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level))

    # Write FINEMAP results to the desired file
    with open("{0}/{1}_finemap_clpp_status.txt".format(locus.basedir, locus.gwas_suffix.replace(".", "_")), "a") as a:
            a.write("{0}_{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(locus.chrom, locus.pos, locus.eqtl_suffix, locus.gwas_suffix, locus.gene, len(gwas_probs), finemap_clpp, -1*math.log10(top_hits[0]), -1*math.log10(top_hits[1])))

    return finemap_clpp

# Function purge_tmp_files
# Remove temporary files created during this run of eCAVIAR,
# to free up space.
def purge_tmp_files(locus):
    subprocess.call("rm -rf {4}/ecaviar/{0}/{1}_{2}/{3}/{5}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir, locus.gene), shell=True)
    subprocess.call("rm -rf {4}/plink/{0}/{1}_{2}/{3}/{5}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.tmpdir, locus.gene), shell=True)

def load_and_filter_variants(filename, locus, combined, window, ref_types):
    
    # First, extract nearby variants using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 500 | grep \\#CHROM".format(filename), shell=True).strip().split()
    stream = StringIO(subprocess.check_output("tabix {8} {1}:{6}-{7}".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, locus.pos-window, locus.pos+window, filename), shell=True))

    # For readability, load the header too
    # Load with pandas
    vcf = pd.read_csv(stream, sep="\t", names=header)

    # Remove variants not in the GWAS table
    vcf["POS"] = (vcf["POS"]).astype(int)
    vcf = vcf[vcf["POS"].isin(list(combined["snp_pos"]))]

    # Remove variants with position appearing multiple times
    dup_counts = {}
    for pos in vcf["POS"]:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1
    vcf["dup_counts"] = [dup_counts[pos] for pos in vcf['POS']]
    vcf = vcf[vcf["dup_counts"] == 1]

    # Remove multiallelic variants with only one entry in VCF
    l = lambda x: "," not in x
    vcf = vcf[vcf["REF"].apply(l) & vcf["ALT"].apply(l)]

    # Remove monoallelic variants.
    # Allele frequency might be input as counts or as percentages,
    # so handle this.

    af_id = "AF"
    def fn(x):
        info = [s for s in x.split(";") if s.startswith(af_id + "=")][0]
        af = float(info.split("=")[1])
        return af > 0.01 and 1-af > 0.01 
    vcf = vcf[vcf["INFO"].apply(fn)]

    # Remove variants where alt/ref don't match between GWAS/eQTL and VCF
    # Flipped is okay. A/C and C/A are fine, A/C and A/G not fine.
    merged = pd.merge(combined, vcf, left_on="snp_pos", right_on="POS")
    keep_indices = \
            (((merged['ref_gwas'] == merged['REF']) & (merged['alt_gwas'] == merged['ALT'])) | \
            ((merged['alt_gwas'] == merged['REF']) & (merged['ref_gwas'] == merged['ALT']))) & \
            (((merged['ref_eqtl'] == merged['REF']) & (merged['alt_eqtl'] == merged['ALT'])) | \
            ((merged['alt_eqtl'] == merged['REF']) & (merged['ref_eqtl'] == merged['ALT'])))

    # Now, reverse the z-score for any SNPs whose ref/alt direction doesn't match the direction
    # in the reference genome.
    assert "gwas" in ref_types or "eqtl" in ref_types
    if "gwas" in ref_types:
        merged['reverse_gwas'] = (merged['alt_gwas'] == merged['REF']) & (merged['ref_gwas'] == merged['ALT'])
        merged['ZSCORE_gwas'] = merged['ZSCORE_gwas'] * (1 - (merged['reverse_gwas'] * 2))
    if "eqtl" in ref_types:
        merged['reverse_eqtl'] = (merged['alt_eqtl'] == merged['REF']) & (merged['ref_eqtl'] == merged['ALT'])
        merged['ZSCORE_eqtl'] = merged['ZSCORE_eqtl'] * (1 - (merged['reverse_eqtl'] * 2))

    keep = merged['POS'][keep_indices]
    vcf = vcf[vcf['POS'].isin(list(keep))]
    
    # Subset SNPs down to SNPs present in the reference VCF.
    combined = combined[combined['snp_pos'].isin(list(vcf["POS"]))] 

    # Return list as DataFrame.
    return (vcf, combined)

def intersect_reference_vcfs(ref1, ref2, combined):

    # Subset each reference VCF down to only variants present in the other VCF.
    ref1 = ref1[ref1["POS"].isin(list(ref2["POS"]))]
    ref2 = ref2[ref2["POS"].isin(list(ref1["POS"]))]

    # Subset SNP list to SNPs that still remain in both VCFs.
    combined = combined[combined['snp_pos'].isin(list(ref2["POS"]))]

    return (ref1, ref2, combined)

# Run PLINK on the locus of interest
def compute_ld(input_vcf, locus, data_type):

    # We don't want to modify the input VCF within this function
    vcf = input_vcf.copy()

    # Repeatedly compute LD until we've eliminated all NaNs.
    removal_list = []
    while True:
        if vcf.shape[0] == 0:
            return "Fail"

        # Write VCF to tmp file
        vcf.to_csv('{7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.{6}.vcf'.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir), sep="\t", index=False, header=True)

        # Use PLINK to generate bim bam fam files
        command = '''../bin/plink --vcf {7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}.{6}.vcf --keep-allele-order --make-bed --double-id --out {7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}_plinked > /dev/null'''.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir)
        subprocess.check_call(command, shell=True)

        # Use PLINK to generate LD score
        command = '''../bin/plink -bfile {7}/plink/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}_plinked --r square --out {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6} > /dev/null'''.format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir)
        subprocess.check_call(command, shell=True)

        # See if nans remain. If so, remove the offending lines.
        try:
            lines = [int(n.split(":")[0])-1 for n in subprocess.check_output("grep -n nan {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir), shell=True).strip().split("\n")]
        except:
            break

        # Save IDs of items being removed
        removal_list.extend(list(vcf.iloc[lines]['ID']))
        
        # Remove desired rows (variants)
        vcf = vcf.drop(vcf.index[lines])

    # Get indices of items to remove in original list
    removal_list = [list(input_vcf['ID']).index(id) for id in removal_list]

    # Replace tabs with spaces because FINEMAP requires this.
    subprocess.check_call("cat {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.ld | sed s/\\\\t/\\ /g > {7}/ecaviar/{0}/{1}_{2}/{3}/{4}_fastqtl_level{5}_{6}.fixed.ld".format(locus.gwas_suffix, locus.chrom, locus.pos, locus.eqtl_suffix, locus.gene, locus.conditional_level, data_type, locus.tmpdir), shell=True)

    return removal_list
