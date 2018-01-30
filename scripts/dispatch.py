#!/usr/bin/python
# Author: Mike Gloudemans
# 
# Dispatch colocalization analyses
# for various GWAS and eQTL experiments,
# parameter settings, and colocalization
# methods.
#

# Built-in libraries
import sys
from shutil import copyfile
import datetime
import subprocess
import math
import os

# Custom libraries
import config
import preprocess
import SNP
import argparse
from TestLocus import TestLocus

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gwas_file", help="GWAS summary statistics", action="store", required=True)
    parser.add_argument("-e", "--eqtl_file", help="eQTL summary statistics", action="store", required=True)
    parser.add_argument("-c", "--chrom", help="chromosome", action="store", required=True)
    parser.add_argument("-s", "--snp_pos", help="central SNP coordinate", action="store", required=True)
    parser.add_argument("-v", "--vcf_file", help="reference VCF file for computing LD", action="store", required=True)
    parser.add_argument("-N", "--N", help="number of individuals in reference VCF", action="store", required=True)
    parser.add_argument("-w", "--window", help="window size for tested region in bp", action="store", default = 500000)
    parser.add_argument("--gene", help="eQTL gene of interest", action="store", default = -1)
    parser.add_argument("--eqtl_format", help="format of eQTL summary statistics file", choices=['effect_size', 'chisq', 'tstat'], action="store", default="effect_size")
    parser.add_argument("--gwas_format", help="format of GWAS summary statistics file", choices=['effect_size', 'log_odds_ratio', 'pval_only'], action="store", default="effect_size")
    args = parser.parse_args()

    # Change directory to the script's directory
    os.chdir((os.path.dirname(os.path.realpath(parser.prog))))
    

    settings = {}
    settings["gwas_file"] = args.gwas_file
    settings["eqtl_file"] = args.eqtl_file
    settings["chrom"] = int(args.chrom)
    settings["snp_pos"] = int(args.snp_pos)
    settings["vcf_file"] = args.vcf_file
    settings["N"] = int(args.N)
    settings["gene"] = args.gene
    settings["window"] = int(args.window)
    settings['eqtl_format'] = args.eqtl_format
    settings['gwas_format'] = args.gwas_format

    # Make timestamped results directory, under which all output for this run will be stored.
    now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    base_output_dir = "../output/{0}".format(now)
    base_tmp_dir = "../tmp/{0}".format(now)

    subprocess.call(["mkdir", "-p", base_output_dir])

    # Run FINEMAP pipeline
    analyze_snp(settings, base_output_dir, base_tmp_dir, settings["gene"])

    # Clean up after ourselves
    subprocess.call("rm -r {0}".format(base_tmp_dir), shell=True)


def analyze_snp(settings, base_output_dir, base_tmp_dir, restrict_gene=-1):

    snp = SNP.SNP((settings["chrom"], settings["snp_pos"], 1))

    # Load relevant GWAS and eQTL data.
    gwas_data = preprocess.get_gwas_data(settings, snp, settings["window"]) # Get GWAS data
    eqtl_data = preprocess.get_eqtl_data(settings, snp, settings["window"]) # Get eQTL data

    # Skip it if this entire locus has no genes
    if isinstance(eqtl_data, basestring):
        print "Error: No genes near this locus."
        return

    # Get all genes whose eQTLs we're testing at this locus
    if restrict_gene == -1:
        genes = set(eqtl_data['gene'])
    else:
        genes = [restrict_gene]

    # Write header
    gwas_suffix = settings["gwas_file"].split("/")[-1].replace(".", "_")
    with open("{0}/{1}_finemap_clpp_status.txt".format(base_output_dir, gwas_suffix), "w") as w:
        w.write("snp\teqtl_file\tgwas_file\tgene\tn_snps_tested\tclpp_score\tgwas_max_log_pval\teqtl_max_log_pval\n")

    # Loop through all genes now
    for gene in genes:
        combined = preprocess.combine_summary_statistics(gwas_data, eqtl_data, gene, snp, settings)

        # Skip it if this site is untestable.
        if isinstance(combined, basestring):
            print "Skipped gene {}; no SNPs overlapping between eQTL and GWAS.".format(gene)
            continue

        # Create a TestLocus object using merged GWAS and eQTL,
        # any important metadata about the experiment such as the directory,
        # and the Config object.
        task = TestLocus(combined, base_output_dir, base_tmp_dir, gene, snp, settings)
        task.run()

if __name__ == "__main__":
	main()
