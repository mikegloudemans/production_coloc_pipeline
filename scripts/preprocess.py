# Author: Mike Gloudemans
#
# preprocess.py
#
# Tools for loading summary statistics.
#

import subprocess 
import pandas as pd 
import operator 
import SNP 
from scipy import stats 
import math
import gzip

import sys 
if sys.version_info[0] < 3: 
    from StringIO import StringIO 
else: 
    from io import StringIO


# Load summary statistics for GWAS
def get_gwas_data(settings, snp, window=500000):

    gwas_file = settings["gwas_file"]

    # Get GWAS data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(gwas_file), shell=True)

    raw_gwas = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True) + \
            subprocess.check_output("tabix {0} chr{1}:{2}-{3}".format(gwas_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    gwas_table = pd.read_csv(StringIO(header + raw_gwas), sep="\t")
    gwas_table['snp_pos'] = gwas_table['snp_pos'].astype(int)
    
    gwas_table['ref'] = gwas_table['ref'].apply(lambda x: x.upper())
    gwas_table['alt'] = gwas_table['alt'].apply(lambda x: x.upper())

    #
    # 'gwas_format' must be specified, to make sure the users know what they're doing.
    #
    # Possible settings for 'gwas_format':
    #   - case_control
    #   - effect_size
    #   - pval_only (requires effect direction)
    #

    if settings['gwas_format'] == 'log_odds_ratio':
        assert 'log_or' in gwas_table and 'se' in gwas_table
        gwas_table['ZSCORE'] = gwas_table['log_or'] / gwas_table['se']
    elif settings['gwas_format'] == 'effect_size':
        assert 'beta' in gwas_table
        gwas_table['ZSCORE'] = gwas_table['beta'] / gwas_table['se']
    elif settings['gwas_format'] == 'pval_only':
        assert 'pvalue' in gwas_table and "direction" in gwas_table
        # Need to cap it at z-score of 40 for outrageous p-values (like with AMD / RPE stuff)
        gwas_table['ZSCORE'] = pd.Series([min(x, 40) for x in stats.norm.isf(gwas_table["pvalue"] / 2)]) * (2*(gwas_table["direction"] == "+")-1)
    else:
        return "Improper GWAS format specification"

    return gwas_table

# Load summary statistics for eQTL
def get_eqtl_data(settings, snp, window=500000):
    
    eqtl_file = settings["eqtl_file"]

    # Get eQTL data using tabix
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 1".format(eqtl_file), shell=True)
    raw_eqtls = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(eqtl_file, \
            snp.chrom, snp.pos - window, snp.pos + window), shell=True)
    eqtls = pd.read_csv(StringIO(header + raw_eqtls), sep="\t")

    if eqtls.shape[0] == 0:
        return "Gene desert."

    eqtls['ref'] = eqtls['ref'].apply(lambda x: x.upper())
    eqtls['alt'] = eqtls['alt'].apply(lambda x: x.upper())

    eqtls['snp_pos'] = eqtls['snp_pos'].astype(int)
    
    #
    # 'eqtl_format' must be specified, to make sure the users know what they're doing.
    #
    # Possible settings for 'eqtl_format':
    #   - tstat
    #   - effect_size
    #   - chisq
    #

    if settings['eqtl_format'] == 'tstat':
        assert 't-stat' in eqtls
        eqtls['ZSCORE'] = eqtls['t-stat']
    elif settings['eqtl_format'] == 'effect_size':
        assert 'beta' in eqtls
        eqtls['ZSCORE'] = eqtls['beta'] / eqtls['se']
    elif settings['eqtl_format'] == 'chisq':
        assert "chisq" in eqtls
        # Here we're dealing with RASQUAL data
        # where effect size is given by allelic imbalance percentage pi.
        # Use max function to protect against underflow in chi2 computation
        eqtls['pvalue'] = stats.chi2.sf(eqtls["chisq"],1)
        eqtls['ZSCORE'] = stats.norm.isf(eqtls['pvalue']/2) * (2 * (eqtls["pi"] > 0.5) - 1)
    else:
        return "Improper eQTL format specification"

    return eqtls

# Input: GWAS pandas dataframe, eQTL pandas dataframe, gene name as a string,
#   GWAS SNP as a tuple.
# Returns: a combined table of summary statistics, or None if we need to skip
#   the site due to insufficient data.
def combine_summary_statistics(gwas_data, eqtl_data, gene, snp, settings): 
    
    # Filter SNPs down to the gene of interest.
    eqtl_subset = eqtl_data[eqtl_data['gene'] == gene]

    # Sometimes the GWAS SNP is outside of the range of eQTLs tested for a certain
    # gene, or on the outside fringe of the range. If this is the case, then skip it.
    if snp.pos > max(eqtl_subset['snp_pos']) - 50000 or snp.pos < min(eqtl_subset['snp_pos']) + 50000:
            return "SNP outside range."
 
    # If not explicitly allowing them, remove pvalues with danger
    # of underflow.
    if min(eqtl_subset['pvalue']) < 1e-150:
        print "Adjusted minimum eQTL pvalue to prevent underflow."
        eqtl_subset['pvalue'] = eqtl_subset['pvalue'].apply(lambda x: max(x, 1e-150))

    if min(gwas_data['pvalue']) < 1e-150:
        print "Adjusted minimum GWAS pvalue to prevent underflow."
        gwas_data['pvalue'] = gwas_data['pvalue'].apply(lambda x: max(x, 1e-150))

    combined = pd.merge(gwas_data, eqtl_subset, on="snp_pos", suffixes=("_gwas", "_eqtl"))

    dup_counts = {}
    for pos in combined['snp_pos']:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1

    combined['dup_counts'] = [dup_counts[pos] for pos in combined['snp_pos']]
    combined = combined[combined['dup_counts'] == 1]

    # Check to make sure there are SNPs remaining; if not, just move on
    # to next gene.
    if combined.shape[0] == 0: 
        return "No overlapping SNPs in eQTL and GWAS"

    return combined

