# FINEMAP colocalization pipeline

## Author: Mike Gloudemans

### About

This pipeline fine-maps eQTL and GWAS using [FINEMAP](http://www.christianbenner.com/), an efficient tool developed by Benner et al.

The colocalization test is performed by multiplying posterior causal probabilities 
from the eQTL and GWAS studies, to generate a colocalization posterior probability (CLPP),  
as described in the 2016 [eCAVIAR paper](https://www.ncbi.nlm.nih.gov/pubmed/27866706) by Hormozdiari et al.

The input files for this pipeline are a GWAS file, an eQTL file, and a reference VCF file (for computing LD;
we recommend using a VCF from 1000 Genomes).

### Setup

To run this pipeline, first clone this repository into your desired directory.

You will need to install [PLINK 1.9](https://www.cog-genomics.org/plink2) and [FINEMAP](http://www.christianbenner.com/).
Binary executables for these tools must be placed in the `/bin` directory. The PLINK executable must be named `plink`
and the FINEMAP executable must be named `finemap`.

### Preparing Input Files

#### GWAS and eQTL summary statistics

GWAS and eQTL input files must be tab-delimited lists of SNPs, containing at least the following columns:

* `chr`: Chromosome
* `snp_pos`: Position of SNP in hg19
* `alt`: Effect allele
* `ref`: Non-effect allele
* `beta`: Estimated effect size of the `alt` allele, when compared with the `ref` allele
* `se`: Standard error of `beta`

eQTL files must additionally contain a `gene` column that specifies which gene's expression
is being tested.

The column order does not matter, but they _must_ be specified in the header of the file.

An example GWAS input file:
```
chr     snp_pos alt     ref     beta        se          pvalue
1       751756  C       T       .013006     .017324     .4528019
1       752566  A       G       -.005243    .0157652    .7394597
1       752721  G       A       -.003032    .0156381    .8462652
1       752894  C       T       .00464      .0162377    .7750657
[...]
```

Note that it is okay to have additional columns beyond the required ones; these will be ignored when the program is run.
What's important is that all the required columns are included.

(Note: For the eQTL files, a chi-squared test statistic can be specified under a `chisq`
column as an alternative to the `beta` and `se` columns.)

The GWAS, eQTL, and VCF files must be sorted first by chromosome and then by snp position,
and finally zipped with `bgzip` and indexed using `tabix`. This ensures that the pipeline will run
at a reasonable speed.

The following two commands could be used to index a file named `my_eqtl_file.txt` with the chromosome number
in the second column (corresponds to `-s` parameter) and the chromosome position in the third column (`-b` and `-e` parameters):

```
bgzip -f /home/my_eqtl_file.txt
tabix -f -S 1 -s 2 -b 3 -e 3 /home/my_eqtl_file.txt.gz
```

#### Reference VCF

FINEMAP requires a reference VCF to compute the LD between positions of interest.

We recommend using the 1000 Genomes VCF, which can be obtained freely from the 1000 Genomes
website.


### Running the Pipeline

To run the pipeline, you must specify the input files, the reference genome, and the number
of individuals in the reference genome. You must also specify a SNP position of interest; the pipeline
will analyze all SNPs near this reference SNP. Currently, SNPs on the sex chromosomes are not supported.

By default, the pipeline will test for colocalization with every gene that has eQTLs in the region.
To limit analysis to a single gene, modify the `gene` parameter as described below in the "optional
parameters" section.

The required parameters are `chrom`, `snp_pos`, `eqtl_file`, `gwas_file`, `vcf_file`,  and `N`.

* `chrom`: chromosome of reference SNP of interest
* `snp_pos`: position of reference SNP of interest
* `eqtl_file`: bgzipped and tabixed eQTL summary statistics file
* `gwas_file`: bgzipped and tabixed GWAS summary statistics file
* `vcf_file`: bgzipped and tabixed reference genome VCF, for computing LD
* `N`: number of individuals in the reference genome

The `vcf_file` parameter can be modified to accommodate separate files for each chromosome.
i.e. If the name of the VCF for chromosome is `chr1.vcf.gz`, then you can input `chr{0}.vcf.gz`
and the script will substitute the appropriate chromosome number as necessary.



_Example_:

The following command would test for colocalization between `gwas.txt.gz` and eQTLs
for all genes in `eqtl.txt.gz` that are near the position chr1:10000000. LD between SNPs
would be computed using `1000genomes.chr{0}.vcf.gz`, which includes 2504 individuals.

```
python ./dispatch.py -g gwas.txt.gz -e eqtls.txt.gz -v 1000genomes.chr{0}.vcf.gz \
            -N 2504 -c 1 -s 10000000
```

_Optional parameters_:

* `gene`: limit analysis to this gene only
* `eqtl_format`: can be one of three values: [`effect_size`, `chisq`, `tstat`]. The default value is effect_size, which will be used for most eQTL studies. 
    eQTLs called using RASQUAL may instead report chi-squared statistics, which requires use of the `chisq` option. In this case, the eQTL summary statistics file
    must have a `chisq` column but does not require `beta` and `se`.
* `gwas_format`: can be one of three values: [`effect_size`, `log_odds_ratio`, `pval_only`]. The default value is `effect_size`. `log_odds_ratio` allows the `beta`
        column to be replaced with a column labeled `log_or`. `pval_only` requires a `pvalue` column; with this option the `beta` and `se` columns may be omitted.
* `window`: number of bp to analyze in either direction from the reference SNP. Default value is 500000. Large values (>1MB) may take significantly longer to run.

### Viewing Results

The results of the pipeline are placed in a time-stamped folder in the output directory.

At the top level, the directory contains a file `[gwas_name]_finemap_clpp_status.txt`.
This file contains a row for each locus tested. The sixth column is the CLPP score.
For more information on interpreting the CLPP score, see [Hormozdiari et al](https://www.ncbi.nlm.nih.gov/pubmed/27866706).
We suggest filtering tested loci by the CLPP score and then visualizing the locus for further confirmation, to avoid
false positives.

A good starting filter is to select loci with at least 100 variants tested, a CLPP score of > 0.02,
and a top -log(pvalue) of > 5 for both the eQTL and the GWAS studies. These filters can be made more
stringent as needed if too many false positives are being detected.

"Rainbow" visualizations of the tested loci are available in the `plots` subdirectory. Each locus
contains one scatterplot showing eQTL significance vs. GWAS significance, and a double Manhattan plot
showing the significance of each SNP for the GWAS study on top and the eQTL study
on the bottom. Each SNP is colored the same in all three of the plots. However, the colors provide only
information about physical distance in base pairs; no information is given about LD. To visualize
LD, we recommend instead using a tool like LocusZoom.

The `finemap` directory contains two file for each colocalization test, showing the posterior causal
probabilities for each SNP at the locus. These probabilities are shown in the `snp_prob` column.
