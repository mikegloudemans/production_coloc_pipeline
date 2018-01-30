# Author: Mike Gloudemans
#
# Class: TestLocus
# A locus to test for colocalization, along with all of the
# relevant data needed to conduct the tests.
#

import finemap
import plot_loci as plot
import math

class TestLocus:
    def __init__(self, data, basedir, tmpdir, gene, snp, settings):
        self.basedir = basedir
        self.data = data
        self.settings = settings
        self.gene = gene
        self.chrom = snp.chrom
        self.pos = snp.pos
        self.pval = snp.pval
        self.gwas_file = settings["gwas_file"]
        self.eqtl_file = settings["eqtl_file"]
        self.gwas_suffix = self.gwas_file.split("/")[-1].replace(".", "_")
        self.eqtl_suffix = self.eqtl_file.split("/")[-1].replace(".", "_")
        self.tmpdir = tmpdir
        self.conditional_level = 0      # Currently serves no purpose, but may be implemented later
        self.vcf_file = settings["vcf_file"]
        self.N = settings["N"]

    # Run colocalization tests with FINEMAP
    def run(self):

        print "Analyzing {0} {1} {2} {3} {4} {5}".format(self.gwas_suffix, self.eqtl_suffix, self.gene, self.chrom, self.pos, self.pval)

        clpp = finemap.run_finemap(self)

        if clpp == "Fail":
            print "FAILED Analyzing {0} {1} {2} {3} {4} {5}".format(self.gwas_suffix, self.eqtl_suffix, self.gene, self.chrom, self.pos, self.pval)

        finemap.purge_tmp_files(self)

        # Plot results
        plot.locus_zoom_plot(self, clpp)
        plot.pvalue_plot(self, clpp)


