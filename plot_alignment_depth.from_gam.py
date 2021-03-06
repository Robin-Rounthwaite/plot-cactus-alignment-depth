#%%
import os
import subprocess
import matplotlib.pyplot as plt
from types import SimpleNamespace
import pandas as pd
from argparse import ArgumentParser
from plot_alignment_depth_shared_fxns import *

def get_mapping_bed(gam, xg, chrom_name, intermediate_dir, outBed, vg_dir):
    bam = intermediate_dir + "/" + (".".join(gam.split(".")[:-1]) + ".bam").split("/")[-1]
    with open(bam, "w") as outf:
        subprocess.run([vg_dir, "surject", "-b", "-p", chrom_name, "-x", xg, gam], stdout=outf)
    unsorted_bed = intermediate_dir + "/" + (".".join(outBed.split(".")[:-1]) + ".unsorted.bam").split("/")[-1]
    with open(unsorted_bed, "w") as outf:
        subprocess.run(["bamToBed", "-i", bam], stdout=outf)
    with open(outBed, "w") as outf:
        subprocess.run(["sortBed", "-i", unsorted_bed], stdout=outf)
    return outBed

    # return unsortedBed

def plot_gam_depths(gam, xg, xg_chrom, chrom, reads_name, full_chrom_length, interval_length, intermediate_dir, output_file, options, ylim=(0,1.1), vg_dir="/home/robin/paten_lab/vg_binary"):
    """
    Current goal: plot depths for first 500Kb of chr20.
    """
    bed = intermediate_dir + "/" + (".".join(gam.split(".")[:-1]) + ".bed").split("/")[-1]
    bed = get_mapping_bed(gam, xg, xg_chrom, intermediate_dir, bed, vg_dir)

    fig, ax = plt.subplots()

    plot_depth(ax, bed, reads_name, chrom, interval_length, full_chrom_length, reads_name, intermediate_dir, options)

    ax.set_xlabel("Chromosome location (binned into " + human_format(interval_length) + " subregions)")
    ax.set_ylabel("Average mapping depth in each " + human_format(interval_length) + " subregion", )
    ax.set_title("Mapping depth of " + reads_name + " on " + chrom)
    if ylim is not None:
        bottom, top = plt.ylim()
        ax.set_ylim(ylim)

    # plt.savefig(output_dir + "/" + reads_name + "_mapped_to_" + chrom + "." + human_format(interval_length) + "_interval.mapping_depth.png")
    plt.savefig(output_file)
    
def debug_gam_depths(test_bed, chrom, reads_name, full_chrom_length, interval_length, intermediate_dir, output_file, options, ylim=(0,1.1), vg_dir="/home/robin/paten_lab/vg_binary"):
    fig, ax = plt.subplots()
    
    plot_depth(ax, test_bed, reads_name, chrom, interval_length, full_chrom_length, reads_name, intermediate_dir, options)

    ax.set_xlabel("Chromosome location (binned into " + human_format(interval_length) + " subregions)")
    ax.set_ylabel("Average mapping depth in each " + human_format(interval_length) + " subregion", )
    ax.set_title("Mapping depth of " + reads_name + " on " + chrom)
    if ylim is not None:
        # bottom, top = plt.ylim()
        ax.set_ylim(ylim)

    # plt.savefig(output_dir + "/" + reads_name + "_mapped_to_" + chrom + "." + human_format(interval_length) + "_interval.mapping_depth.png")
    plt.savefig(output_file)


def main():
    # parser = ArgumentParser()
    # parser.add_argument('chrom', help='Name of the chromosome being subdivided into intervals.', type=str)
    # parser.add_argument('full_asm_length', help='Length of chrom.', type=int)
    # parser.add_argument('interval_length', help='Length of the interval.', type=int)
    # parser.add_argument('mapping_bed', help='Name of the mapping bedfile which indicates the coverage on the chrom of interest.', type=str)
    # parser.add_argument('depth_file', help='Name of the output file.', type=str)
    # parser.add_argument('--intermediate_dir', help='Name of the output file.', default="intermediate_files/", type=str)
    # options = parser.parse_args()


    ## chr20 test:
    # for running on Robin's computer:
    # python plot_alignment_depth.from_gam.py /home/robin/paten_lab/cactus_projects/analyze_chr20_cactus_alignments/mapping_pileups/first_100.gam /home/robin/paten_lab/cactus_projects/analyze_chr20_cactus_alignments/mapping_pileups/lc2019_12ont-hg38.cactus.minimap2_star-all-to-ref-fatanc-no-secondary-july-8.xg hg38_chr20.chr20 HG002-glenn_giabfeb26 64444167 1000000 chr20_test/plot_alignment_depth.from_gam/intermediate_files chr20_test/plot_alignment_depth.from_gam/HG002-glenn_giabfeb26.giraffe39k15wPaths.mapped_to_chr20.1M_interval.mapping_depths.png
    # for debugging on Robin's computer:
    # python plot_alignment_depth.from_gam.py hg38_chr20.chr20 HG002-glenn_giabfeb26 64444167 32222000 chr20_test/debug_gam_depths/intermediate_files chr20_test/debug_gam_depths/debug_gam_depths.png --run_test --test_bed debug_gam_depths.bed 
    # simplified:
    # python plot_alignment_depth.from_gam.py gam.gam xg.xg hg38_chr20.chr20 HG002-glenn_giabfeb26 64444167 1000000 intermediate_files HG002-glenn_giabfeb26.giraffe39k15wPaths.mapped_to_chr20.1M_interval.mapping_depths.png

    parser = ArgumentParser()
    #todo: fix it so I don't have to  comment out gam + xg when running --run_test
    parser.add_argument("gam", help='', type=str)
    parser.add_argument("xg", help='', type=str)
    parser.add_argument("xg_chrom", help='', type=str)
    parser.add_argument("reads_name", help='', type=str)
    parser.add_argument("full_chrom_length", help='', type=int) 
    parser.add_argument("interval_length", help='', type=int)
    parser.add_argument("intermediate_dir", help='', type=str)
    parser.add_argument("output_file", help='', type=str)
    parser.add_argument("--vg_dir", help='', default='vg', type=str)
    parser.add_argument("--run_test", help='', action='store_true')
    parser.add_argument("--test_bed", help='', type=str)
    options = parser.parse_args()

    options.chrom = options.xg_chrom

    # ## chr20 test:
    # options = SimpleNamespace()
    # options.gam = "/home/robin/paten_lab/cactus_projects/analyze_chr20_cactus_alignments/mapping_pileups/first_100.gam"
    # options.xg = "/home/robin/paten_lab/cactus_projects/analyze_chr20_cactus_alignments/mapping_pileups/lc2019_12ont-hg38.cactus.minimap2_star-all-to-ref-fatanc-no-secondary-july-8.xg"
    # options.xg_chrom = "hg38_chr20.chr20"
    # options.chrom = options.xg_chrom
    # # options.chrom = "chr20"
    # options.reads_name = "HG002-glenn_giabfeb26"
    # options.full_chrom_length = 64444167
    # options.interval_length = 1000000 # reasonable number of bins
    # options.intermediate_dir = "chr20_test/plot_alignment_depth.from_gam/intermediate_files"
    # options.output_file = "chr20_test/plot_alignment_depth.from_gam/"

    ## cleaning input dirs:
    #ensure intermediate_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.intermediate_dir) and os.path.exists(options.intermediate_dir):
        raise ValueError("--intermediate_dir is a file, not a directory.")
    elif not os.path.isdir(options.intermediate_dir) and not os.path.exists(options.intermediate_dir):
        os.makedirs(options.intermediate_dir)

    if options.intermediate_dir[-1] == "/":
        options.intermediate_dir = options.intermediate_dir[:-1]

    if options.run_test:
        debug_gam_depths(options.test_bed, options.chrom, options.reads_name, options.full_chrom_length, options.interval_length, options.intermediate_dir, options.output_file, options, ylim=None, vg_dir=options.vg_dir)
    else:
        plot_gam_depths(options.gam, options.xg, options.xg_chrom, options.chrom, options.reads_name, options.full_chrom_length, options.interval_length, options.intermediate_dir, options.output_file, options, ylim=None, vg_dir=options.vg_dir)
if __name__ == "__main__":
    main()