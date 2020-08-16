#%%
import os
import subprocess
import matplotlib.pyplot as plt
from types import SimpleNamespace
import pandas as pd
from plot_alignment_depth_shared_fxns import *

def plot_depths_for_all_asms(all_asm_liftovers, all_asm_dipcalls, chrom, full_asm_length, interval_length, intermediate_dir, output_dir, ylim=(0,1.1)):
    """
    Current goal: plot depths for first 500Kb of chr20.
    """
    #todo: make sure that I only pass bin_stop that's less than full_asm_length
    for asm in all_asm_liftovers:
        fig, ax = plt.subplots()

        plot_depth(ax, all_asm_dipcalls[asm], asm, chrom, interval_length, full_asm_length, "dipcall", intermediate_dir)
        plot_depth(ax, all_asm_liftovers[asm], asm, chrom, interval_length, full_asm_length, "cactus", intermediate_dir)

        ax.set_xlabel("Chromosome location (binned into " + human_format(interval_length) + " subregions)")
        ax.set_ylabel("Average mapping depth in each " + human_format(interval_length) + " subregion", )
        ax.set_title("Mapping depth in cactus vs. dipcall for " + asm)
        bottom, top = plt.ylim()
        ax.set_ylim(ylim)

        plt.savefig(output_dir + "/" + asm + "_mapped_to_" + chrom + "." + human_format(interval_length) + "_interval.mapping_depth.png")
    

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
    options = SimpleNamespace()
    options.liftover_dir = "../cactus-connectivity/chr20_bed_overlaps/connectivity_beds.sorted/"
    options.dipcall_dir = "../cactus-connectivity/chr20_bed_overlaps/dipbeds.sorted/"
    options.chrom = "chr20"
    options.full_asm_length = 64444167
    options.interval_length = 1000000 # reasonable number of bins
    # options.interval_length = 30000000
    # options.interval_length = 500000
    # #for plotting coverage:
    # options.intermediate_dir = "chr20_test/chr20_ref_mapping_depth_output/intermediate_files"
    # options.output_dir = "chr20_test/chr20_ref_mapping_depth_output/"
    # options.comparison_type = "coverage"
    #for plotting depth (less useful for comparison with dipcall):
    options.intermediate_dir = "chr20_test/chr20_ref_mapping_depth_output_using_genome_cov/intermediate_files"
    options.output_dir = "chr20_test/chr20_ref_mapping_depth_output_using_genome_cov/"
    options.comparison_type = "depth"

    ## cleaning input dirs:
    #ensure dipcall_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.dipcall_dir) and os.path.exists(options.dipcall_dir):
        raise ValueError("--dipcall_dir is a file, not a directory.")
    elif not os.path.isdir(options.dipcall_dir) and not os.path.exists(options.dipcall_dir):
        os.makedirs(options.dipcall_dir)
    #ensure liftover_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.liftover_dir) and os.path.exists(options.liftover_dir):
        raise ValueError("--liftover_dir is a file, not a directory.")
    elif not os.path.isdir(options.liftover_dir) and not os.path.exists(options.liftover_dir):
        os.makedirs(options.liftover_dir)
    #ensure output_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.output_dir) and os.path.exists(options.output_dir):
        raise ValueError("--output_dir is a file, not a directory.")
    elif not os.path.isdir(options.output_dir) and not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    #ensure intermediate_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.intermediate_dir) and os.path.exists(options.intermediate_dir):
        raise ValueError("--intermediate_dir is a file, not a directory.")
    elif not os.path.isdir(options.intermediate_dir) and not os.path.exists(options.intermediate_dir):
        os.makedirs(options.intermediate_dir)

    if options.dipcall_dir[-1] == "/":
        options.dipcall_dir = options.dipcall_dir[:-1]
    if options.liftover_dir[-1] == "/":
        options.liftover_dir = options.liftover_dir[:-1]
    if options.output_dir[-1] == "/":
        options.output_dir = options.output_dir[:-1]
    if options.intermediate_dir[-1] == "/":
        options.intermediate_dir = options.intermediate_dir[:-1]

    ## finding locations of liftovers and dipcalls in dir:
    all_asm_liftovers = dict()
    liftover_names = os.listdir(options.liftover_dir)
    for name in liftover_names:
        asm = name.split(".")[0] # requires that asm be the first component of filename
        all_asm_liftovers[asm] = options.liftover_dir + "/" + name

    all_asm_dipcalls = dict()
    dipcall_names = os.listdir(options.dipcall_dir)
    for name in dipcall_names:
        asm = name.split(".")[0] # requires that asm be the first component of filename
        all_asm_dipcalls[asm] = options.dipcall_dir + "/" + name

    if options.comparison_type == "coverage":
        plot_mappings_for_all_asms(all_asm_liftovers, all_asm_dipcalls, options.chrom, options.full_asm_length, options.interval_length, options.intermediate_dir, options.output_dir)
    elif options.comparison_type == "depth":
        plot_depths_for_all_asms(all_asm_liftovers, all_asm_dipcalls, options.chrom, options.full_asm_length, options.interval_length, options.intermediate_dir, options.output_dir)
    else:
        raise ValueError("comparison_type must be 'depth' or 'coverage'")
    
if __name__ == "__main__":
    main()

#%%
output_dir = "chr20_test/chr20_ref_mapping_depth_output_using_genome_cov/"

print(output_dir)