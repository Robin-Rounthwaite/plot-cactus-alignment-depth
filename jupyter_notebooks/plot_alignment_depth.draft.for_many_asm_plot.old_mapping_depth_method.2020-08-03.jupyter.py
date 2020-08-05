#%%
#todo: make it so that the plot shows comparison line between cactus and dipcall
#todo: make it so that I can run get_alignment_depth on a directory of assemblies. (Do it in the python code - this would be a common use of this code).

from argparse import ArgumentParser
import subprocess
from types import SimpleNamespace
import os
import pandas as pd
import matplotlib.pyplot as plt
def human_format(num):
    """
    rtaft's answer from https://stackoverflow.com/questions/579310/formatting-long-numbers-as-strings-in-python
    """
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])

    
def write_intervals(chrom, full_asm_length, interval_length, interval_bed):
    with open(interval_bed, "w+") as outf:
        cur_pos = int()
        while cur_pos != full_asm_length:
            #todo: make sure no off-by-one errors here. Bed is 0-based. What about my length calculation? If it's one-based, need to subtract one for compatibility.
            if cur_pos+interval_length > full_asm_length:
                outf.write(chrom + "\t" + str(cur_pos) + "\t" + str(full_asm_length) + "\n")
                cur_pos = full_asm_length
            else:
                outf.write(chrom + "\t" + str(cur_pos) + "\t" + str(cur_pos + interval_length) + "\n")
                cur_pos += interval_length

    return interval_bed


def write_depths(interval_bed, mapping_bed, depth_file):
    with open(depth_file, "w+") as outf:
        subprocess.run(["coverageBed", "-a", interval_bed, "-b", mapping_bed], stdout=outf)
    return depth_file

def plot_mapping_coverage_by_proportion_seq_mapped(ax, depth_file, assembly_name, legend_label, interval_length):
    """
    Most useful for comparing reference sequence coverage between two graphs. This is a good metric for comparing to dipcall.
    """
    alignment_depth = pd.read_csv(depth_file, sep="\t", names=["chrom", "start", "stop", "num_reads", "bases_covered", "bases_total", "coverage_ratio"])
    alignment_depth.plot.line(x="start", y="coverage_ratio", ax=ax, label=legend_label)
    ax.set_xlabel("Chromosome location (binned into " + human_format(interval_length) + " subregions)")
    ax.set_ylabel("Proportion of ref seq covered by mapping\nin each " + human_format(interval_length) + " subregion", )
    ax.set_title("Mapping coverage in cactus vs. dipcall for " + assembly_name, )

    # ax.set_yscale('log')
    return ax

def plot_mapping_coverage_by_mapping_depth(ax, depth_file, assembly_name, legend_label, interval_length):
    """
    Most useful for comparing between two graphs, to see how many distinct mapping events are overlapping a given region.
    """
    alignment_depth = pd.read_csv(depth_file, sep="\t", names=["chrom", "start", "stop", "num_reads", "bases_covered", "bases_total", "coverage_ratio"])
    alignment_depth.plot.line(x="start", y="num_reads", ax=ax, label=legend_label)
    ax.set_xlabel("Chromosome location (binned into " + human_format(interval_length) + " subregions)")
    ax.set_ylabel("Number of mapping events between asm and ref\nin each " + human_format(interval_length) + " subregion", )
    ax.set_title("Mapping depth in cactus vs. dipcall for " + assembly_name, )

    ax.set_yscale('log')
    return ax

def plot_mappings_for_all_asms(all_asm_liftovers, all_asm_dipcalls, chrom, full_asm_length, interval_length, intermediate_dir, output_dir, comparison_type="coverage"):
    """
    Compares liftover mappings and dipcall mappings to a given reference chrom.
    #Note: Could replace dipcall mappings with mapppings from a different graph or asm to 
    #compare mapping coverage between those.
    
    comparison_type options: "depth" or "coverage". "coverage" is default, as it's the
    most useful metric for comparison to dipcall.

    #todo: If we ever work with a reference with >1 chrom, generalize function? Could make subplots for each chrom...
    """
    if intermediate_dir[-1] == "/":
        intermediate_dir = intermediate_dir[:-1]

    interval_bed = intermediate_dir + "/" + chrom + ".subdivided_to_" + human_format(interval_length) + "_intervals.bed"
    write_intervals(chrom, full_asm_length, interval_length, interval_bed)

    for asm in all_asm_dipcalls.keys():
        dipcall_depth_file = intermediate_dir + "/" + asm + "_mapped_to_" + chrom + ".dipcall." + human_format(interval_length) + "_interval.tsv"
        write_depths(interval_bed, all_asm_dipcalls[asm], dipcall_depth_file)

        liftover_depth_file = intermediate_dir + "/" + asm + "_mapped_to_" + chrom + ".cactus_liftover." + human_format(interval_length) + "_interval.tsv"
        write_depths(interval_bed, all_asm_liftovers[asm], liftover_depth_file)

        fig, ax = plt.subplots()
        if comparison_type == "coverage":
            # for plotting Mapping Coverage (GOOD FOR DIPCALL):
            plot_mapping_coverage_by_proportion_seq_mapped(ax, dipcall_depth_file, asm, "dipcall", interval_length)
            plot_mapping_coverage_by_proportion_seq_mapped(ax, liftover_depth_file, asm, "cactus", interval_length)
        elif comparison_type == "depth":
            # for plotting Mapping Coverage (unclear comparison with dipcall):
            plot_mapping_coverage_by_mapping_depth(ax, dipcall_depth_file, asm, "dipcall", interval_length)
            plot_mapping_coverage_by_mapping_depth(ax, liftover_depth_file, asm, "cactus", interval_length)
        else:
            raise ValueError("comparison_type must be 'depth' or 'coverage'")

        if output_dir[-1] == "/":
            output_dir = output_dir[:-1]
        plt.savefig(output_dir + "/" + asm + "_mapped_to_" + chrom + "." + human_format(interval_length) + "_interval.mapping_" + comparison_type + ".png")
        

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
    options.interval_length = 500000
    #for plotting coverage:
    options.intermediate_dir = "chr20_test/chr20_ref_mapping_depth_output/intermediate_files"
    options.output_dir = "chr20_test/chr20_ref_mapping_depth_output/"
    options.comparison_type = "coverage"
    # #for plotting depth (less useful for comparison with dipcall):
    # options.intermediate_dir = "chr20_test/chr20_ref_mapping_depth_output/intermediate_files"
    # options.output_dir = "chr20_test/chr20_ref_mapping_depth_output/"
    # options.comparison_type = "depth"


    if options.liftover_dir[-1] != "/":
        options.liftover_dir += "/"
    if options.dipcall_dir[-1] != "/":
        options.dipcall_dir += "/"

    all_asm_liftovers = dict()
    liftover_names = os.listdir(options.liftover_dir)
    for name in liftover_names:
        asm = name.split(".")[0] # requires that asm be the first component of filename
        all_asm_liftovers[asm] = options.liftover_dir + name

    all_asm_dipcalls = dict()
    dipcall_names = os.listdir(options.dipcall_dir)
    for name in dipcall_names:
        asm = name.split(".")[0] # requires that asm be the first component of filename
        all_asm_dipcalls[asm] = options.dipcall_dir + name

    #ensure output_dir exists and has a "/" at the end:
    if not os.path.isdir(options.output_dir) and os.path.exists(options.output_dir):
        raise ValueError("--intermediate_dir is a file, not a directory.")
    elif not os.path.isdir(options.output_dir) and not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    if options.output_dir[-1] != "/":
        options.output_dir += "/"
    #ensure intermediate_dir exists and has a "/" at the end:
    if not os.path.isdir(options.intermediate_dir) and os.path.exists(options.intermediate_dir):
        raise ValueError("--intermediate_dir is a file, not a directory.")
    elif not os.path.isdir(options.intermediate_dir) and not os.path.exists(options.intermediate_dir):
        os.mkdir(options.intermediate_dir)
    if options.intermediate_dir[-1] != "/":
        options.intermediate_dir += "/"

    plot_mappings_for_all_asms(all_asm_liftovers, all_asm_dipcalls, options.chrom, options.full_asm_length, options.interval_length, options.intermediate_dir, options.output_dir, comparison_type=options.comparison_type)


if __name__ == "__main__":
    main()