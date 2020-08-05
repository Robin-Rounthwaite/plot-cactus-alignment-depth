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


def get_depths(chrom, full_asm_length, interval_length, intermediate_dir, mapping_bed, depth_file):
    if intermediate_dir[-1] == "/":
        interval_bed = intermediate_dir + "intervals.bed"
    else:
        interval_bed = intermediate_dir + "/intervals.bed"

    write_intervals(chrom, full_asm_length, interval_length, interval_bed)

    with open(depth_file, "w+") as outf:
        subprocess.run(["coverageBed", "-a", interval_bed, "-b", mapping_bed], stdout=outf)

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

def main():
    # parser = ArgumentParser()
    # parser.add_argument('chrom', help='Name of the chromosome being subdivided into intervals.', type=str)
    # parser.add_argument('full_asm_length', help='Length of chrom.', type=str)
    # parser.add_argument('interval_length', help='Length of the interval.', type=str)
    # parser.add_argument('mapping_bed', help='Name of the mapping bedfile which indicates the coverage on the chrom of interest.', type=str)
    # parser.add_argument('depth_file', help='Name of the output file.', type=str)
    # parser.add_argument('--intermediate_dir', help='Name of the output file.', default="intermediate_files/", type=str)
    # options = parser.parse_args()

    #test:
    options = SimpleNamespace()
    options.chrom = "chr20"
    options.full_asm_length = 64444167
    options.interval_length = 500000
    # options.interval_length = 100000 #What I use for plotting with num_reads as line in plot_depths 
    options.liftover_bed = "CHM13_test/CHM13_paf_chr20.sorted.liftover.bed"
    options.dipcall_bed = "CHM13_test/CHM13_paf_chr20.sorted.dip.bed"
    options.depth_file = "CHM13_test/CHM13_paf_chr20.100k_interval.mapping_depth.tsv"
    options.intermediate_dir = "CHM13_test/CHM13_paf_chr20_intermediate_files/"
    options.assembly_name = "CHM13_paf_chr20"
    if not os.path.isdir(options.intermediate_dir) and os.path.exists(options.intermediate_dir):
        raise ValueError("--intermediate_dir is a file, not a directory.")
    elif not os.path.isdir(options.intermediate_dir) and not os.path.exists(options.intermediate_dir):
        os.mkdir(options.intermediate_dir)

    fig, ax = plt.subplots()
    
    # for plotting Mapping Coverage (GOOD FOR DIPCALL):
    # plot the dipcall line
    get_depths(options.chrom, int(options.full_asm_length), int(options.interval_length), options.intermediate_dir, options.dipcall_bed, options.depth_file)
    plot_mapping_coverage_by_proportion_seq_mapped(ax, options.depth_file, options.assembly_name, "dipcall", options.interval_length)

    # plot the liftover line
    get_depths(options.chrom, int(options.full_asm_length), int(options.interval_length), options.intermediate_dir, options.liftover_bed, options.depth_file)
    plot_mapping_coverage_by_proportion_seq_mapped(ax, options.depth_file, options.assembly_name, "cactus", options.interval_length)

    # # for plotting Mapping Coverage (unclear comparison with dipcall):
    # # plot the dipcall line
    # get_depths(options.chrom, int(options.full_asm_length), int(options.interval_length), options.intermediate_dir, options.dipcall_bed, options.depth_file)
    # plot_mapping_coverage_by_mapping_depth(ax, options.depth_file, options.assembly_name, "dipcall", options.interval_length)

    # # plot the liftover line
    # get_depths(options.chrom, int(options.full_asm_length), int(options.interval_length), options.intermediate_dir, options.liftover_bed, options.depth_file)
    # plot_mapping_coverage_by_mapping_depth(ax, options.depth_file, options.assembly_name, "cactus", options.interval_length)


    plt.savefig(options.assembly_name+".coverage.png")



if __name__ == "__main__":
    main()