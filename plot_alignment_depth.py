#%%
import os
import subprocess
import matplotlib.pyplot as plt
from types import SimpleNamespace
import pandas as pd

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

def subset_file(to_subset, chrom, start, stop, intermediate_dir, outfile):
    if not os.path.exists(intermediate_dir + "/for_subsetting_inputs/"):
        os.mkdir(intermediate_dir + "/for_subsetting_inputs/")

    #make the subsetter file
    chrom_region_file = intermediate_dir + "/for_subsetting_inputs/" + chrom + "_" + str(start) + "_" + str(stop) + ".bed"
    with open(chrom_region_file, "w") as region:
        region.write(chrom + "\t" + str(start) + "\t" + str(stop)) 
        
    #subset the file
    outfile_name = outfile.split("/")[-1]
    unsorted_outfile = intermediate_dir + "/for_subsetting_inputs/" + ".".join(outfile_name.split(".")[:-1]) + ".unsorted.bed"
    with open(unsorted_outfile, "w") as outf:
        subprocess.run(["intersectBed", "-a", chrom_region_file, "-b", to_subset], stdout=outf)

    #sort the subset file:
    with open(outfile, "w") as outf:
        subprocess.run(["sortBed", "-i", unsorted_outfile], stdout=outf)

    return outfile

def get_bin_depths(subset_bedfile, chrom, bin_start, bin_stop, intermediate_dir):
    if not os.path.exists(intermediate_dir + "/for_get_bin_depths/"):
        os.mkdir(intermediate_dir + "/for_get_bin_depths/")
    
    genome_file = intermediate_dir + "/for_get_bin_depths/" + chrom + "_to_" + str(bin_stop) + ".genome.txt"
    with open(genome_file, "w") as outf:
        outf.write(chrom + "\t" + str(bin_stop))

    #subset_name = name of file without directory info or file extension:
    subset_name = ".".join((subset_bedfile.split("/")[-1]).split(".")[:-1])
    depths_file = intermediate_dir + "/for_get_bin_depths/" + subset_name + ".depths.txt"
    with open(depths_file, "w") as outf:
        subprocess.run(["genomeCoverageBed", "-i", subset_bedfile, "-g", genome_file], stdout=outf)

    depths = pd.read_csv(depths_file, "\t", names=["chrom", "depth", "num_bases", "chrom_size", "coverage_proportion"])

    #to adjust depths output so that depths only calculates within bin_start to bin_end:
    for i in range(len(depths)):
        depths.loc[i, 'chrom_size'] -= bin_start
        if depths.loc[i, 'depth'] == 0:
            depths.loc[i, 'num_bases'] -= bin_start
        depths.loc[i, 'coverage_proportion'] = depths.loc[i, 'num_bases']/depths.loc[i, 'chrom_size']

    return depths

def get_avg_depth(depths, chrom):
    depths_for_chrom = depths[depths["chrom"] == chrom]
    chrom_size = depths_for_chrom.loc[0, "chrom_size"]
    avg_depth = sum(depths_for_chrom["depth"].to_numpy() * depths_for_chrom["num_bases"].to_numpy()) / chrom_size
    return avg_depth
        
def plot_depth(ax, bedfile, asm, chrom, interval_length, full_asm_length, legend_label, intermediate_dir,):
    if not os.path.exists(intermediate_dir + "/subsets/"):
        os.mkdir(intermediate_dir + "/subsets/")
    
    bedfile_name = bedfile.split("/")[-1]
    bedfile_root_name = ".".join(bedfile_name.split(".")[:-1])

    #calculate all bin_start, bin_stop bin ranges for the plot.
    bins = list()
    bin_start = int()
    while bin_start != full_asm_length:
        if bin_start + interval_length <= full_asm_length:
            bins.append((bin_start, bin_start+interval_length))
            bin_start += interval_length
        else:
            bins.append((bin_start, full_asm_length))
            bin_start = full_asm_length

    x = [plt_bin[0] for plt_bin in bins]
    y = list() # to be filled with avg_depth values
    for (bin_start, bin_stop) in bins:
        #extract the region of bedfile dictated by bin_start and bin_stop
        subset = subset_file(bedfile, chrom, bin_start, bin_stop, intermediate_dir, intermediate_dir+"/subsets/"+bedfile_root_name+"_subset_"+str(bin_start)+"_"+str(bin_stop)+".bed")

        depths = get_bin_depths(subset, chrom, bin_start, bin_stop, intermediate_dir)

        y.append(get_avg_depth(depths, chrom))
        print("bin", bin_start, "finished for", bedfile)

    ax.plot(x, y, label=legend_label)
    ax.legend()        

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
        os.mkdir(options.dipcall_dir)
    #ensure liftover_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.liftover_dir) and os.path.exists(options.liftover_dir):
        raise ValueError("--liftover_dir is a file, not a directory.")
    elif not os.path.isdir(options.liftover_dir) and not os.path.exists(options.liftover_dir):
        os.mkdir(options.liftover_dir)
    #ensure output_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.output_dir) and os.path.exists(options.output_dir):
        raise ValueError("--output_dir is a file, not a directory.")
    elif not os.path.isdir(options.output_dir) and not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    #ensure intermediate_dir exists and hasn't a "/" at the end:
    if not os.path.isdir(options.intermediate_dir) and os.path.exists(options.intermediate_dir):
        raise ValueError("--intermediate_dir is a file, not a directory.")
    elif not os.path.isdir(options.intermediate_dir) and not os.path.exists(options.intermediate_dir):
        os.mkdir(options.intermediate_dir)

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