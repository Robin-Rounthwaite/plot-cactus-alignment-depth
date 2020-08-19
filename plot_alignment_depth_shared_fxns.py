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

    #todo: consider restructuring code somehow so I don't have to run this genomeCoverageBed many times? It slows with each iteration.
    # print("start_coverage_bed")
    #subset_name = name of file without directory info or file extension:
    subset_name = ".".join((subset_bedfile.split("/")[-1]).split(".")[:-1])
    depths_file = intermediate_dir + "/for_get_bin_depths/" + subset_name + ".depths.txt"
    with open(depths_file, "w") as outf:
        subprocess.run(["genomeCoverageBed", "-i", subset_bedfile, "-g", genome_file], stdout=outf)
    # print("end_coverage_bed")

    depths = pd.read_csv(depths_file, "\t", names=["chrom", "depth", "num_bases", "chrom_size", "coverage_proportion"])

    #to adjust depths output so that depths only calculates within bin_start to bin_end:
    for i in range(len(depths)):
        depths.loc[i, 'chrom_size'] -= bin_start
        if depths.loc[i, 'depth'] == 0:
            depths.loc[i, 'num_bases'] -= bin_start
        depths.loc[i, 'coverage_proportion'] = depths.loc[i, 'num_bases']/depths.loc[i, 'chrom_size']

    return depths

def get_avg_depth(depths, chrom, options):
    depths_for_chrom = depths[depths["chrom"] == chrom]
    if options.run_test:
        print("depths", depths)
    chrom_size = depths_for_chrom.loc[0, "chrom_size"]
    avg_depth = sum(depths_for_chrom["depth"].to_numpy() * depths_for_chrom["num_bases"].to_numpy()) / chrom_size
    return avg_depth
        
def plot_depth(ax, bedfile, asm, chrom, interval_length, full_asm_length, legend_label, intermediate_dir, options):
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

        y.append(get_avg_depth(depths, chrom, options))
        print("bin", bin_start, "finished for", bedfile)
    if options.run_test:
        print("x", x)
        print("y", y)
    ax.plot(x, y, label=legend_label)
    ax.legend()        
