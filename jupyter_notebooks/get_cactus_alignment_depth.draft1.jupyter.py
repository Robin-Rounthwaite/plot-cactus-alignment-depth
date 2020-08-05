#%%
from make_bed_intervals import write_intervals

from argparse import ArgumentParser
import subprocess

def get_depths(chrom, full_asm_length, interval_length, intermediate_dir, mapping_bed, outbed):
    if intermediate_dir[:-1] == "/":
        interval_bed = intermediate_dir + "intervals.bed"
    else:
        interval_bed = intermediate_dir + "/intervals.bed"

    write_intervals(chrom, full_asm_length, interval_length, interval_bed)

    with open(outbed, "w") as outf:
        subprocess.run(["coverageBed", "-a", interval_bed, "-b", mapping_bed], stdout=outf)


def main():
    parser = ArgumentParser()
    parser.add_argument('chrom', help='Name of the chromosome being subdivided into intervals.', type=str)
    parser.add_argument('full_asm_length', help='Length of chrom.', type=str)
    parser.add_argument('interval_length', help='Length of the interval.', type=str)
    parser.add_argument('mapping_bed', help='Name of the mapping bedfile which indicates the coverage on the chrom of interest.', type=str)
    parser.add_argument('outbed', help='Name of the output file.', type=str)
    parser.add_argument('--intermediate_dir', help='Name of the output file.', default="intermediate_files/", type=str)
    options = parser.parse_args()

    #test:
    options.chrom = "chr20"
    options.full_asm_length = 64444167
    options.interval_length = 100000
    options.mapping_bed = "test/CHM13_paf_chr20.sorted.liftover.bed"
    options.outbed = "test/CHM13_paf_chr20.100k_interval.mapping_depth.bed"
    options.outbed = "test/CHM13_paf_chr20_intermediate_files/"

    get_depths(options.chrom, int(options.full_asm_length), int(options.interval_length), options.intermediate_dir, options.mapping_bed, options.outbed)

    # write_intervals(options.chrom, int(options.full_asm_length), int(options.interval_length), options.intermediate_dir, options.mapping_bed, options.outbed)

if __name__ == "__main__":
    main()