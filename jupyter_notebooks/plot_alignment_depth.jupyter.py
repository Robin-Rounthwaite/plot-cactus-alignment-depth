#%%
# read in alignment depth file (output of coverageBed)
import pandas as pd
import matplotlib.pyplot as plt
depth_file = "test/CHM13_paf_chr20.100k_interval.mapping_depth.bed"

alignment_depth = pd.read_csv(depth_file, sep="\t", names=["chrom", "start", "stop", "num_reads", "bases_covered", "bases_total", "coverage_ratio"])


#%%
print(len(alignment_depth))
alignment_depth.plot.bar(y="num_reads")
plt.show()
#%%
ax = alignment_depth.plot.line(x="start", y="num_reads")
ax.set_xlabel("chromsome region (binned into 100k subregions)")
ax.set_ylabel("number of assembly mappings", )
ax.set_title("mappings in cactus between reference and assembly X", )



# from get_alignment_depth import get_alignment_depth

# def plot_fig():

# def main():
#     get_alignment_depth()
#     plot_fig()
    
# if __name__ == "__main__":
#     main()