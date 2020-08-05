# plot-cactus-alignment-depth
Plotting was originally built for comparing cactus and dipcall alignments between the reference and a set of given assemblies. Edit and run plot_alignment_depth.py to make plots. 

Goal of plot_alignment_depth:
For each subregion "bin" of the reference, plot the average number of times that subregion
is aligned to a chosen assembly. Scale for every possible asm.

Input:
    * liftover file, generated srcGenome [asm], srcBed [entire expanse of asm], tgtGenome [ref].
    This means that the liftover bed is in terms of reference coordinates.

    * dipcall file, between the same asm and ref.

    * number of bins to subdivide the reference into

    * total length of reference

Output:
    * Bar plot of average mapping depth for each subregion of the reference.
