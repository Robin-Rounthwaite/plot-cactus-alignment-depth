seqdir=small_chr21
echo "seqdir: $seqdir"

#todo: just realized that output is in terms of chr21. As in, we find out how much of chr21 is covered by the asm.
#todo:  This means I *did* do the halLiftovers correctly originally - I need to flip those arguments back, and also pass contig_lengths for ref instead of the asm.

#indexing ref:
samtools faidx $seqdir/hg38_chr21.fa -o $seqdir/hg38_chr21.fa.fai
for prefix in HG03098_paf_chr21 HG03492_paf_chr21
do
    # to run dipcall:
    ./run-dipcall $prefix $seqdir/hg38_chr21.fa $seqdir/$prefix.fa $seqdir/$prefix.fa > $prefix.mak
    make -f $prefix.mak

    #non bedmap approach to extracting coverage info:
    echo "coverage from dipcall from file $prefix.fa"
    bedops --merge $prefix.dip.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done

#for merging any non-overlapped bed intervals, then getting the length of the intervals: 
# bedops --merge A | bedmap --echo --echo-map-size hg19.extents.bed - > answer.bed

#for getting extents of file "prefix"
# cat $seqdir/$prefix.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t0\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $prefix.extents.bed 

    #bedmap approach to extracting bed coverage info: NOTE: confusing bedmap output with semicolon list. Also doesn't seem to see any of the intervals?
    # cat $seqdir/hg38_chr21.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t0\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > hg38_chr21.extents.bed 
    # bedops --merge $prefix.dip.bed | bedmap --echo --echo-map-size hg38_chr21.extents.bed - > $prefix.total_coverage.bed
