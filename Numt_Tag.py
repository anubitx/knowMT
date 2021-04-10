import pysam
import os
import simplesam
from intervaltree_bio import GenomeIntervalTree

##Next step pysam pileup function to identify mitochondrial mutations

# path_bed = os.path.abspath("Reference_Numts_with_chrM.bed")
# gr_bed = pr.read_bed(path_bed)
# print(gr_bed)

## Checking for the dna fragment instead of individual reads for overlap

##Conducted a test to cross check if the operation that consumes the query is accounted for or not. 
##It is accounted for, hence we use simplesam package.
# y = simplesam.Sam(cigar = '5M5I5N5=')
# print(len(y))

##Next step pysam pileup function to identify mitochondrial mutations

ref_numt = GenomeIntervalTree()

for line in open("Reference_Numts_with_chrM.bed"):
    fields = line.rstrip().split("\t")
    ref_numt[fields[0]].addi(int(fields[1]), int(fields[2]), (fields[3], int(fields[2]) - int(fields[1])))

def numt_Overlap(chr, start, end):
    overlap = 0
    overlapid = "NA"
    overlapcontained = 0
    overlap_count = 0
    nlen = -1
    refoverlap = ref_numt[chr].search(start,end)

    if len(refoverlap) > 0:
        overlap = 1
        overlap_count = len(refoverlap)
        refoverlap = refoverlap.pop()
        overlapid = refoverlap[2][0]
        nlen = refoverlap[2][1]
        if refoverlap[0] <= start and refoverlap[1] >= end:
            overlapcontained = 1
        elif refoverlap[0] >= start and refoverlap[1] <= end:
            overlapcontained = 2
    return str(overlap), str(overlapcontained), overlapid, str(nlen) , str(overlap_count)

 
samfile = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.bam', 'rb')
outfile = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.numt_tag.bam', 'wb', header = samfile.header)

for read in samfile:
    if (not(read.is_unmapped)):
        self_info = numt_Overlap(read.reference_name, read.reference_start, read.reference_end)
    else:
        self_info = ("0", "0", "NA", "-1" , "0")
    if (not( read.mate_is_unmapped)):
        x = simplesam.Sam(cigar = read.get_tag("MC"))
        mate_info = numt_Overlap(read.next_reference_name, read.next_reference_start, read.next_reference_start + len(x))
    else:
       mate_info = ("0", "0", "NA", "-1" , "0")
    if self_info[0] == "1" or mate_info[0] == "1":
        read.set_tag("Ni", ','.join(self_info + mate_info), value_type = "Z", replace = False)
        outfile.write(read)

# gr_read = pr.PyRanges(chromosomes = [read.reference_name], starts= [read.reference_start], ends = [read.reference_end])
# gr_mate = pr.PyRanges(chromosomes = [read.next_reference_name], starts= [read.next_reference_start], ends = [read.next_reference_start + len(x)])
    
