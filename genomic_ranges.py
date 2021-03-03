# import pybedtools
# a = pybedtools.BedTool('ReferenceNumts.bed')

# for numt in a:
#     print(numt)

# import pyranges as pr
# import os
# from interval import interval

# path = os.path.abspath("ReferenceNumts.bed")
# numt_df = pr.read_bed(path, as_df = True)
# print(numt_df)
#print(numt_df[['Start','End']])


import HTSeq
bed_file = HTSeq.BED_Reader("Reference_Numts_edited.bed")
alignment_file = HTSeq.BAM_Reader('ERR1019039.MT.csort.alignment.bam')

# f = open("Reference_Numts_edited.bed", mode = 'w')

# for line in open( "ReferenceNumts.bed" ):
#     fields = line.split( "\t" )
#     if 'Dayama Numt\n' in fields[3]:
#         fields[3] = 'Dayama_Numt\n'
#     f.write('\t'.join(fields))

####Create genomic interval of the mate using alignment file information:
for ref_numt in bed_file:
    for aln in alignment_file:
        for cigop in aln.cigar:
            try:
                read_iv = HTSeq.GenomicInterval( cigop.ref_iv.chrom, cigop.ref_iv.start + 1 , cigop.ref_iv.end + 1 , cigop.ref_iv.strand )
            except Exception:
                pass        
        try:
            mate_cigop = HTSeq.parse_cigar(aln.optional_field("MC"), aln.mate_start.pos + 1, aln.mate_start.chrom, aln.mate_start.strand)
        except Exception:
            pass 
        #print(mate_cigop)
        for i in mate_cigop:
            try:
                mate_iv = HTSeq.GenomicInterval( aln.mate_start.chrom, aln.mate_start.pos + 1 , i.ref_iv.end , aln.mate_start.strand )
            except Exception:
                pass
            if cigop.ref_iv.is_contained_in(ref_numt.iv) == True or mate_iv.is_contained_in(ref_numt.iv) == True:
                #print("read_iv", cigop.ref_iv, ref_numt.iv, "mate_iv",mate_iv, ref_numt.iv)
                print(aln.get_sam_line())
            elif cigop.ref_iv.overlaps(ref_numt.iv) == True or mate_iv.overlaps(ref_numt.iv) == True:
                #print("read_iv", cigop.ref_iv, ref_numt.iv, "mate_iv",mate_iv, ref_numt.iv)
                print(aln.get_sam_line())
            