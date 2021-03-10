import HTSeq
bed_file = HTSeq.BED_Reader("Reference_Numts_edited.bed")
alignment_file = HTSeq.BAM_Reader('ERR1019039.MT.csort.alignment.bam')

for ref_numt in bed_file:
    for aln in alignment_file:
        for cigop in aln.cigar:
            #print(cigop.ref_iv)
            # print(cigop)
            try:
                read_iv = HTSeq.GenomicInterval( cigop.ref_iv.chrom, cigop.ref_iv.start + 1 , cigop.ref_iv.end + 1 , cigop.ref_iv.strand )
            except Exception:
                pass
            #print(read_iv)       
        try:
            mate_cigop = HTSeq.parse_cigar(aln.optional_field("MC"), aln.mate_start.pos + 1, aln.mate_start.chrom, aln.mate_start.strand)
        except Exception:
            pass 
        for i in mate_cigop:
            try:
                mate_iv = HTSeq.GenomicInterval( aln.mate_start.chrom, aln.mate_start.pos + 1 , i.ref_iv.end , aln.mate_start.strand )
            except Exception:
                pass
            try:
                if read_iv.is_contained_in(ref_numt.iv) == True:
                    if mate_iv.is_contained_in(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MCN," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    elif mate_iv.overlaps(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MON," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    elif aln.mate_start.chrom == "chrM":
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MAM," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    else:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MNO," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")

                elif read_iv.overlaps(ref_numt.iv) == True:
                    if mate_iv.is_contained_in(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MCN," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    elif mate_iv.overlaps(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MON," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    elif aln.mate_start.chrom == "chrM":
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MAM," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    else:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MNO," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")

                elif aln.iv.chrom == "chrM":
                    if mate_iv.is_contained_in(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MCN," + str(-1) + ";")
                    elif mate_iv.overlaps(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MON," + str(-1) + ";")
                    elif aln.mate_start.chrom == "chrM":
                        print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MAM," + str(-1) + ";")
                    else:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MNO," + str(-1) + ";")
                
                else:
                    if mate_iv.is_contained_in(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SNO,MCN," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    elif mate_iv.overlaps(ref_numt.iv) == True:
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SNO,MON," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                    elif aln.mate_start.chrom == "chrM":
                        print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SNO,MAM," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
            except Exception:
                pass     
                
               