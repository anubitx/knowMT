import pysam
import HTSeq
import pandas as pd

#import collections
bed_file = HTSeq.BED_Reader("Reference_Numts_edited.bed")
alignment_file = HTSeq.BAM_Reader('ERR1019039.MT.csort.alignment.bam')
bam_writer = HTSeq.BAM_Writer.from_BAM_Reader("Potential_Numts_With_Added_Tags_1.bam", alignment_file)
#f = open("PN.sam", "a")
#f1 = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.bam', 'wb')
#bam_writer = HTSeq.BAM_Writer("Potential_Numts_With_Added_Tags_1.bam")


#SCN_MCN_count = 0 
max_count = [] 
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
        for i in mate_cigop:
            try:
                mate_iv = HTSeq.GenomicInterval( aln.mate_start.chrom, aln.mate_start.pos + 1 , i.ref_iv.end , aln.mate_start.strand )
            except Exception:
                pass
            try:

                if read_iv.is_contained_in(ref_numt.iv) == True:
                    if mate_iv.is_contained_in(ref_numt.iv) == True:
                        SCN_MCN = aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MCN," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";"
                        data = SCN_MCN.split("\t")
                        print(data)
                        #df = pd.DataFrame(SCN_MCN.split("\t")).transpose()
                        #print(df[12])
                        #print(df.loc['ERR1019039.621090469'])
                        #print(data[0], data[-1])
                        # df = pd.DataFrame(SCN_MCN.split("\t")).transpose()
                        #a = df.set_index(0)
                        #print(a)
                        # print(a[16])
                        #print(a[15].fillna(value=0, inplace=True))
                        #a.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'MD', 'AS', 'XS', 'Ni' ]
                        #print(max_count[:100])
                        #print(pd.DataFrame(data=data,columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'MD', 'MC', 'AS', 'XS', 'SA','XA', 'Ni' ]).transpose())
                        # print(pd.DataFrame(data, ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'MD', 'MC', 'AS', 'XS', 'XA', 'Ni' ])
                        #print(SCN_MCN.split("\t")[-1])
                        #bam_writer.write(SCN_MCN)
                        # f.write(SCN_MCN)
                    # elif mate_iv.overlaps(ref_numt.iv) == True:
                    #      SCN_MON = aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MON," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";"
                    # elif aln.mate_start.chrom == "chrM":
                    #      SCN_MAM = aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MAM," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";"
                    # else:
                    #      SCN_MNO = aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SCN,MNO," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";"

                #  elif read_iv.overlaps(ref_numt.iv) == True:
                #      if mate_iv.is_contained_in(ref_numt.iv) == True:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MCN," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                #      elif mate_iv.overlaps(ref_numt.iv) == True:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MON," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                #      elif aln.mate_start.chrom == "chrM":
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MAM," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                #      else:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SON,MNO," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")

                #  elif aln.iv.chrom == "chrM":
                #      if mate_iv.is_contained_in(ref_numt.iv) == True:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MCN," + str(-1) + ";")
                #      elif mate_iv.overlaps(ref_numt.iv) == True:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MON," + str(-1) + ";")
                #      elif aln.mate_start.chrom == "chrM":
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MAM," + str(-1) + ";")
                #      else:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:chrM,SAM,MNO," + str(-1) + ";")
                
                #  else:
                #      if mate_iv.is_contained_in(ref_numt.iv) == True:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SNO,MCN," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                #      elif mate_iv.overlaps(ref_numt.iv) == True:
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SNO,MON," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")
                #      elif aln.mate_start.chrom == "chrM":
                #          print(aln.get_sam_line() + "\t" + "Ni:Z:" + ref_numt.name + ",SNO,MAM," + str(ref_numt.iv.end - ref_numt.iv.start - 1) + ";")


            except Exception:
                pass