import vcf
import pysam
from intervaltree_bio import GenomeIntervalTree
import re

##To make a bed file of the chrM variants
## chrM 73 73 G/A This means that my sample has a G instead of an A at position 73 of its mitochondria
# vcf_reader = vcf.Reader(open('ERR1019039_chrM_1000x.vcf.gz', 'rb'))

##we convert from 1-based, closed [start, end] Variant Call Format v4 (VCF) to sorted, 0-based, half-open [start-1, end).
# with open('chrM_variants_1000x.bed', 'w') as f:
#     print('CHROM', 'Start', 'End', 'REF/ALT', file=f)
#     for record in vcf_reader:
#         print(record.CHROM,record.POS - 1 ,record.POS ,str(record.REF) + "/" + str(record.ALT[0]), file = f)

##Check if my reads overlap this position of the chrM:
##If they do we want to ignore those reads.
mt_var = []
chrM_variants =  GenomeIntervalTree()
for line in open("chrM_variants_1000x.bed"):
    fields = line.rstrip().split(" ")
    if fields[0].startswith('chrM'):
        # print(fields)
        chrM_variants[fields[0]].addi(int(fields[1]), int(fields[2]), fields[3])
        mt_var.append(int(fields[1]))

def variant_overlap(chrm, begin, end):
    var_ovlp = chrM_variants[chrm].overlaps(begin,end)
    ovlp_position  = chrM_variants[chrm].search(begin,end)
    return var_ovlp, ovlp_position

bamfile = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.numt_tag_new.bam', 'rb')
# outfile = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.numt_tag_no_chrM.bam', 'wb', header = bamfile.header)

def get_index(MD_tag,reference_start, Mt_snp_pos):

    split_MD = re.split('(\d+)', MD_tag)
    clean_MD = [x.strip() for x in split_MD if x.strip()]
    alpnum = [int(x) if x.isdigit() else x for x in clean_MD]
    # print(alpnum)
    # if int(clean_MD[0]) + reference_start == Mt_snp_pos:
    #     read_snp_pos = reference_start + int(clean_MD[0])
    #     read_snp = int(clean_MD[1])
        
    #     print(read_snp)
    # # else:
        # a = []
        # read_snp_pos = reference_start + 

    # return read_snp_pos, read_snp
    return alpnum

noise = []
remaining_reads = []

# import pysam
# samfile = pysam.Samfile(bam, "rb")
# for pileupcolumn in samfile.pileup( '1', 14693, 14694):
#     for pileupread in pileupcolumn.pileups:
#         if pileupcolumn.pos == 14693:
#             base = pileupread.alignment.query_sequence[pileupread.query_position]
#             print base
c = 0
f = []
for pileupcolumn in bamfile.pileup('chrM', 0, 16568):
    # if i in mt_var:
    if pileupcolumn.reference_pos in mt_var:
        print(pileupcolumn.get_num_aligned(), pileupcolumn.reference_pos)
    # print(pileupcolumn.get_query_positions())
    # print()
        for pileupread in pileupcolumn.pileups:
            # print(pileupcolumn.reference_pos, pileupread.query_position, pileupcolumn.nsegments,pileupread.alignment.query_sequence[pileupread.query_position])
            # print(pileupread.alignment)
            # if pileupread.alignment.query_sequence[pileupread.query_position] != 'G':
            print(pileupread.alignment.query_sequence[pileupread.query_position])
            ###look for read and mate here pileupread.alignment is pysam.AlignedSegment object 
            # remove both read and mate if anyone overlaps and contains ALT allele
            ##is_reverse == True ?? look at that or see if pileup has taken care of it already

        c += 1
        if c == 1:
            break
    #     print(pileupread.alignment)
    #     print(pileupread.alignment.query_sequence[pileupread.query_position])
# print(len(f))