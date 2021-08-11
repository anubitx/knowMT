import pysam 
import os
import hardcoded_path
import matplotlib.pyplot as plt
import subprocess
import re
from tqdm import tqdm

# Check number of num nuts with only 1 supporting grade in num nuts tag file 


class alignment_approach():
    # def __init__(self, max_ksize, min_ksize):
    def __init__(self, input_bam:str, output_dir:str, threads:int, ksize:int):
        self.input_bam = input_bam
        self.output_dir = output_dir
        self.threads = threads
        self.ksize = ksize
        bambase = os.path.basename(self.input_bam)
        self.bambase = bambase
        self.input_numts_reads = os.path.join(self.output_dir,os.path.splitext(bambase)[0] + ".numts.reads.bam")
        self.fastq_file = os.path.join(self.output_dir,os.path.splitext(bambase)[0] + f".maskedN.k{ksize}.numt.reads.fastq")
        self.minimap_out = os.path.join(self.output_dir,os.path.splitext(bambase)[0] + f".spadesk{ksize}.minimap.sorted.bam")
    
    def replace_sequence(self, query_sequence, res : list):
        # print(res)
        if len(res) == 0:
            return query_sequence
        else:
            for i in res:
                if i is not None:    
                    start = i[0]
                    end = i[1]
                    query_sequence = query_sequence.replace(query_sequence[start:end], len(query_sequence[start:end]) * "N")
            return query_sequence


    def mask_matched_bases(self,all_matched_coordinates: list,ksize: int,aligned_pairs: list, query_sequence:str):
        res = []
        for j in all_matched_coordinates:
            # print(j)
            # query_start = None
            # query_end = None
            if j[1] - j[0] >= ksize:
                # query_coords = self.get_query_coordinates(all_matched_coordinates, aligned_pairs,query_sequence)
                for i in aligned_pairs:
                    if i[1] == j[0]:
                        query_start = i[0]
                #I have subtracted 1 becasue I am using a number not a range.
                    if i[1] == j[1] - 1: 
                        query_end = i[0] + 1
                if query_start != None and query_end != None:
                    # print((query_start, query_end))
                    res.append((query_start, query_end))
        return self.replace_sequence(query_sequence, res)

# why query_end + 1?                    
#For example 
#34S34M32S
#Without + 1
#(34, 67) 67-34 = 33M 
# whereas there are 34 M in the cigarstring so to take care of that (34, 68)
    def write_fastq_files(self):
        numt_bam = pysam.AlignmentFile(self.input_numts_reads, 'rb', threads = self.threads)
        fastq_file = open(self.fastq_file, 'w')
        no_sup = set()
        masked_query_seq = ""

        for i in numt_bam:
            if (i.is_read1 == True) and (i.query_name + '/1' not in no_sup):
                masked_query_seq = self.mask_matched_bases(i.get_blocks(),self.ksize,i.get_aligned_pairs(matches_only = True), i.query_sequence)
                runs_of_ATGC = re.findall(r"[ATGC]{35,}", masked_query_seq)
                if len(runs_of_ATGC) > 0:
                    no_sup.add('@' + i.query_name + '/1' + '\n' + masked_query_seq + '\n' + '+' + '\n' + 'I'*len(masked_query_seq) + '\n')
            elif (i.is_read2 == True) and (i.query_name + '/2' not in no_sup):
                # if i.cigarstring != "100M":
                masked_query_seq = self.mask_matched_bases(i.get_blocks(),self.ksize,i.get_aligned_pairs(matches_only = True), i.query_sequence)
                runs_of_ATGC = re.findall(r"[ATGC]{35,}", masked_query_seq)
                if len(runs_of_ATGC) > 0:
                    no_sup.add('@' + i.query_name + '/1' + '\n' + masked_query_seq + '\n' + '+' + '\n' + 'I'*len(masked_query_seq) + '\n')
        fastq_reads = []
        for i in no_sup:
            fastq_reads.append(i)
            fastq_file.write(i)
        fastq_file.close()

    def assemble_masked_reads(self, ref_mmi, path_to_spades):
        spades_out = os.path.join(os.getcwd(), 'spades_output')
        contigs_out = os.path.join(os.getcwd(), 'spades_output', 'contigs.fasta')
        spades_cmd = f"{path_to_spades} --cov-cutoff 2 -s {self.fastq_file} --phred-offset 33 -o {spades_out}"
        s = subprocess.run([spades_cmd], shell = True)
        minimap_cmd = f"minimap2 -a {ref_mmi} --split-prefix temp_file {contigs_out} > {self.minimap_out} "
        m = subprocess.run([minimap_cmd], shell= True)


# nm.append(i.get_tag("NM"))
# plt.hist(nm, bins = 7)
# plt.savefig("Alignment_approach_dist_of_nm.png")
    
# aln_app = alignment_approach(19,17)
# aln_app.write_fastq_files()
# aln_app.assemble_masked_reads()