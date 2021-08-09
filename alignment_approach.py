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
    def __init__(self):
        # self.min_match_size = min_match_size
        # self.MAX_KSIZE = max_ksize
        # self.MIN_KSIZE = min_ksize
        self.input_nmt_reads = sorted(os.listdir(hardcoded_path.numt_reads))
        self.kmer_masked_reads = []
        self.out_ksize = []
    
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
    def parse_name_masked_reads(self, ksize,name):
        out_name = f".maskedN.k{ksize}.numt.reads.fastq"
        s = name.split('.')
        s = s[0] + out_name
        self.out_ksize.append(ksize)
        return s

    def write_fastq_files(self,ksize):
    # min_match_size = 31
        for idx, name in enumerate(self.input_nmt_reads):
            sub_folder = name.split(".")[0]
            # while(ksize <= self.MAX_KSIZE):
            # for ksize in tqdm(range(self.MIN_KSIZE, self.MAX_KSIZE+2, 2)):
            numt_bam = pysam.AlignmentFile(os.path.join(hardcoded_path.numt_reads, name), 'rb')
            masked_reads = self.parse_name_masked_reads(ksize, name)
            if not os.path.exists(os.path.join(hardcoded_path.kmer_fastq_files, sub_folder)):
                os.mkdir(os.path.join(hardcoded_path.kmer_fastq_files, sub_folder))
            if not os.path.exists(os.path.join(hardcoded_path.Minimap_output, sub_folder)):
                os.makedirs(os.path.join(hardcoded_path.Minimap_output, sub_folder))
            fastq_file = open(os.path.join(hardcoded_path.kmer_fastq_files, sub_folder, masked_reads), 'w')
            no_sup = set()
            masked_query_seq = ""
            # not_mapped_reads = [] 
            c = 0
            for i in numt_bam:
                if (i.is_read1 == True) and (i.query_name + '/1' not in no_sup):
                    masked_query_seq = self.mask_matched_bases(i.get_blocks(),ksize,i.get_aligned_pairs(matches_only = True), i.query_sequence)
                    runs_of_ATGC = re.findall(r"[ATGC]{35,}", masked_query_seq)
                    if len(runs_of_ATGC) > 0:
                        no_sup.add('@' + i.query_name + '/1' + '\n' + masked_query_seq + '\n' + '+' + '\n' + 'I'*len(masked_query_seq) + '\n')
                elif (i.is_read2 == True) and (i.query_name + '/2' not in no_sup):
                    # if i.cigarstring != "100M":
                    masked_query_seq = self.mask_matched_bases(i.get_blocks(),ksize,i.get_aligned_pairs(matches_only = True), i.query_sequence)
                    runs_of_ATGC = re.findall(r"[ATGC]{35,}", masked_query_seq)
                    if len(runs_of_ATGC) > 0:
                        no_sup.add('@' + i.query_name + '/1' + '\n' + masked_query_seq + '\n' + '+' + '\n' + 'I'*len(masked_query_seq) + '\n')
            fastq_reads = []
            for i in no_sup:
                fastq_reads.append(i)
                fastq_file.write(i)
            fastq_file.close()
            
                
        self.kmer_masked_reads = sorted(os.listdir(hardcoded_path.kmer_fastq_files))


    def parse_name_assembled_reads(self,name,ksize):
        out_name = f".spades.minimap.sorted.bam"
        # out_list = []
        # dirs = self.kmer_masked_reads
        idx = 0
        print(name)
        s = name.split(".")
        print(s)
        s = s[0] + "." + str(ksize) + out_name
        return s

    def assemble_masked_reads(self):
        dirs = self.kmer_masked_reads
        
        idx = 0
        for sub_dir in dirs:
            s_dir1 = os.listdir(os.path.join(hardcoded_path.kmer_fastq_files,sub_dir))
            s_dir1 = sorted(s_dir1)
            for file in s_dir1:
                print(self.out_ksize)
                assembled = self.parse_name_assembled_reads(file,self.out_ksize[idx])
                in_path = os.path.join(hardcoded_path.kmer_fastq_files, sub_dir, file)
                out_path = os.path.join(hardcoded_path.Minimap_output, sub_dir, assembled)
                print(in_path, out_path)
                idx += 1
                spades_cmd = f"/home/simarpreet/SPAdes-3.12.0-Linux/bin/spades.py --cov-cutoff 2 -s {in_path} --phred-offset 33 -o /home/simarpreet/Desktop/Remote_pc_files/spades_output"
                s = subprocess.run([spades_cmd], shell = True)
                minimap_cmd = f"minimap2 -a /home/simarpreet/Desktop/Remote_pc_files/ref.mmi --split-prefix temp_file /home/simarpreet/Desktop/Remote_pc_files/spades_output/contigs.fasta > {out_path} "
                m = subprocess.run([minimap_cmd], shell= True)


# nm.append(i.get_tag("NM"))
# plt.hist(nm, bins = 7)
# plt.savefig("Alignment_approach_dist_of_nm.png")
    
# aln_app = alignment_approach(19,17)
# aln_app.write_fastq_files()
# aln_app.assemble_masked_reads()