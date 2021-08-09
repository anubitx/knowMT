from posix import ST_RELATIME
import pysam
from collections import defaultdict
import statistics
import matplotlib.pyplot as plt
import math
import hardcoded_path
import os
from tqdm import tqdm
import subprocess

class remove_mitochondrial_noise():
    # def __init__(self, Step3_output_dirpath:str) -> None:
    def __init__(self) -> None:

        # self.Step3_output_dirpath = Step3_output_dirpath
        # self.Step3_output_dirpath = hardcoded_path.step3_output_dir
        self.out_file = sorted(os.listdir(hardcoded_path.out_dir))
        self.step3_out = self.parse_bam_name_for_step3_outfile()
        # self.step3_output_filename = self.get_step3_output()
        # self.step3_output_filename = sorted(os.listdir(self.Step3_output_dirpath))
        self.step3_output_filename = sorted(os.listdir(hardcoded_path.step3_output_dir))
        self.counts_dict = defaultdict(int)
        self.stepf_outfile_name = self.parse_name_for_Stepf_outfile()
        self.paired = set()
        self.get_input = self.get_input_for_insert()
        self.final_mt = self.parse_name_for_MT_final_outfile()

        
        self.numt_reads = self.parse_name_numt_reads()
        # print(os.path.join(hardcoded_path.step3_output_dir, self.con_Mt_filename))

    def parse_bam_name_for_step3_outfile(self):
        out_name = ".Genome.numt_tag.consensus_Mt.bam"
        outfile_list = []
        for bam_name in self.out_file:
            s = bam_name.split('.')
            s = f"{s[0]}.{s[1]}{out_name}"
            outfile_list.append(s)
        return outfile_list

    def step3_for_all_files(self):
        cwd = os.getcwd()
        p = os.path.join(cwd, "Step3.sh")
        # if not os.path.exists(hardcoded_path.step3_output_dir):
            # os.mkdir(hardcoded_path.step3_output_dir)
        #If step3_output_dir doesnot exist create it
        for idx, file_name in enumerate(self.out_file):
            in_path = os.path.join(hardcoded_path.out_dir, file_name)
            out_path = os.path.join(hardcoded_path.step3_output_dir, self.step3_out[idx])
            print(in_path,out_path)
            cmd = f" -o={out_path} -i={in_path}"
            cmd = p+cmd
            print(cmd)
            s = subprocess.Popen([cmd], shell = True)
            out, err = s.communicate()

    # def get_step3_output(self):
        # return sorted(os.listdir(self.Step3_output_dirpath))

    def parse_name_for_Stepf_outfile(self):
        ## pp = proper paired and fl = full length
        out_name = "pp_fl_temp"
        out_list = []
        for name in self.step3_output_filename:
            s = name.split('.')
            s.insert(-1, out_name)
            out_list.append('.'.join(s))
        # print(out_list)
        return out_list

    def getreads_mapped_in_full_length(self,sam_file, out_file):
        # for idx,i in enumerate(self.step3_output_filename):
        #     sam_file = pysam.AlignmentFile(os.path.join(self.Step3_output_dirpath,i), 'rb')
        #     out_file = pysam.AlignmentFile(os.path.join(hardcoded_path.stepf_outfile, self.stepf_outfile_name[idx]), 'wb', header = sam_file.header)
        for i in sam_file:
            if (not(i.is_unmapped)):
                if i.is_proper_pair:
                    if i.get_tag('NM') == 0 and i.get_tag('MD') == '100':
                        self.counts_dict[i.query_name] += 1
                        # print(i)
                        ##How do I write this out_file to stepf_outfile folder?
                        out_file.write(i)

    def get_counts_from_dict(self):
        for k,v in self.counts_dict.items():
            if v == 2:
                self.paired.add(k)        
            
    def get_input_for_insert(self):
        return os.listdir(hardcoded_path.stepf_outfile)

    def parse_name_numt_reads(self):
        out_name = "only_numt_reads"
        out_list = []
        subdir = sorted(os.listdir(hardcoded_path.bam_dir))
        for name in subdir:
            s = name.split('.')
            s.insert(-1, out_name)
            out_list.append('.'.join(s))
        # print(out_list)
        return out_list

    def get_insert_stats(self):
        for idx, name in enumerate(tqdm(self.step3_output_filename)):
            insert = []
            self.paired = set()
            self.counts_dict = defaultdict(int)
            # sam_file = pysam.AlignmentFile(os.path.join(self.Step3_output_dirpath,self.step3_output_filename[idx]), 'rb')
            sam_file = pysam.AlignmentFile(os.path.join(hardcoded_path.step3_output_dir,name), 'rb')
            out_file =  pysam.AlignmentFile(os.path.join(hardcoded_path.stepf_outfile, self.stepf_outfile_name[idx]), 'wb', header = sam_file.header)
            self.getreads_mapped_in_full_length(sam_file,out_file)
            out_file.close()
            sam_file.close()
            in_bam = pysam.AlignmentFile(os.path.join(hardcoded_path.stepf_outfile,self.stepf_outfile_name[idx]), 'rb')
            self.get_counts_from_dict()
            for j in in_bam:
                if j.query_name in self.paired:
                    if j.flag == 99 or j.flag == 83:
                        insert_size = abs(j.template_length) - 200 
                        insert.append(insert_size)
            # print(insert)
            mean_insert_size = statistics.mean(insert)
            SD_insert_size = statistics.stdev(insert)
            pvar_insert_size =  statistics.pvariance(insert, mu= mean_insert_size)
            SD_99_low =  mean_insert_size - 3 * SD_insert_size 
            SD_99_upp = mean_insert_size + 3 * SD_insert_size
            plt.hist(insert)
            # plt.hist(insert, bins = 832)
            plt.savefig(f"Insert Size Distribution No bins{idx}.png")
            plt.close()
            print(f"Mean Insert Size:{mean_insert_size},SD Insert Size:{SD_insert_size},Population var Insert size:{pvar_insert_size}")

            Mitochondrial_reads_only =  pysam.AlignmentFile(os.path.join(hardcoded_path.Mt_reads_dir, self.final_mt[idx]), 'wb', header = in_bam.header)
            in_bam = pysam.AlignmentFile(os.path.join(hardcoded_path.stepf_outfile,self.stepf_outfile_name[idx]), 'rb')
            mt_origin_set = self.paired_end_reasonable_distance_reads([SD_99_low,SD_99_upp],Mitochondrial_reads_only ,in_bam)
            print("Total number of reads detected as Mitochondrial in origin:", len(mt_origin_set) * 2 )
            # Mitochondrial_reads_only.close()
            ##Instead of opening numt_file (outfile directory), open the file in Step3_output directory (MT.Genome.numt_tag.consensus_Mt.bam)
            ##and write only the reads that are not in mt_origin_set to the folder Only_numt_reads 
            ##ERR1019039.MT.csort.alignment.only_numt_reads.bam
            # numt_file = pysam.AlignmentFile(os.path.join(hardcoded_path.out_dir, self.out_file[idx]), 'rb')
            # sam_file = pysam.AlignmentFile(os.path.join(self.Step3_output_dirpath,name), 'rb')
            sam_file = pysam.AlignmentFile(os.path.join(hardcoded_path.step3_output_dir,name), 'rb')
            numt_reads = pysam.AlignmentFile(os.path.join(hardcoded_path.numt_reads,self.numt_reads[idx]), 'wb', header = sam_file.header)
            candidate_numt_reads = []
            for n in sam_file:
                # print(n)
                if n.query_name not in mt_origin_set:
                    candidate_numt_reads.append(n.query_name)
                    numt_reads.write(n)
            print("Total number of candidate numt reads:", len(candidate_numt_reads))
            numt_reads.close()

    def parse_name_for_MT_final_outfile(self):
        ## pp = proper paired and fl = full length
        out_name = "Mt_origin_reads"
        out_list = []
        for name in self.step3_output_filename:
            s = name.split('.')
            s.insert(-1, out_name)
            out_list.append('.'.join(s))
        return out_list

    def paired_end_reasonable_distance_reads(self, SD_range, Mitochondrial_reads_only, in_bam):
        mt_origin_qn = set()
        print(round(SD_range[0]), round(SD_range[1]))
        for j in in_bam:
            if j.query_name in self.paired:
                insert_size = abs(j.template_length) - 200
                if insert_size in range(round(SD_range[0]), round(SD_range[1])):
                    mt_origin_qn.add(j.query_name)
                    Mitochondrial_reads_only.write(j)
                    
        return mt_origin_qn

c = remove_mitochondrial_noise()
# c.getreads_mapped_in_full_length()
c.step3_for_all_files()
c.get_insert_stats()


        
# Mitochondrial_reads_only.close()  
# Mitochonrial_reads_only = pysam.AlignmentFile('ERR1019039.MT.Ref.nmt_tag.con_MT.Mt_origin_reads.bam', 'rb')

# mt_origin_qn = set()
# for i in Mitochonrial_reads_only:
#     mt_origin_qn.add(i.query_name)

# numt_file = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.numt_tag_new.bam', 'rb')
# numt_reads = pysam.AlignmentFile('ERR1019039.MT.csort.alignment.numt_reads.bam', 'wb', header = numt_file.header)

# for j in numt_file:
    # if j.query_name not in mt_origin_qn:
        # numt_reads.write(j)

#original_files/outfile/ERR1019039.MT.csort.alignment.numt_tag.bam


##eg when you have to call the class:
# 
# in_bam = pysam.AlignmentFile(os.path.join(hardcoded_path.stepf_outfile, os.listdir(hardcoded_path.stepf_outfile)[0]), 'rb')