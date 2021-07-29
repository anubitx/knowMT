import pysam
import os
import simplesam
import argparse
from intervaltree_bio import GenomeIntervalTree
import subprocess
from tqdm import tqdm

## Checking for the dna fragment instead of individual reads for overlap
## import hardcoded_path

##Conducted a test to cross check if the operation that consumes the query is accounted for or not. 
##It is accounted for, hence we use simplesam package.
# y = simplesam.Sam(cigar = '5M5I5N5=')
# print(len(y))
    
class numt_chrM_overlapping_reads():
    def __init__(self, numtbed:str, inputbam:str, outputdir:str, threads:int, reference:str) -> None:
        
        self.ref_numt = GenomeIntervalTree()
        self.bed = numtbed
        self.inputbam = inputbam
        self.threads = threads
        self.reference = reference
        self.read_bed()
        self.outputdir = outputdir
        # self.bam_file_name = self.get_bam()
        bambase = os.path.basename(self.inputbam)
        self.outputbam = os.path.join(self.outputdir,os.path.splitext(bambase)[0] + ".numt_tag.bam")
        #self.step3_out = self.parse_bam_name_for_step3_outfile()
        self.generate_numt_tag()
        # self.step3_for_all_files()
        # print(self.outfiles)

    # def parse_bam_name_for_outfile(self):
    #     out_name = "numt_tag.bam"
    #     # outfile_list = []
    #     # for bam_name in self.bam_file_name:
    #     s = self.inputbam.split('.')
    #     s[:-1] + out_name
    #     # outfile_list.append('.'.join(s))
    #     # return outfile_list
    #     return '.'.join(s)
    def parse_bam_name_for_step3_outfile(self):
        out_name = ".Genome.numt_tag.consensus_Mt.bam"
        outfile_list = []
        for bam_name in self.outfiles:
            s = bam_name.split('.')
            s = f"{s[0]}.{s[1]}{out_name}"
            outfile_list.append(s)
        return outfile_list

    def read_bed(self):
        for line in open(self.bed):
            fields = line.rstrip().split("\t")
            self.ref_numt[fields[0]].addi(int(fields[1]), int(fields[2]), (fields[3], int(fields[2]) - int(fields[1])))

    # def get_bam(self):
    #     return os.listdir(self.bam_folder)

    def numt_Overlap(self, chr:str, start, end):
        overlap = 0
        overlapid = "Mapped"
        overlapcontained = 0
        overlap_count = 0
        nlen = -1
        refoverlap = self.ref_numt[chr].search(start,end)
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
    
    def robust_get_tag(self, read, tag_name, default_value):
        try:  
            return read.get_tag(tag_name)
        except KeyError: 
            return default_value

    def generate_numt_tag(self):
        # if self.inputbam.endswith(".bam"):
        #     samfile = pysam.AlignmentFile(self.inputbam, 'rb')
        # elif self.inputbam.endswith(".sam"):
        #     samfile = pysam.AlignmentFile(self.inputbam, 'r')
        # elif self.inputbam.endswith(".cram"):
        #     samfile = pysam.AlignmentFile(self.inputbam, 'rc')
        samfile = pysam.AlignmentFile(self.inputbam, threads = self.threads, reference_filename = self.reference)
        outfile = pysam.AlignmentFile(self.outputbam, 'wb', header = samfile.header)
        reads_retained = 0
        counts = 0
        for read in samfile:
            if (not(read.is_unmapped)):
                self_info = self.numt_Overlap(read.reference_name, read.reference_start, read.reference_end)
            else:
                self_info = ("0", "0", "Unmapped", "-1" , "0")
            if (not( read.mate_is_unmapped)):
#                x = simplesam.Sam(cigar = read.get_tag("MC"))
                mctag = self.robust_get_tag(read, "MC", "NoMC")
                if mctag != "NoMC":
                    x = simplesam.Sam(cigar = mctag)
                    mate_info = self.numt_Overlap(read.next_reference_name, read.next_reference_start, read.next_reference_start + len(x))
                else:
                    print("WARN: NoMC =>", read)
            else:
                mate_info = ("0", "0", "Unmapped", "-1" , "0")
            if self_info[0] == "1" or mate_info[0] == "1":
                reads_retained += 1
                read.set_tag("Ni", ','.join(self_info + mate_info), value_type = "Z", replace = False)
                outfile.write(read)
            counts += 1
            if counts % 10000 == 0:
                print("Total number of reads processed for the file:", counts)

        print("Total number of reads that overlap numts/MT:",reads_retained)
        print("Total number of reads processed for the file:", counts)
        outfile.close()  

    def step3_for_all_files(self):
        # p = "/mnt/c/Users/simar/JCSMR_numtS_Project/ANU_Numts/nmt/Step3.sh"
        cwd = os.getcwd()
        p = os.path.join(cwd, "Step3.sh")
        for idx, file_name in enumerate(self.outfiles):
            in_path = os.path.join(hardcoded_path.out_dir, file_name)
            out_path = os.path.join(hardcoded_path.step3_output_dir, self.step3_out[idx])
            print(in_path,out_path)
            cmd = f" -o={out_path} -i={in_path}"
            cmd = p+cmd
            print(cmd)
            s = subprocess.run([cmd], shell = True)
            # s = subprocess.Popen([cmd], shell = True)
            # out, err = s.communicate()

